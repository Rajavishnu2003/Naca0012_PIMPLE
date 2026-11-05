// Main solver: sets up the case, runs the time/iteration loop and writes the output.
// The SIMPLE, PISO and PIMPLE steps themselves are in simple_solver, piso_solver and pimple_solver.
#include "cfd/coupling_solver.h"
#include "cfd/output_writer.h"
#include <filesystem>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cstdio>
namespace cfd {
    namespace {
        Config validated(Config c) {
            c.validate();
            return c;
        }
        // airfoil uses the domain from the case file, the cavity is the unit square
        StaggeredGrid makeGrid(const Config& c) {
            if(c.case_type=="external")return StaggeredGrid(c.nx,c.ny,c.x_min,c.y_min,c.x_max-c.x_min,c.y_max-c.y_min);
            return StaggeredGrid(c.nx,c.ny);
        }
    }
    CouplingSolver::CouplingSolver(Config c,bool output):c_(validated(std::move(c))),g_(makeGrid(c_)),body_(buildBodyMask(g_,c_)),bc_(BoundaryConditions::forCase(c_)),
    f_(g_),momentum_(g_,body_),residual_workspace_(g_,body_),pressure_(g_,body_),output_(output) {
        c_.validate();
        // external flow starts from the uniform free stream (the body faces are set to zero just below)
        if(c_.case_type=="external")std::fill(f_.u.begin(),f_.u.end(),bc_.inlet_speed);
        bc_.enforceNormalVelocities(g_,body_,f_);
    }
    // momentum predictor: build and solve u and v with the current pressure
    void CouplingSolver::momentumPredictor(const FlowFields& old_time,double dt,double alpha) {
        momentum_.assemble(f_,old_time,bc_,c_,dt,alpha);
        auto u=packU(g_,f_),v=packV(g_,f_);
        auto ru=solveMomentum(momentum_.u.A,momentum_.u.b,u,c_);
        auto rv=solveMomentum(momentum_.v.A,momentum_.v.b,v,c_);
        unpackU(g_,u,f_);
        unpackV(g_,v,f_);
        // put back the known boundary values and fix the outlet so mass in = mass out
        bc_.enforceNormalVelocities(g_,body_,f_);
        bc_.updateOutflow(g_,f_);
        summary_.momentum_iterations+=ru.iterations+rv.iterations;
        summary_.momentum_solves+=2;
        summary_.momentum_solve_seconds+=ru.seconds+rv.seconds;
    }
    // one pressure correction: solve for p' from the mass imbalance, then correct u, v and p
    void CouplingSolver::pressureCorrect(PressureLinearSolver& solver,double alpha_p) {
        pressure_.assembleRHS(f_);
        std::vector<double> pp(g_.np(),0);
        const auto r=solver.solve(pressure_.b,pp);
        summary_.pressure_iterations+=r.iterations;
        ++summary_.pressure_solves;
        summary_.pressure_solve_seconds+=r.seconds;
        correctPressureVelocity(g_,momentum_,pp,alpha_p,f_);
        bc_.enforceNormalVelocities(g_,body_,f_);
    }
    RunSummary CouplingSolver::run() {
        const auto start=std::chrono::steady_clock::now();
        const bool external=c_.case_type=="external";
        std::ofstream history,force_history;
        // list of VTK frames written so far (file name + time), for the .series file
        std::vector<std::pair<std::string,double>> vtk_frames;
        const bool snapshots=output_&&c_.vtk_every>0;
        auto snapshot=[&](int step) {
            char name[32];
            std::snprintf(name,sizeof name,"vtk/flow_%06d.vtk",step);
            writeVTKSnapshot(g_,body_,f_,c_.output+"/"+name,summary_.time);
            vtk_frames.emplace_back(name,summary_.time);
            writeVTKSeries(c_.output+"/flow.vtk.series",vtk_frames);  // rewritten each frame so you can open a run that is still going
        };
        // open the output files: input copy, residual history, force history
        if(output_) {
            std::filesystem::create_directories(c_.output);
            writeConfig(c_,c_.output+"/input.cfg");
            history.open(c_.output+"/residuals.csv");
            if(!history)throw std::runtime_error("Cannot open convergence history.");
            if(external) {
                force_history.open(c_.output+"/forces.csv");
                if(!force_history)throw std::runtime_error("Cannot open force history.");
                force_history<<"step,time,cd,cl,cd_pressure,cd_viscous\n"<<std::setprecision(17);
            }
            if(snapshots) {
                std::filesystem::create_directories(c_.output+"/vtk");
                if(c_.vtk_start_time<=0)snapshot(0);  // initial condition
            }
            history<<"step,time,u_momentum_inf,v_momentum_inf,divergence_inf,velocity_change,courant,cds_faces,upwind_faces,max_abs_peclet,momentum_iterations,pressure_iterations\n"<<std::setprecision(17);
        }
        // main loop: one SIMPLE iteration or one PISO/PIMPLE time step per pass
        for(int step=1;step<=c_.max_steps;++step) {
            const FlowFields previous=f_;
            // SIMPLE is steady so dt = 0, and when running to end_time the last step is cut to land on it
            double dt=c_.algorithm=="simple"?0.0:c_.dt;
            if(!c_.stop_at_steady)dt=std::min(dt,c_.end_time-summary_.time);
            if(c_.algorithm=="simple")simpleIteration(previous);
            else if(c_.algorithm=="piso")pisoTimeStep(previous,dt);
            else pimpleTimeStep(previous,dt);
            summary_.time+=dt;
            summary_.steps=step;
            summary_.residuals=evaluateResiduals(g_,f_,previous,bc_,c_,residual_workspace_,dt);
            const auto& r=summary_.residuals;
            // lift and drag on the airfoil, written every step to forces.csv
            if(external) {
                summary_.forces=computeBodyForces(g_,body_,f_,c_.reynolds);
                const auto& F=summary_.forces;
                if(force_history)force_history<<step<<','<<summary_.time<<','<<F.cd<<','<<F.cl<<','<<2*F.fx_pressure<<','<<2*F.fx_viscous<<'\n';
            }
            if(history)history<<step<<','<<summary_.time<<','<<r.u_momentum<<','<<r.v_momentum<<','<<r.divergence<<','<<r.velocity_change<<','<<r.courant<<','<<r.schemes.cds_faces<<','<<r.schemes.upwind_faces<<','<<r.schemes.max_abs_peclet<<','<<summary_.momentum_iterations<<','<<summary_.pressure_iterations<<'\n';
            // never call it converged in the first two steps
            if(step>2&&steadyConverged(r,c_,dt))summary_.converged=true;
            if(snapshots&&step%c_.vtk_every==0&&summary_.time>=c_.vtk_start_time-1e-9)snapshot(step);
            // progress line in the terminal / log file
            if(output_&&(step==1||step%c_.log_every==0||summary_.converged)) {
                std::cout<<"step="<<step<<" t="<<summary_.time<<" mom="<<std::max(r.u_momentum,r.v_momentum)<<" div="<<r.divergence<<" change="<<r.velocity_change<<" upwind="<<r.schemes.upwind_faces;
                if(external)std::cout<<" Cd="<<summary_.forces.cd<<" Cl="<<summary_.forces.cl;
                std::cout<<'\n';
            }
            // stop when steady (stop_at_steady = true) or when end_time is reached
            if((c_.stop_at_steady&&summary_.converged)||(!c_.stop_at_steady&&summary_.time>=c_.end_time-1e-12)) {
                summary_.completed=true;
                break;
            }
            if(!c_.stop_at_steady)summary_.converged=false;
        }
        summary_.wall_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        // final fields, summary.json and the final flow.vtk
        if(output_)writeResults(g_,body_,f_,c_,summary_,momentum_,pressure_);
        return summary_;
    }
}
