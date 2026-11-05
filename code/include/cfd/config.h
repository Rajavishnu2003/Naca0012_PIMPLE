// all the case settings that can be set in the .cfg file
#pragma once
#include <string>
namespace cfd {
    // everything is non dimensional: rho = 1, speed = 1, length = 1, so nu = 1/Re.
    // cavity: unit square with lid speed 1.
    // external: uniform inflow U = 1 from the left past a body with chord (or diameter) 1.
    struct Config {
        std::string case_type="cavity";
        // cavity, external
        int nx=16, ny=16;
        double reynolds=100.0;
        std::string algorithm="simple";
        // simple, piso, pimple
        std::string pressure_solver="cg_ic0";
        // cg_ic0, cg, sd_ic0, bicgstab, multigrid
        std::string convection="peclet";
        // peclet (hybrid), cds, upwind. the last two are only for comparison
        double alpha_u=0.7, alpha_p=0.3;  // under relaxation for velocity and pressure
        double dt=0.01;
        int pressure_correctors=2, outer_correctors=3;  // PISO / PIMPLE loops per time step
        int max_steps=20000, max_linear_iterations=3000;
        double momentum_tolerance=1e-6, continuity_tolerance=1e-8, change_tolerance=1e-8;
        double linear_rtol=1e-10, linear_atol=1e-13;
        // tolerance for the pressure correction relative to its rhs (the mass error), 0 means use linear_rtol.
        // momentum always uses linear_rtol. its rhs has the big relaxation term in it, so a loose
        // tolerance there would skip the solve and SIMPLE would get stuck
        double pressure_rtol=0.0;
        bool stop_at_steady=true;  // true: stop when steady, false: run until end_time
        double end_time=1.0;
        std::string mg_cycle="V";
        int mg_pre=2, mg_post=2;
        bool mg_fmg_start=false;
        int log_every=100;
        // write a binary vtk snapshot every vtk_every steps into <output>/vtk/ (0 = only the final field).
        // <output>/flow.vtk.series lists the frames and times, so ParaView opens them as one animation
        int vtk_every=0;
        // snapshots start when t >= vtk_start_time, so we can skip the start up part for the videos
        double vtk_start_time=0.0;
        std::string output="results/cavity";
        bool dump_matrices=false;
        // only for external flow. staircase body: a cell is solid if its centre is inside
        std::string body="naca";
        // naca (leading edge at origin), cylinder (centred at origin), none (empty channel)
        std::string naca="0012";
        // NACA 4-digit code
        double aoa_deg=0.0;
        // angle of attack in degrees, the airfoil is turned nose up about its quarter chord
        double x_min=-2.0, x_max=4.0, y_min=-2.0, y_max=2.0;  // domain size in chords
        void validate() const;
    };
    Config readConfig(const std::string& path);
    void writeConfig(const Config& c, const std::string& path);
}
// namespace cfd
