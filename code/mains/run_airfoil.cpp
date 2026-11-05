// main program for the airfoil / external flow cases
// usage: ./airfoil case.cfg, it reads the case, runs the solver and prints Cd and Cl
#include "cfd/coupling_solver.h"
#include <iostream>
#include <exception>
int main(int argc,char** argv) {
    if(argc!=2) {
        std::cerr<<"Usage: airfoil path/to/case.cfg   (case_type = external)\n";
        return 1;
    }
    try {
        // read the case file, this program only does external flow
        auto c=cfd::readConfig(argv[1]);
        if(c.case_type!="external")throw std::invalid_argument("This program needs case_type = external; use ./cavity for the cavity.");
        std::cout<<"CSR staircase external flow: body="<<c.body<<(c.body=="naca"?" NACA "+c.naca:"")<<", AoA="<<c.aoa_deg
        <<" deg, "<<c.algorithm<<", "<<c.nx<<"x"<<c.ny<<", Re="<<c.reynolds<<", convection="<<c.convection<<'\n';
        // build grid, body mask and fields, then march in time
        cfd::CouplingSolver solver(c);
        std::cout<<"solid cells="<<solver.body().solid_cells<<" of "<<solver.grid().np()<<", h="<<solver.grid().dx<<" x "<<solver.grid().dy<<'\n';
        auto s=solver.run();
        std::cout<<"completed="<<s.completed<<" steady_converged="<<s.converged<<" steps="<<s.steps<<" wall_seconds="<<s.wall_seconds<<'\n';
        // Cd = 2*Fx because 0.5*rho*U^2*c = 0.5
        std::cout<<"Cd="<<s.forces.cd<<" (pressure "<<2*s.forces.fx_pressure<<", viscous "<<2*s.forces.fx_viscous<<")  Cl="<<s.forces.cl<<'\n';
        return s.completed?0:2;  // 2 means it hit max_steps before finishing
    }
    catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        return 1;
    }
}
