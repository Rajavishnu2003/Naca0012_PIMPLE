// Cd and Cl from the pressure and shear on the staircase surface
#include "cfd/body_forces.h"
namespace cfd {
    BodyForces computeBodyForces(const StaggeredGrid& g,const BodyMask& body,const FlowFields& f,double reynolds) {
        BodyForces F;
        const double mu=1.0/reynolds;
        // loop over fluid cells, only the ones next to the body add something
        for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            if(body.cell(g,i,j)) continue;
            const double p=f.p[g.p(i,j)];
            // velocity at the cell centre, which is half a cell away from the wall
            const double uc=0.5*(f.u[g.uf(i,j)]+f.u[g.uf(i+1,j)]);
            const double vc=0.5*(f.v[g.vf(i,j)]+f.v[g.vf(i,j+1)]);
            // wall shear tau = mu * tangential velocity / (half cell), it pulls the body along with the flow
            const double tau_x=mu*uc/(0.5*g.dy),tau_y=mu*vc/(0.5*g.dx);
            // pressure pushes on the body from the fluid side (-p*n with n the outward body normal)
            if(i+1<g.nx&&body.cell(g,i+1,j)) {F.fx_pressure+=p*g.dy; F.fy_viscous+=tau_y*g.dy;}  // body on the east
            if(i>0&&body.cell(g,i-1,j))      {F.fx_pressure-=p*g.dy; F.fy_viscous+=tau_y*g.dy;}  // body on the west
            if(j+1<g.ny&&body.cell(g,i,j+1)) {F.fy_pressure+=p*g.dx; F.fx_viscous+=tau_x*g.dx;}  // body to the north
            if(j>0&&body.cell(g,i,j-1))      {F.fy_pressure-=p*g.dx; F.fx_viscous+=tau_x*g.dx;}  // body to the south
        }
        const double q=0.5;  // 0.5*rho*U^2*c
        // force = pressure part + viscous part, then Cd = 2*Fx and Cl = 2*Fy
        F.cd=(F.fx_pressure+F.fx_viscous)/q;
        F.cl=(F.fy_pressure+F.fy_viscous)/q;
        return F;
    }
}
