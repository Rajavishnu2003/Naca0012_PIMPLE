// Builds the u and v momentum equations on the staggered grid.
// Each face velocity gets one row: a_P u_P = sum(a_nb u_nb) + b, the standard FV form.
#include "cfd/momentum_assembler.h"
#include "cfd/csr_tools.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd {
    // u has (nx-1) x ny unknowns and v has nx x (ny-1), both with a 5 point stencil
    MomentumSystem::MomentumSystem(int nx,int ny):A(fivePointPattern(nx,ny)),b(nx*ny,0),d(nx*ny,0){}
    MomentumAssembler::MomentumAssembler(const StaggeredGrid& g,const BodyMask& body):u(g.nx-1,g.ny),v(g.nx,g.ny-1),g_(g),body_(body){}
    namespace {
        // clear the matrix values but keep the sparsity pattern
        void reset(MomentumSystem& s) {
            zeroValues(s.A);
            std::fill(s.b.begin(),s.b.end(),0);
            s.stats={};
        }
        // add one face to row r: if the neighbour is an unknown it goes in the matrix,
        // if it is a known value (wall, inlet, body) it moves to the right hand side
        void face(MomentumSystem& s,int r,int neighbour,double F,double D,double known,const std::string& scheme) {
            const auto a=faceCoefficients(F,D,scheme);
            s.A.values[s.A.diagonalIndex[r]]+=a.diagonal;
            if(neighbour>=0)addCoefficient(s.A,r,neighbour,a.neighbour);
            else s.b[r]-=a.neighbour*known;
            if(a.cds)++s.stats.cds_faces;
            else ++s.stats.upwind_faces;
            s.stats.max_abs_peclet=std::max(s.stats.max_abs_peclet,std::abs(a.peclet));
        }
        // zero gradient outlet face: the face value equals the cell value, so the outgoing
        // convection only adds F to the diagonal and there is no diffusion
        // backflow (F < 0) would bring in unknown momentum, so we treat it as zero
        void outflowFace(MomentumSystem& s,int r,double F) {
            s.A.values[s.A.diagonalIndex[r]]+=std::max(F,0.0);
        }
        // blocking off: a face inside the body gets the equation 1*u = 0 and no
        // pressure correction coefficient (d = 0), so it stays exactly zero
        void blockOff(MomentumSystem& s,int r) {
            s.A.values[s.A.diagonalIndex[r]]=1.0;
            s.b[r]=0.0;
            s.d[r]=0.0;
        }
        // last part of each row: time term, under-relaxation and the d coefficient
        void finish(MomentumSystem& s,int r,double area,double mass_dt,double old_time,double old_iter,double alpha) {
            auto& diag=s.A.values[s.A.diagonalIndex[r]];
            // implicit Euler time term: (V/dt) u on the left, (V/dt) u_old on the right
            diag+=mass_dt;
            s.b[r]+=mass_dt*old_time;
            const double unrelaxed=diag;
            if(!(unrelaxed>0))throw std::runtime_error("Nonpositive momentum diagonal. Check CDS/cell Peclet and relaxation.");
            // Patankar relaxation: a_P/alpha on the left and (1-alpha)/alpha * a_P * u_prev on the right
            diag=unrelaxed/alpha;
            s.b[r]+=(1.0-alpha)*diag*old_iter;
            // d = face area / a_P, this links the velocity correction to p' later
            s.d[r]=area/diag;
            // alpha comes in only once, through the relaxed diagonal
        }
    }
    // neighbour rules (written for u, v is the same with x and y swapped):
    // in the flow direction (E/W) a body face or domain edge is a known normal velocity
    // at the full distance dx, so conductance D with the known value
    // across the flow (N/S) a no-slip wall sits on the cell face, half a cell away,
    // so conductance 2D with the known wall speed (same as the cavity lid)
    void MomentumAssembler::assemble(const FlowFields& f,const FlowFields& old,
    const BoundaryConditions& bc,const Config& c,double dt,double alpha) {
        if(alpha<=0||alpha>1||dt<0)throw std::invalid_argument("Invalid assembly controls.");
        const auto& g=g_;
        const auto& body=body_;
        const auto& scheme=c.convection;
        // rho = 1 and U = 1, so the viscosity is just 1/Re
        const double mu=1.0/c.reynolds;
        const double Dx=mu*g.dy/g.dx,Dy=mu*g.dx/g.dy;
        const double mass_dt=dt>0?g.volume()/dt:0;
        const bool outlet=bc.east==BoundaryType::Outlet;
        reset(u);
        reset(v);
        // u momentum, one control volume around every inner u face
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i) {
            int r=g.u(i,j);
            if(body.uFace(g,i,j)) {
                blockOff(u,r);
                continue;
            }
            // mass fluxes through the 4 sides of the u control volume (outward positive)
            double E=0.5*(f.u[g.uf(i,j)]+f.u[g.uf(i+1,j)])*g.dy;
            double W=-0.5*(f.u[g.uf(i-1,j)]+f.u[g.uf(i,j)])*g.dy;
            double N=0.5*(f.v[g.vf(i-1,j+1)]+f.v[g.vf(i,j+1)])*g.dx;
            double S=-0.5*(f.v[g.vf(i-1,j)]+f.v[g.vf(i,j)])*g.dx;
            // east / west
            if(i<g.nx-1&&!body.uFace(g,i+1,j))face(u,r,g.u(i+1,j),E,Dx,0,scheme);
            else if(i==g.nx-1&&outlet)outflowFace(u,r,E);
            else face(u,r,-1,E,Dx,f.u[g.uf(i+1,j)],scheme);
            if(i>1&&!body.uFace(g,i-1,j))face(u,r,g.u(i-1,j),W,Dx,0,scheme);
            else face(u,r,-1,W,Dx,f.u[g.uf(i-1,j)],scheme);
            // north / south (a slip boundary adds nothing here, so no shear on it)
            if(j<g.ny-1&&!body.uFace(g,i,j+1))face(u,r,g.u(i,j+1),N,Dy,0,scheme);
            else if(j<g.ny-1)face(u,r,-1,N,2*Dy,0,scheme);
            else if(bc.north!=BoundaryType::Slip)face(u,r,-1,N,2*Dy,bc.lid_speed,scheme);
            if(j>0&&!body.uFace(g,i,j-1))face(u,r,g.u(i,j-1),S,Dy,0,scheme);
            else if(j>0||bc.south!=BoundaryType::Slip)face(u,r,-1,S,2*Dy,0,scheme);
            // pressure force: (p_west - p_east) times the face area
            u.b[r]+=(f.p[g.p(i-1,j)]-f.p[g.p(i,j)])*g.dy;
            finish(u,r,g.dy,mass_dt,old.u[g.uf(i,j)],f.u[g.uf(i,j)],alpha);
        }
        // v momentum, same idea around every inner v face
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            int r=g.v(i,j);
            if(body.vFace(g,i,j)) {
                blockOff(v,r);
                continue;
            }
            double E=0.5*(f.u[g.uf(i+1,j-1)]+f.u[g.uf(i+1,j)])*g.dy;
            double W=-0.5*(f.u[g.uf(i,j-1)]+f.u[g.uf(i,j)])*g.dy;
            double N=0.5*(f.v[g.vf(i,j)]+f.v[g.vf(i,j+1)])*g.dx;
            double S=-0.5*(f.v[g.vf(i,j-1)]+f.v[g.vf(i,j)])*g.dx;
            // east / west (across the flow direction of v), the inlet has v = 0 like a wall
            if(i<g.nx-1&&!body.vFace(g,i+1,j))face(v,r,g.v(i+1,j),E,Dx,0,scheme);
            else if(i==g.nx-1&&outlet)outflowFace(v,r,E);
            else face(v,r,-1,E,2*Dx,0,scheme);
            if(i>0&&!body.vFace(g,i-1,j))face(v,r,g.v(i-1,j),W,Dx,0,scheme);
            else face(v,r,-1,W,2*Dx,0,scheme);
            // north / south
            if(j<g.ny-1&&!body.vFace(g,i,j+1))face(v,r,g.v(i,j+1),N,Dy,0,scheme);
            else face(v,r,-1,N,Dy,f.v[g.vf(i,j+1)],scheme);
            if(j>1&&!body.vFace(g,i,j-1))face(v,r,g.v(i,j-1),S,Dy,0,scheme);
            else face(v,r,-1,S,Dy,f.v[g.vf(i,j-1)],scheme);
            v.b[r]+=(f.p[g.p(i,j-1)]-f.p[g.p(i,j)])*g.dx;
            finish(v,r,g.dx,mass_dt,old.v[g.vf(i,j)],f.v[g.vf(i,j)],alpha);
        }
    }
    // PISO 2nd (and later) corrector: keep the matrix, but recompute
    // u = (b + pressure change - sum(a_nb u_nb)) / a_P with the latest neighbour velocities
    void updateMomentumH(const StaggeredGrid& g,const MomentumAssembler& m,const std::vector<double>& p0,FlowFields& f) {
        auto ux=packU(g,f),vx=packV(g,f),newu=ux,newv=vx;
        // sum of the off-diagonal terms of row r, i.e. sum(a_nb u_nb)
        auto off=[](const CSRMatrix& A,const std::vector<double>& x,int r) {
            double s=0;
            for(int k=A.rowPtr[r];k<A.rowPtr[r+1];++k)if(A.colIndex[k]!=r)s+=A.values[k]*x[A.colIndex[k]];
            return s;
        };
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i) {
            int r=g.u(i,j),l=g.p(i-1,j),e=g.p(i,j);
            if(m.u.d[r]==0)continue;  // body face stays zero
            // b already has the pressure the matrix was built with (p0), so only add the change since then
            double dp=(f.p[l]-f.p[e])-(p0[l]-p0[e]);
            newu[r]=(m.u.b[r]+g.dy*dp-off(m.u.A,ux,r))/m.u.A.values[m.u.A.diagonalIndex[r]];
        }
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            int r=g.v(i,j),s=g.p(i,j-1),n=g.p(i,j);
            if(m.v.d[r]==0)continue;  // body face stays zero
            double dp=(f.p[s]-f.p[n])-(p0[s]-p0[n]);
            newv[r]=(m.v.b[r]+g.dx*dp-off(m.v.A,vx,r))/m.v.A.values[m.v.A.diagonalIndex[r]];
        }
        unpackU(g,newu,f);
        unpackV(g,newv,f);
    }
}
