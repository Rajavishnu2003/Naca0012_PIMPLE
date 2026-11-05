// all output files: fields.csv, flow.vtk, face csv files, summary.json and the vtk snapshots
#include "cfd/output_writer.h"
#include "cfd/csr_tools.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <cstring>
namespace cfd {
    namespace {
        std::ofstream open(const std::string& p) {
            std::ofstream f(p);
            if(!f)throw std::runtime_error("Cannot write "+p);
            f<<std::setprecision(17);
            return f;
        }
        // velocity at the cell centre = average of the two faces of the cell
        void cellVelocities(const StaggeredGrid& g,const FlowFields& f,std::vector<double>& uc,std::vector<double>& vc) {
            uc.assign(g.np(),0);
            vc.assign(g.np(),0);
            for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
                uc[g.p(i,j)]=0.5*(f.u[g.uf(i,j)]+f.u[g.uf(i+1,j)]);
                vc[g.p(i,j)]=0.5*(f.v[g.vf(i,j)]+f.v[g.vf(i,j+1)]);
            }
        }
        // vorticity w = dv/dx - du/dy at cell centres. central difference inside, one sided at the domain edges
        std::vector<double> vorticity(const StaggeredGrid& g,const std::vector<double>& uc,const std::vector<double>& vc) {
            std::vector<double> w(g.np(),0);
            for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
                const int il=std::max(i-1,0),ir=std::min(i+1,g.nx-1),jb=std::max(j-1,0),jt=std::min(j+1,g.ny-1);
                const double dvdx=(vc[g.p(ir,j)]-vc[g.p(il,j)])/((ir-il)*g.dx);
                const double dudy=(uc[g.p(i,jt)]-uc[g.p(i,jb)])/((jt-jb)*g.dy);
                w[g.p(i,j)]=dvdx-dudy;
            }
            return w;
        }
        // old style binary vtk wants big endian, so swap the bytes on a normal (little endian) PC
        void putBE(std::ostream& o,float x) {
            unsigned char b[4];
            std::memcpy(b,&x,4);
            const std::uint16_t one=1;
            if(*reinterpret_cast<const unsigned char*>(&one)==1)std::swap(b[0],b[3]),std::swap(b[1],b[2]);
            o.write(reinterpret_cast<const char*>(b),4);
        }
    }
    void writeVTKSnapshot(const StaggeredGrid& g,const BodyMask& body,const FlowFields& f,const std::string& path,double time) {
        std::ofstream vtk(path,std::ios::binary);
        if(!vtk)throw std::runtime_error("Cannot write "+path);
        std::vector<double> uc,vc;
        cellVelocities(g,f,uc,vc);
        const auto w=vorticity(g,uc,vc);
        vtk<<std::setprecision(17)<<"# vtk DataFile Version 3.0\nCSR staggered external flow t="<<time<<"\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS "<<g.nx<<' '<<g.ny<<" 1\nORIGIN "<<g.xc(0)<<' '<<g.yc(0)<<" 0\nSPACING "<<g.dx<<' '<<g.dy<<" 1\n";
        // store the time in the file too, then the fields one after the other
        vtk<<"FIELD FieldData 1\nTIME 1 1 float\n";
        putBE(vtk,float(time));
        vtk<<"\nPOINT_DATA "<<g.np()<<"\nVECTORS velocity float\n";
        for(int k=0;k<g.np();++k)putBE(vtk,float(uc[k])),putBE(vtk,float(vc[k])),putBE(vtk,0.f);
        vtk<<"\nSCALARS pressure float 1\nLOOKUP_TABLE default\n";
        for(double a:f.p)putBE(vtk,float(a));
        vtk<<"\nSCALARS vorticity float 1\nLOOKUP_TABLE default\n";
        for(double a:w)putBE(vtk,float(a));
        vtk<<"\nSCALARS solid float 1\nLOOKUP_TABLE default\n";
        for(char a:body.p)putBE(vtk,float(a));
        vtk<<'\n';
    }
    void writeVTKSeries(const std::string& path,const std::vector<std::pair<std::string,double>>& frames) {
        auto s=open(path);
        s<<"{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n";
        for(std::size_t k=0;k<frames.size();++k)
            s<<"    { \"name\" : \""<<frames[k].first<<"\", \"time\" : "<<frames[k].second<<" }"<<(k+1<frames.size()?",":"")<<'\n';
        s<<"  ]\n}\n";
    }
    void writeResults(const StaggeredGrid& g,const BodyMask& body,const FlowFields& f,const Config& c,const RunSummary& s,const MomentumAssembler& m,const PressureAssembler& p) {
        const bool external=c.case_type=="external";
        // fields.csv and flow.vtk have the same cell centre values
        auto field=open(c.output+"/fields.csv");
        field<<(external?"x,y,u,v,p,solid\n":"x,y,u,v,p\n");
        auto vtk=open(c.output+(external?"/flow.vtk":"/cavity.vtk"));
        vtk<<"# vtk DataFile Version 3.0\nCSR staggered "<<(external?"external flow":"cavity")<<" (cell-centred output)\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "<<g.nx<<' '<<g.ny<<" 1\nORIGIN "<<g.xc(0)<<' '<<g.yc(0)<<" 0\nSPACING "<<g.dx<<' '<<g.dy<<" 1\nPOINT_DATA "<<g.np()<<"\nVECTORS velocity double\n";
        for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            double u=0.5*(f.u[g.uf(i,j)]+f.u[g.uf(i+1,j)]),v=0.5*(f.v[g.vf(i,j)]+f.v[g.vf(i,j+1)]);
            field<<g.xc(i)<<','<<g.yc(j)<<','<<u<<','<<v<<','<<f.p[g.p(i,j)];
            if(external)field<<','<<int(body.cell(g,i,j));
            field<<'\n';
            vtk<<u<<' '<<v<<" 0\n";
        }
        vtk<<"SCALARS pressure double 1\nLOOKUP_TABLE default\n";
        for(double a:f.p)vtk<<a<<'\n';
        {
            std::vector<double> uc,vc;
            cellVelocities(g,f,uc,vc);
            vtk<<"SCALARS vorticity double 1\nLOOKUP_TABLE default\n";
            for(double a:vorticity(g,uc,vc))vtk<<a<<'\n';
        }
        if(external) {
            vtk<<"SCALARS solid int 1\nLOOKUP_TABLE default\n";
            for(char a:body.p)vtk<<int(a)<<'\n';
        }
        // raw face values before averaging
        auto uf=open(c.output+"/u_faces.csv");
        uf<<"i,j,x,y,u\n";
        for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i)uf<<i<<','<<j<<','<<g.xf(i)<<','<<g.yc(j)<<','<<f.u[g.uf(i,j)]<<'\n';
        auto vf=open(c.output+"/v_faces.csv");
        vf<<"i,j,x,y,v\n";
        for(int j=0;j<=g.ny;++j)for(int i=0;i<g.nx;++i)vf<<i<<','<<j<<','<<g.xc(i)<<','<<g.yf(j)<<','<<f.v[g.vf(i,j)]<<'\n';
        // cavity only: centre line profiles to compare with the reference data
        if(!external) {
            auto uc=open(c.output+"/centerline_u.csv");
            uc<<"y,u\n0,0\n";
            int il=g.nx/2;
            double ax=0.5*g.nx-il;
            for(int j=0;j<g.ny;++j)uc<<(j+.5)*g.dy<<','<<(1-ax)*f.u[g.uf(il,j)]+ax*f.u[g.uf(il+1,j)]<<'\n';
            uc<<"1,1\n";
            auto vc=open(c.output+"/centerline_v.csv");
            vc<<"x,v\n0,0\n";
            int jl=g.ny/2;
            double ay=0.5*g.ny-jl;
            for(int i=0;i<g.nx;++i)vc<<(i+.5)*g.dx<<','<<(1-ay)*f.v[g.vf(i,jl)]+ay*f.v[g.vf(i,jl+1)]<<'\n';
            vc<<"1,0\n";
        }
        // summary.json has the final numbers, the python scripts read this
        auto json=open(c.output+"/summary.json");
        json<<std::boolalpha<<"{\n";
        json<<"  \"algorithm\": \""<<c.algorithm<<"\",\n  \"pressure_solver\": \""<<c.pressure_solver<<"\",\n  \"convection\": \""<<c.convection<<"\",\n";
#define J(key,val) json<<"  \""<<key<<"\": "<<(val)<<",\n";
        J("nx",g.nx) J("ny",g.ny) J("Re",c.reynolds) J("dt",c.algorithm=="simple"?0:c.dt)
        J("converged",s.converged) J("completed",s.completed) J("steps",s.steps) J("time",s.time)
        J("wall_seconds",s.wall_seconds) J("momentum_solve_seconds",s.momentum_solve_seconds) J("pressure_solve_seconds",s.pressure_solve_seconds) J("pressure_setup_seconds",s.pressure_setup_seconds)
        J("momentum_iterations",s.momentum_iterations) J("pressure_iterations",s.pressure_iterations) J("pressure_solves",s.pressure_solves) J("momentum_solves",s.momentum_solves)
        J("u_momentum_inf",s.residuals.u_momentum) J("v_momentum_inf",s.residuals.v_momentum) J("divergence_inf",s.residuals.divergence)
        J("velocity_change",s.residuals.velocity_change) J("courant",s.residuals.courant)
        J("cds_faces",s.residuals.schemes.cds_faces) J("upwind_faces",s.residuals.schemes.upwind_faces) J("max_abs_peclet",s.residuals.schemes.max_abs_peclet)
        J("u_unknowns",g.nu()) J("v_unknowns",g.nv()) J("p_unknowns",g.np()) J("u_nnz",m.u.A.nnz()) J("v_nnz",m.v.A.nnz())
        if(external) {
            json<<"  \"body\": \""<<c.body<<"\",\n  \"naca\": \""<<c.naca<<"\",\n";
            J("aoa_deg",c.aoa_deg) J("solid_cells",body.solid_cells)
            J("cd",s.forces.cd) J("cl",s.forces.cl)
            // times 2 because the coefficient is force / 0.5
            J("cd_pressure",2*s.forces.fx_pressure) J("cd_viscous",2*s.forces.fx_viscous)
            J("cl_pressure",2*s.forces.fy_pressure) J("cl_viscous",2*s.forces.fy_viscous)
        }
#undef J
        json<<"  \"p_nnz\": "<<p.A.nnz()<<"\n}\n";
        if(c.dump_matrices) {
            dumpCSR(m.u.A,m.u.b,c.output+"/u_predictor");
            dumpCSR(m.v.A,m.v.b,c.output+"/v_predictor");
            dumpCSR(p.A,p.b,c.output+"/pressure_correction");
        }
    }
}
