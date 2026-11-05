// makes the solid mask: which cells and faces are inside the airfoil
#include "cfd/body_geometry.h"
#include <cmath>
#include <stdexcept>
namespace cfd {
    BodyMask::BodyMask(const StaggeredGrid& g):p(g.np(),0),u((g.nx+1)*g.ny,0),v(g.nx*(g.ny+1),0){}
    bool insideNaca4(const std::string& code,double x,double y) {
        if(x<0||x>1) return false;
        const double m=(code[0]-'0')/100.0;   // maximum camber
        const double pc=(code[1]-'0')/10.0;   // position of maximum camber
        const double t=std::stoi(code.substr(2))/100.0;  // maximum thickness
        // half thickness, we use -0.1036 instead of -0.1015 so the trailing edge is closed
        const double yt=5*t*(0.2969*std::sqrt(x)-0.1260*x-0.3516*x*x+0.2843*x*x*x-0.1036*x*x*x*x);
        double yc=0;
        if(m>0) yc=x<pc?m/(pc*pc)*(2*pc*x-x*x):m/((1-pc)*(1-pc))*((1-2*pc)+2*pc*x-x*x);
        // thickness taken straight up and down from the camber line (exact for 0012 since it has no camber)
        return std::abs(y-yc)<=yt;
    }
    bool insideBody(const Config& c,double x,double y) {
        if(c.body=="cylinder") return x*x+y*y<=0.25;  // diameter 1, centre at origin
        if(c.body=="naca") {
            // positive AoA means nose up. instead of rotating the airfoil we rotate
            // the point the other way about the quarter chord (0.25,0) and check it in the airfoil frame
            const double a=c.aoa_deg*std::acos(-1.0)/180.0;
            const double rx=x-0.25,ry=y;
            const double xa=0.25+rx*std::cos(a)-ry*std::sin(a);
            const double ya=rx*std::sin(a)+ry*std::cos(a);
            return insideNaca4(c.naca,xa,ya);
        }
        return false;
    }
    BodyMask buildBodyMask(const StaggeredGrid& g,const Config& c) {
        BodyMask b(g);
        if(c.case_type!="external"||c.body=="none") return b;
        // mark solid cells by checking each cell centre
        for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) b.p[g.p(i,j)]=insideBody(c,g.xc(i),g.yc(j));
        if(b.p[g.p(0,0)]) throw std::runtime_error("Body covers the pressure reference cell; enlarge the domain.");
        // flood fill starting from the corner cell. any fluid cell we can not reach
        // (like a pocket stuck near the thin trailing edge) makes the pressure system singular, so we make it solid
        std::vector<char> reached(g.np(),0);
        std::vector<int> stack{g.p(0,0)};
        reached[g.p(0,0)]=1;
        while(!stack.empty()) {
            const int k=stack.back();
            stack.pop_back();
            const int i=k%g.nx,j=k/g.nx;
            const int ni[4]={i+1,i-1,i,i},nj[4]={j,j,j+1,j-1};
            for(int n=0;n<4;++n) {
                if(ni[n]<0||ni[n]>=g.nx||nj[n]<0||nj[n]>=g.ny) continue;
                const int q=g.p(ni[n],nj[n]);
                if(!reached[q]&&!b.p[q]) {
                    reached[q]=1;
                    stack.push_back(q);
                }
            }
        }
        for(int k=0;k<g.np();++k) if(!reached[k]) b.p[k]=1;
        for(int k=0;k<g.np();++k) b.solid_cells+=b.p[k];
        if(b.solid_cells==0) throw std::runtime_error("Body is smaller than one cell; refine the grid.");
        // a face touching a solid cell is solid, the velocity there is set to zero
        for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i)
        b.u[g.uf(i,j)]=(i>0&&b.cell(g,i-1,j))||(i<g.nx&&b.cell(g,i,j));
        for(int j=0;j<=g.ny;++j)for(int i=0;i<g.nx;++i)
        b.v[g.vf(i,j)]=(j>0&&b.cell(g,i,j-1))||(j<g.ny&&b.cell(g,i,j));
        return b;
    }
}
