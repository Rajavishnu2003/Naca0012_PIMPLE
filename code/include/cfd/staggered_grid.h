// staggered grid: p at cell centres, u on vertical faces, v on horizontal faces
#pragma once
namespace cfd {
    struct StaggeredGrid {
        int nx,ny;
        double x0,y0;  // lower-left corner of the domain
        double dx,dy;
        // default domain is the unit square (that is for the cavity)
        StaggeredGrid(int nx_,int ny_,double x0_=0.0,double y0_=0.0,double lx=1.0,double ly=1.0);
        // number of p cells, u unknowns and v unknowns
        int np() const {
            return nx*ny;
        }
        int nu() const {
            return (nx-1)*ny;
        }
        int nv() const {
            return nx*(ny-1);
        }
        // pressure index, i = 0..nx-1 and j = 0..ny-1, the cell centre is at (xc(i), yc(j))
        int p(int i,int j) const {
            return j*nx+i;
        }
        // index in the full u and v arrays, these also hold the fixed velocities on the boundary
        int uf(int i,int j) const {
            return j*(nx+1)+i;
        }
        int vf(int i,int j) const {
            return j*nx+i;
        }
        // index in the matrix (unknowns only), the boundary normal velocities are left out
        int u(int i,int j) const {
            return j*(nx-1)+(i-1);
        }
        // i=1..nx-1, j=0..ny-1
        int v(int i,int j) const {
            return (j-1)*nx+i;
        }
        // i=0..nx-1, j=1..ny-1
        double volume() const {
            return dx*dy;
        }
        // coordinates of cell centres (xc, yc) and faces (xf, yf)
        double xc(int i) const {
            return x0+(i+0.5)*dx;
        }
        double yc(int j) const {
            return y0+(j+0.5)*dy;
        }
        double xf(int i) const {
            return x0+i*dx;
        }
        double yf(int j) const {
            return y0+j*dy;
        }
    };
}
