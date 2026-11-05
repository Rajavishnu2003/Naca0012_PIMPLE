// multigrid for the pressure equation, coarse grids are half the size each time
#include "cfd/pressure_multigrid.h"
#include "cfd/csr_tools.h"
#include "conjugate_gradient.h"
#include "ic0_preconditioner.h"
#include "linear_algebra.h"
#include "smoother.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
namespace cfd {
    // prolongation P: bilinear interpolation from coarse cell centres to fine cell centres
    CSRMatrix cellCentredProlongation(int fx,int fy,int cx,int cy,bool pin) {
        CSRMatrix P(fx*fy,cx*cy);
        auto weights=[](int i,int nf,int nc) {
            double s=(i+0.5)*nc/nf-0.5;
            int lo=static_cast<int>(std::floor(s));
            double a=s-lo;
            std::map<int,double> w;
            w[std::clamp(lo,0,nc-1)]+=1-a;
            w[std::clamp(lo+1,0,nc-1)]+=a;
            return w;
        };
        for(int j=0;j<fy;++j)for(int i=0;i<fx;++i) {
            int r=j*fx+i;
            P.rowPtr[r]=P.nnz();
            if(pin&&r==0)continue;  // the pinned reference cell gets no correction
            auto wx=weights(i,fx,cx),wy=weights(j,fy,cy);
            std::map<int,double> row;
            for(auto y:wy)for(auto x:wx)row[y.first*cx+x.first]+=y.second*x.second;
            for(auto kv:row)if(kv.second!=0) {
                P.colIndex.push_back(kv.first);
                P.values.push_back(kv.second);
            }
        }
        P.rowPtr[P.rows]=P.nnz();
        return P;
    }
    namespace {
        // Gauss Seidel going from the last row to the first, used after the coarse correction
        void backwardGS(const CSRMatrix& A,std::vector<double>& x,const std::vector<double>& b,int sweeps) {
            for(int t=0;t<sweeps;++t)for(int i=A.rows-1;i>=0;--i) {
                double sum=b[i];
                for(int k=A.rowPtr[i];k<A.rowPtr[i+1];++k)if(A.colIndex[k]!=i)sum-=A.values[k]*x[A.colIndex[k]];
                x[i]=sum/A.values[A.diagonalIndex[i]];
            }
        }
    }
    PressureMultigrid::PressureMultigrid(const CSRMatrix& A,int nx,int ny,int pre,int post):pre_(pre),post_(post) {
        if(nx*ny!=A.rows||A.rows!=A.cols||pre<1||post<1)throw std::invalid_argument("Invalid MG setup.");
        levels_.push_back( {
            nx,ny,A,{},{}
        });
        // keep halving the grid until it is small, restriction R is just P transposed
        while(nx*ny>16&&nx>=4&&ny>=4) {
            const int cx=nx/2,cy=ny/2;
            const std::size_t l=levels_.size()-1;
            auto P=cellCentredProlongation(nx,ny,cx,cy,l==0);
            auto R=transposeCSR(P);
            auto coarse=productCSR(R,productCSR(levels_[l].A,P));
            validateCSR(coarse);
            levels_[l].P=std::move(P);
            levels_[l].R=std::move(R);
            levels_.push_back( {
                cx,cy,std::move(coarse),{}, {}
            });
            nx=cx;
            ny=cy;
        }
    }
    void PressureMultigrid::cycle(int level,std::vector<double>& x,const std::vector<double>& b,char kind) const {
        const auto& L=levels_.at(level);
        if(level==static_cast<int>(levels_.size())-1) {
            // on the coarsest grid just solve it with CG, it is small and symmetric
            conjugate_gradient_solver coarse(1e-12,2000);
            coarse.solve(L.A,b,x);
            const double bn=vectorNorm(b),rn=vectorNorm(computeResidual(L.A,x,b));
            if(rn>1e-10*std::max(bn,1e-30))throw std::runtime_error("MG coarse solve did not converge.");
            return;
        }
        // smooth, restrict the residual, solve the error on the coarse grid, then add it back and smooth again
        gaussSeidel(L.A,x,b,pre_);
        const auto residual=computeResidual(L.A,x,b);
        const auto rhs=multiply(L.R,residual);
        std::vector<double> correction(rhs.size(),0);
        cycle(level+1,correction,rhs,kind);
        // W cycle visits the coarse level twice, F cycle does one more V
        if(kind=='W')cycle(level+1,correction,rhs,'W');
        else if(kind=='F')cycle(level+1,correction,rhs,'V');
        const auto ef=multiply(L.P,correction);
        axpy(x,1.0,ef);
        backwardGS(L.A,x,b,post_);
    }
    void PressureMultigrid::applyCycle(std::vector<double>& x,const std::vector<double>& b,char kind) const {
        if(kind!='V'&&kind!='W'&&kind!='F')throw std::invalid_argument("Invalid MG cycle.");
        cycle(0,x,b,kind);
    }
    // full multigrid: start on the coarsest grid and work up, gives a good first guess
    void PressureMultigrid::applyFMG(std::vector<double>& x,const std::vector<double>& b,char kind) const {
        std::vector<std::vector<double>> rhs(levels_.size());
        rhs[0]=b;
        for(std::size_t l=0;l+1<levels_.size();++l)rhs[l+1]=multiply(levels_[l].R,rhs[l]);
        std::vector<double> current(rhs.back().size(),0);
        cycle(static_cast<int>(levels_.size())-1,current,rhs.back(),kind);
        for(int l=static_cast<int>(levels_.size())-2;l>=0;--l) {
            current=multiply(levels_[l].P,current);
            cycle(l,current,rhs[l],kind);
        }
        x=std::move(current);
    }
}
