// geometric multigrid for the pressure equation (only for the cavity, not with a body)
#pragma once
#include "csr_matrix.h"
#include "cfd/config.h"
namespace cfd {
    struct PressureMGLevel {
        int nx=0,ny=0;
        CSRMatrix A,P,R;
    };
    // cell centred multigrid with V, W or F cycles, adapted from the V/W/F cycle code we were given.
    // coarse matrix A_coarse = P^T A_fine P, everything stays in CSR. it is its own version, not the old heat equation one
    class PressureMultigrid {
        public:
        PressureMultigrid(const CSRMatrix& A,int nx,int ny,int pre,int post);
        void applyCycle(std::vector<double>& x,const std::vector<double>& b,char kind='V') const;
        void applyFMG(std::vector<double>& x,const std::vector<double>& b,char kind='V') const;
        const std::vector<PressureMGLevel>& levels() const {
            return levels_;
        }
        private:
        std::vector<PressureMGLevel> levels_;
        int pre_,post_;
        void cycle(int level,std::vector<double>& x,const std::vector<double>& b,char kind) const;
    };
    CSRMatrix cellCentredProlongation(int nfx,int nfy,int ncx,int ncy,bool zero_reference_row);
}
