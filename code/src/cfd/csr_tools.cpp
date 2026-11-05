// helper functions for the CSR matrices used by the momentum and pressure equations
#include "cfd/csr_tools.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <map>
#include <fstream>
#include <iomanip>
namespace cfd {
    CSRMatrix fivePointPattern(int nx,int ny) {
        if(nx<1||ny<1) throw std::invalid_argument("Invalid stencil dimensions.");
        CSRMatrix A(nx*ny,nx*ny);
        auto push=[&](int r,int c) {
            if(r==c)A.diagonalIndex[r]=A.nnz();
            A.colIndex.push_back(c);
            A.values.push_back(0);
        };
        // entries are pushed in column order: south, west, P, east, north
        for(int j=0;j<ny;++j) for(int i=0;i<nx;++i) {
            int r=j*nx+i;
            A.rowPtr[r]=A.nnz();
            if(j>0) push(r,r-nx);
            if(i>0) push(r,r-1);
            push(r,r);
            if(i<nx-1) push(r,r+1);
            if(j<ny-1) push(r,r+nx);
        }
        A.rowPtr[A.rows]=A.nnz();
        return A;
    }
    // find where (row, col) sits in the values array, columns are sorted so binary search works
    int coefficientPosition(const CSRMatrix& A,int r,int c) {
        const auto a=A.colIndex.begin()+A.rowPtr.at(r), b=A.colIndex.begin()+A.rowPtr.at(r+1);
        const auto it=std::lower_bound(a,b,c);
        if(it==b||*it!=c) throw std::logic_error("Coefficient missing from CSR pattern.");
        return static_cast<int>(it-A.colIndex.begin());
    }
    void addCoefficient(CSRMatrix& A,int r,int c,double x) {
        A.values[coefficientPosition(A,r,c)]+=x;
    }
    void zeroValues(CSRMatrix& A) {
        std::fill(A.values.begin(),A.values.end(),0.0);
    }
    // sanity check of the CSR arrays
    void validateCSR(const CSRMatrix& A) {
        if(A.rows<1||A.cols<1||static_cast<int>(A.rowPtr.size())!=A.rows+1||A.rowPtr.front()!=0||A.rowPtr.back()!=A.nnz()||A.values.size()!=A.colIndex.size())
        throw std::invalid_argument("Invalid CSR dimensions.");
        for(int r=0;r<A.rows;++r) {
            if(A.rowPtr[r]<0 || A.rowPtr[r+1]>A.nnz() || A.rowPtr[r]>A.rowPtr[r+1])throw std::invalid_argument("Nonmonotone rowPtr.");
            int last=-1;
            for(int k=A.rowPtr[r];k<A.rowPtr[r+1];++k) {
                if(A.colIndex[k]<=last||A.colIndex[k]>=A.cols||!std::isfinite(A.values[k]))throw std::invalid_argument("Unsorted/nonfinite CSR entry.");
                last=A.colIndex[k];
            }
        }
    }
    // the pressure equation only knows p up to a constant, so we fix p = 0 in one cell
    void pinZeroReference(CSRMatrix& A,std::vector<double>& b,int ref) {
        if(A.rows!=A.cols||ref<0||ref>=A.rows||static_cast<int>(b.size())!=A.rows)throw std::invalid_argument("Invalid reference.");
        for(int r=0;r<A.rows;++r)for(int k=A.rowPtr[r];k<A.rowPtr[r+1];++k)
        if(r==ref||A.colIndex[k]==ref)A.values[k]=(r==ref&&A.colIndex[k]==ref)?1.0:0.0;
        b[ref]=0.0;
    }
    // transpose: count entries per column first, then fill
    CSRMatrix transposeCSR(const CSRMatrix& A) {
        CSRMatrix T(A.cols,A.rows);
        T.values.resize(A.nnz());
        T.colIndex.resize(A.nnz());
        for(int c:A.colIndex)++T.rowPtr[c+1];
        for(int r=0;r<T.rows;++r)T.rowPtr[r+1]+=T.rowPtr[r];
        auto next=T.rowPtr;
        for(int r=0;r<A.rows;++r)for(int k=A.rowPtr[r];k<A.rowPtr[r+1];++k) {
            int c=A.colIndex[k],q=next[c]++;
            T.colIndex[q]=r;
            T.values[q]=A.values[k];
            if(c==r)T.diagonalIndex[c]=q;
        }
        return T;
    }
    CSRMatrix productCSR(const CSRMatrix& A,const CSRMatrix& B) {
        if(A.cols!=B.rows)throw std::invalid_argument("Sparse product shape mismatch.");
        CSRMatrix C(A.rows,B.cols);
        // build one row at a time with a small map, so we never make a dense matrix
        for(int r=0;r<A.rows;++r) {
            C.rowPtr[r]=C.nnz();
            std::map<int,double> row;
            for(int k=A.rowPtr[r];k<A.rowPtr[r+1];++k) {
                if(A.values[k]==0)continue;
                int j=A.colIndex[k];
                for(int q=B.rowPtr[j];q<B.rowPtr[j+1];++q)if(B.values[q]!=0)row[B.colIndex[q]]+=A.values[k]*B.values[q];
            }
            for(const auto& kv:row)if(kv.second!=0||kv.first==r) {
                if(kv.first==r)C.diagonalIndex[r]=C.nnz();
                C.colIndex.push_back(kv.first);
                C.values.push_back(kv.second);
            }
        }
        C.rowPtr[C.rows]=C.nnz();
        return C;
    }
    // write a matrix and rhs to a text file (for checking by hand)
    void dumpCSR(const CSRMatrix& A,const std::vector<double>& b,const std::string& prefix) {
        std::ofstream f(prefix+".csr.txt");
        if(!f)throw std::runtime_error("Cannot dump CSR.");
        f<<std::setprecision(17)<<"rows "<<A.rows<<" cols "<<A.cols<<" nnz "<<A.nnz()<<"\nrowPtr\n";
        for(int a:A.rowPtr)f<<a<<' ';
        f<<"\ncolIndex\n";
        for(int a:A.colIndex)f<<a<<' ';
        f<<"\nvalues\n";
        for(double a:A.values)f<<a<<' ';
        f<<"\ndiagonalIndex\n";
        for(int a:A.diagonalIndex)f<<a<<' ';
        f<<"\nrhs\n";
        for(double a:b)f<<a<<' ';
        f<<'\n';
    }
    // max abs value, also stops the run if we get NaN or inf
    double normInf(const std::vector<double>& x) {
        double n=0;
        for(double a:x) {
            if(!std::isfinite(a))throw std::runtime_error("Nonfinite field/residual.");
            n=std::max(n,std::abs(a));
        }
        return n;
    }
}
