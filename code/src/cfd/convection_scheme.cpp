// Face coefficients for convection + diffusion (central, upwind or hybrid by Peclet number).
#include "cfd/convection_scheme.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace cfd {
    // add up the face counts from u and v so we can print how many faces went upwind
    void SchemeStatistics::merge(const SchemeStatistics& s) {
        cds_faces+=s.cds_faces;
        upwind_faces+=s.upwind_faces;
        max_abs_peclet=std::max(max_abs_peclet,s.max_abs_peclet);
    }
    FaceCoefficients faceCoefficients(double F,double D,const std::string& scheme) {
        if(!std::isfinite(F)||!std::isfinite(D)||D<=0)throw std::invalid_argument("Invalid face flux/conductance.");
        if(scheme!="peclet"&&scheme!="cds"&&scheme!="upwind")throw std::invalid_argument("Invalid scheme.");
        // cell Peclet number = convection / diffusion at this face
        const double Pe=F/D;
        // "peclet" (hybrid) uses central differencing only when |Pe| <= 2, above that CDS can give negative coefficients
        const bool cds=scheme=="cds"||(scheme=="peclet"&&std::abs(Pe)<=2.0);
        // central: face value is the average of the two cells, so F/2 goes to each side
        if(cds)return {
            D+0.5*F,-D+0.5*F,true,Pe
        };
        // upwind: face value is taken from the cell the flow comes from
        return {
            D+std::max(F,0.0),-D+std::min(F,0.0),false,Pe
        };
    }
}
