// lift and drag on the body
#pragma once
#include "cfd/flow_fields.h"
#include "cfd/body_geometry.h"
namespace cfd {
    // force on the body per unit span with rho = 1. free stream goes along +x so drag = fx and lift = fy
    struct BodyForces {
        double fx_pressure=0,fx_viscous=0,fy_pressure=0,fy_viscous=0;
        double cd=0,cl=0;
        // coefficient = force / (0.5*rho*U^2*c), here rho = U = c = 1
    };
    // adds up pressure and wall shear on every face between a fluid cell and a solid cell
    BodyForces computeBodyForces(const StaggeredGrid& g,const BodyMask& body,const FlowFields& f,double reynolds);
}
