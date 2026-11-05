// Convection-diffusion coefficients for one face of a control volume.
#pragma once
#include <string>
#include <cstdint>
namespace cfd {
    // counts how many faces used central differencing and how many upwind
    struct SchemeStatistics {
        std::uint64_t cds_faces=0,upwind_faces=0;
        double max_abs_peclet=0;
        void merge(const SchemeStatistics& s);
    };
    // what this face adds to the diagonal and to the neighbour coefficient
    struct FaceCoefficients {
        double diagonal,neighbour;
        bool cds;
        double peclet;
    };
    // F is the mass flux through the face, positive going out of this control volume
    // D = mu * area / distance, the diffusion part is kept in full for both central and upwind
    FaceCoefficients faceCoefficients(double F,double D,const std::string& scheme="peclet");
}
