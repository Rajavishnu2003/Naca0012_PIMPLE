// writing results: csv, vtk for ParaView and summary.json
#pragma once
#include "cfd/coupling_solver.h"
#include <string>
#include <utility>
#include <vector>
namespace cfd {
    // final output at the end of the run
    void writeResults(const StaggeredGrid& g,const BodyMask& body,const FlowFields& f,const Config& c,const RunSummary& s,
    const MomentumAssembler& m,const PressureAssembler& p);
    // one frame for the animation: binary vtk (float) with velocity, pressure, vorticity and the solid mask
    void writeVTKSnapshot(const StaggeredGrid& g,const BodyMask& body,const FlowFields& f,const std::string& path,double time);
    // small json file for ParaView that lists each frame with its time
    void writeVTKSeries(const std::string& path,const std::vector<std::pair<std::string,double>>& frames);
}
