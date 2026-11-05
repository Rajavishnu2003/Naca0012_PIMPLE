// reads the key = value case file, checks the values and writes a copy into the output folder
#include "cfd/config.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <set>
#include <iomanip>
namespace cfd {
    namespace {
        std::string trim(const std::string& s) {
            const auto a=s.find_first_not_of(" \t\r\n");
            if(a==std::string::npos) return "";
            return s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);
        }
        // read one value, gives an error if there is anything left over
        template<class T> void parse(const std::string& s,T& x) {
            std::istringstream in(s);
            char extra;
            if(!(in>>x) || (in>>extra)) throw std::invalid_argument("Invalid value: "+s);
        }
        void parse(const std::string& s,bool& x) {
            if(s=="true"||s=="1") x=true;
            else if(s=="false"||s=="0") x=false;
            else throw std::invalid_argument("Expected true/false: "+s);
        }
        void parse(const std::string& s,std::string& x) {
            x=s;
        }
    }
    // stop early with a clear message if a setting does not make sense
    void Config::validate() const {
        if(nx<3||ny<3||nx>2048||ny>2048) throw std::invalid_argument("nx and ny must be in [3,2048].");
        auto positive=[](double x) {
            return std::isfinite(x)&&x>0;
        };
        if(!positive(reynolds)||!positive(dt)||!positive(end_time)||
        !positive(momentum_tolerance)||!positive(continuity_tolerance)||!positive(change_tolerance)||
        !positive(linear_rtol)||!positive(linear_atol)) throw std::invalid_argument("Physical values and tolerances must be finite and positive.");
        if(!std::isfinite(vtk_start_time)||vtk_start_time<0) throw std::invalid_argument("vtk_start_time must be >= 0.");
        if(!std::isfinite(pressure_rtol)||pressure_rtol<0||pressure_rtol>=1) throw std::invalid_argument("pressure_rtol must be in [0,1).");
        if(!positive(alpha_u)||alpha_u>1||!positive(alpha_p)||alpha_p>1)
        throw std::invalid_argument("Relaxation factors must be in (0,1].");
        if(algorithm!="simple"&&algorithm!="piso"&&algorithm!="pimple") throw std::invalid_argument("Unknown coupling algorithm.");
        if(convection!="peclet"&&convection!="cds"&&convection!="upwind") throw std::invalid_argument("Unknown convection scheme.");
        if(pressure_solver!="cg_ic0"&&pressure_solver!="cg"&&pressure_solver!="sd_ic0"&&pressure_solver!="bicgstab"&&pressure_solver!="multigrid")
        throw std::invalid_argument("Unknown pressure solver.");
        if(mg_cycle!="V"&&mg_cycle!="W"&&mg_cycle!="F") throw std::invalid_argument("MG cycle must be V, W or F.");
        if(max_steps<1||max_linear_iterations<1||pressure_correctors<2||outer_correctors<1||mg_pre<1||mg_post<1||log_every<1||vtk_every<0)
        throw std::invalid_argument("Invalid iteration or smoothing count; PISO requires >=2 correctors.");
        if(algorithm=="simple"&&!stop_at_steady) throw std::invalid_argument("SIMPLE is steady; use PISO/PIMPLE for end_time runs.");
        if(output.empty()) throw std::invalid_argument("Output directory is empty.");
        if(case_type!="cavity"&&case_type!="external") throw std::invalid_argument("case_type must be cavity or external.");
        if(case_type=="external") {
            if(body!="naca"&&body!="cylinder"&&body!="none") throw std::invalid_argument("body must be naca, cylinder or none.");
            if(naca.size()!=4||naca.find_first_not_of("0123456789")!=std::string::npos||naca.substr(2)=="00")
            throw std::invalid_argument("naca must be a 4-digit code with nonzero thickness, e.g. 0012 or 2412.");
            if(naca[0]!='0'&&naca[1]=='0') throw std::invalid_argument("Cambered NACA code needs a nonzero camber position (second digit).");
            if(!std::isfinite(aoa_deg)||std::abs(aoa_deg)>30) throw std::invalid_argument("aoa_deg must be in [-30,30].");
            if(!(x_max>x_min)||!(y_max>y_min)) throw std::invalid_argument("Domain needs x_max>x_min and y_max>y_min.");
            if(pressure_solver=="multigrid") throw std::invalid_argument("Geometric multigrid assumes a full rectangle; use cg_ic0 or bicgstab with a body.");
        }
    }
    Config readConfig(const std::string& path) {
        std::ifstream f(path);
        if(!f) throw std::runtime_error("Cannot open config: "+path);
        Config c;
        std::string line;
        int number=0;
        std::set<std::string> used;
        while(std::getline(f,line)) {
            ++number;
            line=trim(line.substr(0,line.find('#')));  // anything after # is a comment
            if(line.empty()) continue;
            const auto at=line.find('=');
            if(at==std::string::npos) throw std::invalid_argument("Expected key=value on line "+std::to_string(number));
            const auto k=trim(line.substr(0,at)),v=trim(line.substr(at+1));
            if(!used.insert(k).second) throw std::invalid_argument("Duplicate config key: "+k);
            // match the key to the struct member with the same name
#define KEY(name) if(k==#name) {parse(v,c.name);continue;}
            KEY(case_type) KEY(nx) KEY(ny) KEY(reynolds) KEY(algorithm) KEY(pressure_solver) KEY(convection)
            KEY(alpha_u) KEY(alpha_p) KEY(dt) KEY(pressure_correctors) KEY(outer_correctors)
            KEY(max_steps) KEY(max_linear_iterations) KEY(momentum_tolerance) KEY(continuity_tolerance)
            KEY(change_tolerance) KEY(linear_rtol) KEY(pressure_rtol) KEY(linear_atol) KEY(stop_at_steady) KEY(end_time)
            KEY(mg_cycle) KEY(mg_pre) KEY(mg_post) KEY(mg_fmg_start) KEY(log_every) KEY(vtk_every) KEY(vtk_start_time) KEY(output) KEY(dump_matrices)
            KEY(body) KEY(naca) KEY(aoa_deg) KEY(x_min) KEY(x_max) KEY(y_min) KEY(y_max)
#undef KEY
            throw std::invalid_argument("Unknown config key: "+k);
        }
        c.validate();
        return c;
    }
    // saved as input.cfg so every result folder has the exact settings used
    void writeConfig(const Config& c,const std::string& path) {
        std::ofstream f(path);
        if(!f) throw std::runtime_error("Cannot write config: "+path);
        f<<std::setprecision(17)<<std::boolalpha;
#define OUT(name) f<<#name<<" = "<<c.name<<'\n';
        OUT(case_type) OUT(nx) OUT(ny) OUT(reynolds) OUT(algorithm) OUT(pressure_solver) OUT(convection)
        OUT(alpha_u) OUT(alpha_p) OUT(dt) OUT(pressure_correctors) OUT(outer_correctors)
        OUT(max_steps) OUT(max_linear_iterations) OUT(momentum_tolerance) OUT(continuity_tolerance)
        OUT(change_tolerance) OUT(linear_rtol) OUT(pressure_rtol) OUT(linear_atol) OUT(stop_at_steady) OUT(end_time)
        OUT(mg_cycle) OUT(mg_pre) OUT(mg_post) OUT(mg_fmg_start) OUT(log_every) OUT(vtk_every) OUT(vtk_start_time) OUT(output) OUT(dump_matrices)
        OUT(body) OUT(naca) OUT(aoa_deg) OUT(x_min) OUT(x_max) OUT(y_min) OUT(y_max)
#undef OUT
    }
}
