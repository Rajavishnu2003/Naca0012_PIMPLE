cmake_minimum_required(VERSION 3.16)
project(NACA0012_Re1000_PIMPLE LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(ENABLE_SANITIZERS "Address and undefined-behaviour sanitizers (GCC/Clang)" OFF)
if(ENABLE_SANITIZERS AND NOT MSVC)
  add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address,undefined)
endif()
add_library(csr_linalg src/csr_matrix.cpp src/linear_algebra.cpp src/ic0_preconditioner.cpp
    src/conjugate_gradient.cpp src/bicgstab.cpp src/steepest_descent.cpp src/smoother.cpp)
target_include_directories(csr_linalg PUBLIC include)
set(CFD_SOURCES src/cfd/config.cpp src/cfd/staggered_grid.cpp src/cfd/flow_fields.cpp
    src/cfd/boundary_conditions.cpp src/cfd/csr_tools.cpp src/cfd/convection_scheme.cpp
    src/cfd/momentum_assembler.cpp src/cfd/pressure_assembler.cpp src/cfd/pressure_multigrid.cpp
    src/cfd/linear_solver_adapter.cpp src/cfd/convergence_monitor.cpp src/cfd/coupling_solver.cpp
    src/cfd/simple_solver.cpp src/cfd/piso_solver.cpp src/cfd/pimple_solver.cpp src/cfd/output_writer.cpp
    src/cfd/body_geometry.cpp src/cfd/body_forces.cpp)
add_library(cfd_core ${CFD_SOURCES})
target_link_libraries(cfd_core PUBLIC csr_linalg)
if(MSVC)
  target_compile_options(cfd_core PRIVATE /W4)
else()
  target_compile_options(cfd_core PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(airfoil mains/run_airfoil.cpp)
target_link_libraries(airfoil PRIVATE cfd_core)
