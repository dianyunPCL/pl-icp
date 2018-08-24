cmake_minimum_required(VERSION 2.8)

project(pl-icp)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall") 
set(CMAKE_CXX_FLAGS "-std=c++11 -O3" )

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin )
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib )
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)
find_package(G2O REQUIRED)
find_package(CSparse REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    ${G2O_INCLUDE_DIR}
    ${CSPARSE_INCLUDE_DIR}
)

message(STATUS "${CHOLMOD_INCLUDE_DIR}")
message(STATUS "${CSPARSE_INCLUDE_DIR}")
SET(G2O_LIBS
    g2o_cli g2o_core g2o_csparse_extension  g2o_ext_freeglut_minimal
    g2o_incremental g2o_interface g2o_interactive g2o_opengl_helper
    g2o_parser g2o_simulator g2o_solver_cholmod g2o_solver_csparse g2o_solver_dense
    g2o_solver_slam2d_linear g2o_solver_eigen g2o_solver_pcg g2o_solver_structure_only
    g2o_stuff g2o_types_data g2o_types_icp g2o_types_slam2d g2o_types_sba
    cxsparse
)

add_subdirectory(src)

