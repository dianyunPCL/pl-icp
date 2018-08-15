SET(CSM_SRC 
    csm/laser_data.cpp
    csm/math_utils.cpp
    csm/math_utils_gsl.cpp
    csm/utils.cpp
    csm/logging.cpp
    csm/sm_options.cpp
    csm/orientation.cpp
    csm/clustering.cpp
    csm/laser_data_yogo.cpp
)

SET(EGSL_SRC
    egsl/egsl.cpp
    egsl/egsl_conversions.cpp
    egsl/egsl_misc.cpp
    egsl/egsl_ops.cpp
)

SET(GPC_SRC
    gpc/gpc.cpp
    gpc/gpc_utils.cpp
)

SET(ICP_SRC
    icp/icp.cpp
    icp/icp_corr_dumb.cpp
    icp/icp_corr_tricks.cpp
    icp/icp_covariance.cpp
    icp/icp_debug.cpp
    icp/icp_loop.cpp
    icp/icp_outliers.cpp
    icp/pl_icp.cpp
    icp/icp_optimization.cpp
)

SET(OPTIONS_SRC  
    options/options.cpp
    options/options_interface.cpp
)

ADD_LIBRARY(plicp #SHARED
	${CSM_SRC}
    ${EGSL_SRC}
    ${GPC_SRC}
    ${ICP_SRC}
    ${OPTIONS_SRC}
)

add_executable(pl-icp pl-icp.cpp)
target_link_libraries(pl-icp plicp)

#add_executable(plicp-yogocc  plicp-yogocc.cpp)
#target_link_libraries(plicp-yogocc plicp)
