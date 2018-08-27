SET(CSM_SRC 
    csm/laser_data.cpp
    csm/math_utils.cpp
    csm/math_utils_gsl.cpp
    csm/utils.cpp
    csm/logging.cpp
    csm/orientation.cpp
    csm/clustering.cpp
    csm/laser_data_yogo.cpp
    csm/sm_config.cpp
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
    icp/icp_yogo.cpp
    icp/icp_optimization.cpp
)


ADD_LIBRARY(cms_plicp #SHARED
    ${CSM_SRC}
    ${EGSL_SRC}
    ${GPC_SRC}
    ${ICP_SRC}
)

add_executable(pl-icp pl-icp.cpp)
target_link_libraries(pl-icp  cms_plicp)

add_executable(icp_with_optimization icp_with_optimization.cpp)
target_link_libraries(icp_with_optimization
    cms_plicp
    ${CSPARSE_LIBRARY}
    ${G2O_LIBS}
)


#add_executable(plicp-yogocc  plicp-yogocc.cpp)
#target_link_libraries(plicp-yogocc cms_plicp)
