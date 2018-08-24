#ifndef _PL_ICP_
#define _PL_ICP_

#include <vector>

#include "csm/csm_all.h"
#include "csm/utils.h"

/** input data */
LDP set_laser_frame(const int range[], const char flags[],
                    const double odom[], int nrays, int frame_id);

/** set params */
void set_plicp_params(sm_params* params);

/** do pl-icp and return a delta odometry */
//void do_plicp(std::vector<double>& result, Eigen::Matrix3d& covariance);

/** get global position */
void get_global_pose(double* global_pos, const double gp_ref[], const double delta_trans[]);

bool valid_transform(sm_params* params, double* delta_trans);


#endif
