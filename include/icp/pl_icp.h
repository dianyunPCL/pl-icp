#ifndef _PL_ICP_
#define _PL_ICP_

#include <vector>

#include "csm/csm_all.h"
#include "csm/utils.h"


//struct global_position
//{
//    double x;
//    double y;
//    double theta;

//    global_position() {
//        x = y = theta = 0.;
//    }
//};

/** input data */
LDP set_laser_frame(const int range_[], const char flags_[],
                    const int point_num_, const double odom_[]);

/** set params */
void set_plicp_params(LDP laser_ref, LDP laser_curr);

/** do pl-icp and return a delta odometry */
void do_plicp(std::vector<double>& result, Eigen::Matrix3d& covariance);

/** get global position */
std::vector<double> get_global_position(const double* delta_trans);
std::vector<double> get_global_pose(const double* gp_ref,const double* delta_trans);

#endif
