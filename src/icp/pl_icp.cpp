#include "icp/pl_icp.h"
#include "csm/csm_all.h"
#include "csm/math_utils.h"

#include <iostream>

using namespace std;

#define PI 3.1415926

sm_params plicp_param;
sm_result plicp_result;
//struct global_position gp;

/** input data */
LDP set_laser_frame(const int range_[], const char flags_[],
                     const int point_num_, const double odom_[])
{
    // 数据赋值
    LDP laser = ld_alloc_new(point_num_);

    double delta_a = PI / 720;
    int nrays = point_num_ - 1;

    laser->max_theta = 0.5 * nrays * delta_a;
    laser->min_theta = -laser->max_theta;

    for (int i=0; i<point_num_; ++i) {
        laser->readings[i] = (double)range_[i] * 0.001;   // mm to m
        laser->valid[i] = (int)flags_[i];
        laser->theta[i] = laser->min_theta + i * delta_a;
    }

    copy_d(odom_, 3, laser->odometry);
    copy_d(odom_, 3, laser->estimate);

    return laser;
}

/** set params */
void set_plicp_params(LDP laser_ref, LDP laser_curr)
{
    plicp_param.laser_ref = laser_ref;
    plicp_param.laser_sens = laser_curr;
    plicp_param.max_angular_correction_deg = 90.0;
    plicp_param.max_linear_correction = 2.0;
    plicp_param.max_iterations = 50;    // 1000
    plicp_param.epsilon_xy = 0.0001;    // A threshold for stopping icp loop
    plicp_param.epsilon_theta = 0.0001; // A threshold for stopping icp loop
    plicp_param.max_correspondence_dist = 0.5;
    plicp_param.sigma = 0.01;
    plicp_param.use_corr_tricks = 1;
    plicp_param.restart = 1;
    plicp_param.restart_threshold_mean_error = 0.01;
    plicp_param.restart_dt = 0.01;
    plicp_param.restart_dtheta = deg2rad(1.5);
    plicp_param.clustering_threshold = 0.05;
    plicp_param.orientation_neighbourhood = 3;
    plicp_param.use_point_to_line_distance = 1;
    plicp_param.do_alpha_test = 0;
    plicp_param.do_alpha_test_thresholdDeg = 20.0;
    plicp_param.outliers_maxPerc = 0.95;
    plicp_param.outliers_adaptive_mult = 2.0;
    plicp_param.outliers_adaptive_order = 0.7;
    plicp_param.do_visibility_test = 0;
    plicp_param.outliers_remove_doubles = 1;
    plicp_param.do_compute_covariance = 1;
    plicp_param.debug_verify_tricks = 0;
    plicp_param.min_reading = 0.0;
    plicp_param.max_reading = 10.0;
    plicp_param.use_ml_weights = 0;
    plicp_param.use_sigma_weights = 0;
    plicp_param.laser[0] = 0.0;
    plicp_param.laser[1] = 0.0;
    plicp_param.laser[2] = 0.0;
}

/** do pl-icp and return a delta odometry */
void do_plicp(vector<double>& result, Eigen::Matrix3d& covariance)
{
    double odometry[3];
    pose_diff_d(plicp_param.laser_sens->odometry, plicp_param.laser_ref->odometry, odometry);
    double ominus_laser[3], temp[3];
    ominus_d(plicp_param.laser, ominus_laser);
    oplus_d(ominus_laser, odometry, temp);
    oplus_d(temp, plicp_param.laser, plicp_param.first_guess);

    // pl-icp
    clock_t tt = clock();
    sm_icp(&plicp_param, &plicp_result);
    double t = 1000 * (clock() - tt) / (double)CLOCKS_PER_SEC;

//    plicp_result.cost_time = t;

    result.resize(5);
    result[0] = plicp_result.x[0];
    result[1] = plicp_result.x[1];
    result[2] = plicp_result.x[2];
    result[3] = t;
    result[4] = plicp_result.iterations;

    if (plicp_param.do_compute_covariance) {
        printf("-icp- do_compute_covariance.");
        covariance = *plicp_result.cov_x_m;
    }

    return;
}

/** get global position */
std::vector<double> get_global_pose(const double* gp_ref,const double* delta_trans)
{
    double dx = delta_trans[0];
    double dy = delta_trans[1];
    double da = delta_trans[2];
//    printf("delta_trans: %f, %f, %f\n", dx, dy, da);
    double gp_curr[3];
    gp_curr[0] = cos(gp_ref[2])*dx - sin(gp_ref[2])*dy + gp_ref[0];
    gp_curr[1] = sin(gp_ref[2])*dx + cos(gp_ref[2])*dy + gp_ref[1];
    gp_curr[2] = gp_ref[2] + da;
    if (gp_curr[2] > PI) gp_curr[2] -= 2 * PI;
    if (gp_curr[2] < -PI) gp_curr[2] += 2 * PI;

    vector<double> global_pos(3);
    global_pos[0] = gp_curr[0];
    global_pos[1] = gp_curr[1];
    global_pos[2] = gp_curr[2];

    return global_pos;
}

vector<double> get_global_position(const double* delta_trans)
{
    double dx = delta_trans[0];
    double dy = delta_trans[1];
    double da = delta_trans[2];
//    printf("delta_trans: %f, %f, %f\n", dx, dy, da);
    double gp_ref[3], gp_curr[3];
    gp_ref[0] = plicp_param.laser_ref->global_pose[1];
    gp_ref[1] = plicp_param.laser_ref->global_pose[1];
    gp_ref[2] = plicp_param.laser_ref->global_pose[2];
//    copy_d(plicp_param.laser_ref->global_pose, 3, gp_ref);

//    gp.x += cos(gp.theta)*dx - sin(gp.theta)*dy;
//    gp.y += sin(gp.theta)*dx + cos(gp.theta)*dy;
//    gp.theta += da;
    gp_curr[0] = cos(gp_ref[2])*dx - sin(gp_ref[2])*dy;
    gp_curr[1] = sin(gp_ref[2])*dx + cos(gp_ref[2])*dy;
    gp_curr[2] = gp_ref[2] + da;
    if (gp_curr[2] > PI) gp_curr[2] -= 2 * PI;
    if (gp_curr[2] < -PI) gp_curr[2] += 2 * PI;

    vector<double> global_pos(3);
    global_pos[0] = gp_curr[0];
    global_pos[1] = gp_curr[1];
    global_pos[2] = gp_curr[2];

    return global_pos;
}


