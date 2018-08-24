#include "icp/icp_yogo.h"
#include "csm/csm_all.h"
#include "csm/math_utils.h"

#include <iostream>

using namespace std;

#define PI 3.1415926


/** input data */
LDP set_laser_frame(const int range[], const char flags[],
                    const double odom[], int nrays, int frame_id)
{
    double delta_a = PI / 720;

    // 数据赋值
    LDP laser = ld_alloc_new(nrays);
    laser->frame_id = frame_id;
    laser->max_theta = 0.5 * (nrays-1) * delta_a;
    laser->min_theta = -laser->max_theta;

    for (int i=0; i<nrays; ++i) {
        laser->readings[i] = (double)range[i] * 0.001;   // mm to m
        laser->valid[i] = (int)flags[i];
        laser->theta[i] = laser->min_theta + i * delta_a;
    }

    copy_d(odom, 3, laser->odometry);
    copy_d(odom, 3, laser->estimate);

    return laser;
}

/** set params */
void set_plicp_params(sm_params* plicp_param)
{
    plicp_param->max_angular_correction_deg = 45.0;
    plicp_param->max_linear_correction = 0.5;
    plicp_param->max_iterations = 50;    // 1000
    plicp_param->epsilon_xy = 0.0001;    // A threshold for stopping icp loop
    plicp_param->epsilon_theta = 0.005; // A threshold for stopping icp loop
    plicp_param->max_correspondence_dist = 0.25;
    plicp_param->sigma = 0.005;
    plicp_param->use_corr_tricks = 1;
    plicp_param->restart = 1;
    plicp_param->restart_threshold_mean_error = 0.01;
    plicp_param->restart_dt = 0.01;
    plicp_param->restart_dtheta = deg2rad(1.5);
    plicp_param->clustering_threshold = 0.05;
    plicp_param->orientation_neighbourhood = 10;
    plicp_param->use_point_to_line_distance = 1;
    plicp_param->do_alpha_test = 0;
    plicp_param->do_alpha_test_thresholdDeg = 20.0;
    plicp_param->outliers_maxPerc = 0.95;
    plicp_param->outliers_adaptive_mult = 2.0;
    plicp_param->outliers_adaptive_order = 0.7;
    plicp_param->do_visibility_test = 0;
    plicp_param->outliers_remove_doubles = 1;
    plicp_param->do_compute_covariance = 1;
    plicp_param->debug_verify_tricks = 0;
    plicp_param->min_reading = 0.0;
    plicp_param->max_reading = 10.0;
    plicp_param->use_ml_weights = 0;
    plicp_param->use_sigma_weights = 0;
    plicp_param->laser[0] = 0.0;
    plicp_param->laser[1] = 0.0;
    plicp_param->laser[2] = 0.0;
    plicp_param->kf_delta_frame = 1;
    plicp_param->kf_dist_linear = 0.025;
    plicp_param->kf_dist_angular = 0.05;
    plicp_param->pg_max_iterations = 5;
    plicp_param->pg_max_frames = 100;
}

/** get global position */
void get_global_pose(double* global_pos, const double gp_ref[], const double delta_trans[])
{
    double dx = delta_trans[0];
    double dy = delta_trans[1];
    double da = delta_trans[2];
    double gp_curr[3];
    gp_curr[0] = cos(gp_ref[2])*dx - sin(gp_ref[2])*dy + gp_ref[0];
    gp_curr[1] = sin(gp_ref[2])*dx + cos(gp_ref[2])*dy + gp_ref[1];
    gp_curr[2] = gp_ref[2] + da;
    if (gp_curr[2] > PI) gp_curr[2] -= 2 * PI;
    if (gp_curr[2] < -PI) gp_curr[2] += 2 * PI;

    copy_d(gp_curr, 3, global_pos);
}

bool valid_transform(sm_params* params, double* delta_trans)
{
    double dx = delta_trans[0];
    double dy = delta_trans[1];
    double da = delta_trans[2];
    double dsq = dx * dx + dy * dy;

    // deal with too large motion
    double max_dis = params->kf_delta_frame * 0.2;
    if ( dsq > max_dis*max_dis || abs(da) > 1.0) {
        cout << "[icp][warning] too large motion: " << sqrt(dsq)
             << ", abs da: " << abs(da) << ", this may be a wrong match!" << endl;
        return false;
    }

    // deal with too small motion
    // set 0.067 mm and 0.05 degree transform to 0.
    double zero_trans[3] = {0., 0., 0.};
    if ( dsq < 4.5e-04 && abs(da) < 8.73e-04)
        copy_d(zero_trans, 3, delta_trans);

    return true;
}
