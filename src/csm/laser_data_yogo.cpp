#include <errno.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fstream>
#include <dirent.h>
#include <sys/types.h>
#include <iostream>
#include "csm/csm_all.h"
#include "csm/laser_data_yogo.h"

#define PI 3.1415926

using namespace  std;

//double theta[];

LDP ld_from_yogo_stream(const char* file, int nrays, int frame_id) {
    LDP ld = ld_alloc_new(nrays);
	
    double delta_a = PI / 720;
    double delta_odom[3];
    int imu_data;

    // open file
    std::ifstream fin(file);
    if (!fin) {
        sm_error("file %s can not be found!\n", file);
        return nullptr;
    }

    // get first line data
    fin >> delta_odom[0] >> delta_odom[1] >> delta_odom[2] >> imu_data;
    // get ranges data
    double min_angle = -(nrays - 1) * 0.5 * delta_a;
    ld->max_theta = -min_angle;
    ld->min_theta = min_angle;
    ld->odometry[0] = 0.0;
    ld->odometry[1] = 0.0;
    ld->odometry[2] = imu_data * 0.01 * PI / 180;
    ld->frame_id = frame_id;

    int range = 0.;
    bool flag = 0;
    int nonse_data[2];
    int i = 0;
    double t;
    while (!fin.eof()) {
        fin >> range >> flag >> nonse_data[0] >> nonse_data[1];
        ld->readings[i] = range * 0.001;
        ld->valid[i] = flag;
        t = min_angle + i * delta_a;
        ld->theta[i] = t;

        i++;
    }

	return ld;
}


LDP ld_from_keyframe(laser_data laserD) {
    int nrays = laserD.nrays;
    LDP ld = ld_alloc_new(nrays);

    ld->frame_id = laserD.frame_id;
    ld->max_theta = laserD.max_theta;
    ld->min_theta = laserD.min_theta;

    memcpy(ld->valid, laserD.valid, sizeof(int)*nrays);
    memcpy(ld->readings, laserD.readings, sizeof(double)*nrays);
    memcpy(ld->readings_sigma, laserD.readings_sigma, sizeof(double)*nrays);
    memcpy(ld->theta, laserD.theta, sizeof(double)*nrays);
    memcpy(ld->cluster, laserD.cluster, sizeof(int)*nrays);
    memcpy(ld->alpha, laserD.alpha, sizeof(double)*nrays);
    memcpy(ld->cov_alpha, laserD.cov_alpha, sizeof(double)*nrays);
    memcpy(ld->alpha_valid, laserD.alpha_valid, sizeof(int)*nrays);
    memcpy(ld->true_alpha, laserD.true_alpha, sizeof(double)*nrays);
    memcpy(ld->up_bigger, laserD.up_bigger, sizeof(int)*nrays);
    memcpy(ld->up_smaller, laserD.up_smaller, sizeof(int)*nrays);
    memcpy(ld->down_bigger, laserD.down_bigger, sizeof(int)*nrays);
    memcpy(ld->down_smaller, laserD.down_smaller, sizeof(int)*nrays);
    memcpy(ld->corr, laserD.corr, sizeof(struct correspondence)*nrays);

    int i;
    for(i=0; i<nrays; ++i) {
        ld->corr[i].valid = laserD.corr[i].valid;
        ld->corr[i].j1 = laserD.corr[i].j1;
        ld->corr[i].j2 = laserD.corr[i].j2;
        ld->points[i].p[0] = laserD.points[i].p[0];
        ld->points[i].p[1] = laserD.points[i].p[1];
        ld->points[i].rho = laserD.points[i].rho;
        ld->points[i].phi = laserD.points[i].phi;
        ld->points_w[i].p[0] = laserD.points_w[i].p[0];
        ld->points_w[i].p[1] = laserD.points_w[i].p[1];
        ld->points_w[i].rho = laserD.points_w[i].rho;
        ld->points_w[i].phi = laserD.points_w[i].phi;
    }
    for(i=0; i<3; ++i) {
        ld->odometry[i] = laserD.odometry[i];
        ld->estimate[i] = laserD.estimate[i];
        ld->global_pose[i] = laserD.global_pose[i];
        ld->true_pose[i] = laserD.true_pose[i];
    }

    return ld;
}

