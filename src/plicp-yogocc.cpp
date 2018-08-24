#include <time.h>
#include <string.h>
#include <libgen.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <dirent.h>
#include <ctime>
#include <cstdlib>

#include "csm/laser_data_yogo.h"
#include "icp/icp_yogo.h"

using namespace std;

int main(int argc, const char*argv[]) {
    // generate random data for test
    int num = 1081;
    int  range_ref[num];
    char flags_ref[num];
    int  range_curr[num];
    char flags_curr[num];
    double delta_odom_ref[3] = {0., 0., 0.};
    double delta_odom_curr[3] = {.8, .1, .05};
    srand((unsigned)time(0) );
    for (int i=0; i<num; ++i) {
        range_ref[i] = rand()%9000;
        range_curr[i] = range_ref[i] + 1000;
        flags_ref[i] = (char)abs(rand()%2);
        flags_curr[i] = flags_ref[i];
    }

    // input data
    LDP laser_ref = set_laser_frame(range_ref, flags_ref,
                                    delta_odom_ref, num, 0);
    LDP laser_curr = set_laser_frame(range_curr, flags_curr,
                                     delta_odom_curr, num, 1);
    // check data
    if(!ld_valid_fields(laser_ref))  {
        sm_error(" -icp- Invalid laser data in first scan.\n");
        return -2;
    }
    if(!ld_valid_fields(laser_curr))  {
        sm_error(" -icp- Invalid laser data in second scan.\n");
        return -2;
    }
    if(	any_nan(laser_ref->odometry,3) ||
        any_nan(laser_curr->odometry,3) ) {
        printf("odometry NAN.!!\n");
        return -3;
    }

    printf(" -icp- set input data sucessed!\n");

    sm_params params;
    sm_result result;
    set_plicp_params(&params);

    sm_icp(&params, &result);

    // output result
    printf(" -icp- delta transform : [%f, %f, %f]\n", result.x[0], result.x[1], result.x[2]);
    printf(" -icp- iterations: %d\n", result.iterations);

    // free
    ld_free(laser_ref);
    ld_free(laser_curr);

    return 0;
}
