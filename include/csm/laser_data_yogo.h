#ifndef H_LASER_DATA_JSON
#define H_LASER_DATA_JSON

#include "gsl_eigen/gsl_eigen.h"
#include "laser_data.h"
#include "algos.h"


/** Laserdata from yogo data format */
LDP ld_from_yogo_stream(const char* file_path, int nrays, int frame_id);

/** get a LDP from keyframe (laser_data) */
LDP ld_from_keyframe(laser_data laserD);

#endif
