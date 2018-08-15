#ifndef ICP_OPTIMIZATION_H
#define ICP_OPTIMIZATION_H

#include "csm/csm_all.h"
#include "icp/icp.h"
#include "csm/laser_data_yogo.h"

struct icp_vertex {
    int vertex_id;
    double global_pose[3];
};

struct icp_edge {
    int edge_id_ahead, edge_id_behind;
    double transform[3];
    double information[3]; // main value of information matrix
};

void find_new_edges(sm_params& params, std::vector<icp_edge>& edges, std::vector<laser_data>& key_frames);

double compute_weight(double* gp_i, double* gp_j, double* trans);

#endif // ICP_OPTIMIZATION_H
