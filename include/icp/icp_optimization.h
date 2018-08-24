#ifndef ICP_OPTIMIZATION_H
#define ICP_OPTIMIZATION_H

#include "csm/csm_all.h"
#include "icp/icp.h"
#include "csm/laser_data_yogo.h"

using std::vector;

struct icp_vertex {
    int vertex_id;      // local
    int keyframe_id;    // global
    double vertex_pose[3];
};

struct icp_edge {
    int edge_id_to, edge_id_from;
    double transform[3];
    double information[3]; // main value of information matrix
};

void find_new_edges(sm_params* params, std::vector<icp_edge>& edges, std::vector<laser_data>& keyframes);

double compute_weight_new_edge(const double gp_i[], const double gp_j[],
                               const double trans[], int iterations);
double compute_weight_neighbor_edge(int iterations);

bool newKeyframeNeeded(sm_params* param, double* trans);

void add_extra_edges(sm_params* params, std::vector<icp_edge>& edges,
                     std::vector<laser_data>& keyframes, int begin_index, int last_index);

bool needOptimization(sm_params* params, vector<laser_data>& keyframes,
                      int begin_index, int last_index);

#endif // ICP_OPTIMIZATION_H
