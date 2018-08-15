#include "csm/csm_all.h"
#include "icp/icp_optimization.h"
#include <iostream>

using namespace std;

/** @brief find new edges to optimize the icp results
 *  @param[in]      params      icp parameters
 *  @param[in|out]  edges       g2o edges
 *  @param[in]      keyframes   keyframes

 */
void find_new_edges(sm_params& params, std::vector<icp_edge>& edges, std::vector<laser_data>& keyframes)
{
    printf("[pl-icp] finding suitable edges... \n");

    sm_result result;
    icp_edge edge;

    int nDisFilter, nAngFilter, nIterFilter, newEdges, tatal_count;
    nDisFilter = nAngFilter =  nIterFilter = newEdges = tatal_count = 0;
    double da1, da2, da;
    double dx, dy, d;
    LDP kf1 = ld_alloc_new(937);
    LDP kf2 = ld_alloc_new(937);
    int num_frames = keyframes.size();
    for (int i=1; i<num_frames-15; i+=10) {
        kf1 = ld_from_keyframe(keyframes[i]);
        for (int j=i+15; j<num_frames; j+=10) {
            ++tatal_count;
            kf2 = ld_from_keyframe(keyframes[j]);

            // delta_distance < 2 meters to include
            // delta_angule < 30 degree to include
            // delta_frame > 15 frames to include
            // iteration < 30 times to include
            dx = abs(kf1->global_pose[0] - kf2->global_pose[0]);
            dy = abs(kf1->global_pose[1] - kf2->global_pose[1]);
            d = dx * dx + dy * dy;
            if (d > 4) {
                ++nDisFilter;
                free(kf2);
                continue;
            }
            da1 = kf1->global_pose[2];
            da2 = kf2->global_pose[2];
            da = abs(da1 - da2);
            if (da > 0.5236) {
                ++nAngFilter;
                free(kf2);
                continue;
            }

            params.laser_ref = kf1;
            params.laser_sens = kf2;
            double odometry[3], ominus_laser[3], temp[3];
            pose_diff_d(kf2->odometry, kf1->odometry, odometry);
            ominus_d(params.laser, ominus_laser);
            oplus_d(ominus_laser, odometry, temp);
            oplus_d(temp, params.laser, params.first_guess);

            sm_icp(&params, &result);
            if (result.iterations > 30) {
              ++nIterFilter;
              free(kf2);
              continue;
            }

            /** finally a new edge born */
            ++newEdges;
            edge.edge_id_ahead = i;
            edge.edge_id_behind = j;
            copy_d(result.x, 3, edge.transform);
//            double w = compute_weight(kf1.global_pose, kf2.global_pose, result.x);

            // use fixed weight
            double w = 1.0;
            edge.information[0] = 20 * w;
            edge.information[1] = 5 * w;
            edge.information[2] = 286.4789 * w; // 900/pi, 286.4789
            edges.push_back(edge);

            free(kf2);
        }
        free(kf1);
    }
    // print edge selection results
    printf("[pl-icp] finding suitable edges done. tatal count: %d \n", tatal_count);
    printf("[pl-icp] discard edge matching through distance: %d(%f%%)\n",
           nDisFilter, (double)100*nDisFilter/tatal_count);
    printf("[pl-icp] discard edge matching through angle: %d(%f%%)\n",
           nAngFilter, (double)100*nAngFilter/(tatal_count-nDisFilter));
    printf("[pl-icp] discard edge matching through iterations: %d(%f%%)\n",
           nIterFilter, (double)100*nIterFilter/(tatal_count-nDisFilter-nAngFilter));
    printf("[pl-icp] Add tatal new edges: %d\n", newEdges);
}

/** compute weight for edge's information matrix */
double compute_weight(double* gp_i, double* gp_j, double* trans)
{
    double da = gp_j[2] - gp_i[2];
    double dx = gp_j[0] - gp_i[0] * cos(da) + gp_i[1] * sin(da);
    double dy = gp_j[1] - gp_i[0] * sin(da) - gp_i[1] * cos(da);
    double c = (dx - trans[0]) * (dx - trans[0]) +
               (dy - trans[1]) * (dy - trans[1]) +
               (da - trans[2]) * (da - trans[2]);
    double w = 0.5 / (0.5 + c);

    return w;
}
