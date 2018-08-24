#include "csm/csm_all.h"
#include "icp/icp_optimization.h"
#include "icp/icp_yogo.h"
#include <iostream>

using namespace std;

/** find new edges to optimize the icp results (abandon)  */
void find_new_edges(sm_params* params, vector<icp_edge>& edges, vector<laser_data>& keyframes)
{
    cout << "[icp][optimization] finding suitable edges..." << endl;

    sm_result result;
    icp_edge edge;
    params->restart = 0;

    int nDisFilter, nAngFilter, nIterFilter, newEdges, tatal_count;
    nDisFilter = nAngFilter =  nIterFilter = newEdges = tatal_count = 0;
    int df = (int) params->kf_delta_frame / 2;
    double da, dx, dy, d;
    LDP kf1 = ld_alloc_new(keyframes[0].nrays);
    LDP kf2 = ld_alloc_new(keyframes[0].nrays);
    int num_frames = keyframes.size();
    for (int i=0; i<num_frames-df; i+=5) {
        kf1 = ld_from_keyframe(keyframes[i]);
        for (int j=i+df; j<num_frames; j+=df) {
            ++tatal_count;
            kf2 = ld_from_keyframe(keyframes[j]);

            /**
             * delta_distance < 2 meters to include
             * delta_angule < 30 degree to include
             * delta_frame > 15 frames to include
             * iteration < 30 times to include
             **/
//            da = kf2->true_pose[2] - kf1->true_pose[2];
            da = angleDiff(kf2->true_pose[2], kf1->true_pose[2]);
            if (abs(da) > 0.5236) {
                ++nAngFilter;
                ld_free(kf2);
                continue;
            }

            dx = abs(kf1->true_pose[0] - kf2->true_pose[0]);
            dy = abs(kf1->true_pose[1] - kf2->true_pose[1]);
            d = dx * dx + dy * dy;
            if (d > 4) {
                ++nDisFilter;
                ld_free(kf2);
                continue;
            }

            params->laser_ref = kf1;
            params->laser_sens = kf2;
            double odometry[3], ominus_laser[3], temp[3];
            pose_diff_d(kf2->odometry, kf1->odometry, odometry);
            ominus_d(params->laser, ominus_laser);
            oplus_d(ominus_laser, odometry, temp);
            oplus_d(temp, params->laser, params->first_guess);

            sm_icp(params, &result);
            if (result.iterations > 15) {
              ++nIterFilter;
              ld_free(kf2);
              continue;
            }

            /** finally a new edge born */
            ++newEdges;
            edge.edge_id_to = j;
            edge.edge_id_from = i;
            copy_d(result.x, 3, edge.transform);
            double w = compute_weight_new_edge(kf1->true_pose, kf2->true_pose, result.x, result.iterations);

            // use fixed weight
//            double w = 1.0;
            edge.information[0] = 20 * w;
            edge.information[1] = 5 * w;
            edge.information[2] = 286.4789 * w; // 900/pi, 286.4789
            edges.push_back(edge);

            ld_free(kf2);
        }
        ld_free(kf1);
    }
    // print edge selection results
//    printf("[icp] finding suitable edges done. tatal count: %d \n", tatal_count);
//    printf("[icp] discard edge matching through angle: %d(%f%%)\n",
//           nAngFilter, (double)100*nAngFilter/tatal_count);
//    printf("[icp] discard edge matching through distance: %d(%f%%)\n",
//           nDisFilter, (double)100*nDisFilter/(tatal_count-nAngFilter));
//    printf("[icp] discard edge matching through iterations: %d(%f%%)\n",
//           nIterFilter, (double)100*nIterFilter/(tatal_count-nDisFilter-nAngFilter));
    printf("[icp][optimization] Add tatal new edges: %d\n", newEdges);
}

/** for the need of new keyframe judgment */
bool newKeyframeNeeded(sm_params* param, double* trans)
{
    if (param->kf_delta_frame < 2)
        return true;
    if (param->kf_dist_linear == 0 || param->kf_dist_angular == 0)
        return true;
    if (fabs(trans[2]) > param->kf_dist_angular)
        return true;

    double x = trans[0];
    double y = trans[1];
    double kf_dist_linear_sq = param->kf_dist_linear * param->kf_dist_linear;
    if (x*x + y*y > kf_dist_linear_sq)
        return true;

    return false;
}

/** whether need optimization or not */
bool needOptimization(sm_params* params, vector<laser_data>& keyframes,
                      int bi, int li)
{
    if (li - bi < 3) return false;
    if (li - bi >= params->pg_max_frames) {
        cout << "[icp][optimization] need optimization for frames reach to " << li-bi+1 << endl;
        return true;
    }

    double dx = abs(keyframes[li].true_pose[0] - keyframes[bi].true_pose[0]);
    double dy = abs(keyframes[li].true_pose[1] - keyframes[bi].true_pose[1]);
    double d = dx * dx + dy * dy;
    if (d < 1) return false;
    else if (d > 3) {
        cout << "[icp][optimization] need optimization for d = " << d
             << " > 3m, get " << li-bi+1 << " frames to optimize." << endl;
        return true;
    }

    double da = angleDiff(keyframes[li].true_pose[2], keyframes[bi].true_pose[2]);
    if (abs(da) > 1.57) {
        cout << "[icp][optimization] need optimization for da = " << da
             << " > PI/2, index from " << bi << " to " << li
             << ", get " << li-bi+1 << " frames to optimize." << endl;
        return true;
    }

    return false;
}

/** compute weight for edge's information matrix */
double compute_weight_new_edge(const double gp_i[], const double gp_j[],
                               const double trans[], int iterations)
{
    double w = 1.0;
    if (iterations <= 5 ) return w;
    else {
        w = std::max( 1.0-(iterations-5)/50.0, 0.0 );
        if (iterations <= 10) return w;

//        double da =gp_j[2] - gp_i[2];
        double da = angleDiff(gp_j[2], gp_i[2]);
        double dx = gp_j[0] - gp_i[0] * cos(da) + gp_i[1] * sin(da);
        double dy = gp_j[1] - gp_i[0] * sin(da) - gp_i[1] * cos(da);
        double c = (dx - trans[0]) * (dx - trans[0]) +
                   (dy - trans[1]) * (dy - trans[1]) +
                   (da - trans[2]) * (da - trans[2]);
        c = 0.1 / (0.1 + c);
        w = 0.7 * w + 0.3 * c;

        return w;
    }
}

/** compute weight for edge's information matrix */
double compute_weight_neighbor_edge(int iterations)
{
    double w = 1.0;
    if (iterations > 5)
        w = std::max( 1.0-(iterations-5)/50.0, 0.0 );

    return w;
}

/** find extra edges to optimize the icp results */
void add_extra_edges(sm_params* params, vector<icp_edge>& edges,
                     vector<laser_data>& keyframes, int bi, int li)
{
    sm_result result;
    icp_edge edge;
    int nDisFilter, nAngFilter, nIterFilter, newEdges, tatal_count, npgframes;
    nDisFilter = nAngFilter =  nIterFilter = newEdges = tatal_count = 0;
    npgframes = li - bi + 1;

    /**
     * delta_angule < 20 degree to include
     * delta_distance < 0.2*kf_delta_frame meters to include
     * iterations < 25 times to include
     **/
    int nrays = keyframes[0].nrays;
    LDP kf1 = ld_alloc_new(nrays);
    LDP kf2 = ld_alloc_new(nrays);
    for (int i=0; i<npgframes-2; i++) {
        int j = i + 2;
        kf1 = ld_from_keyframe(keyframes[bi+i]);
        kf2 = ld_from_keyframe(keyframes[bi+j]);
        tatal_count++;

        double da = angleDiff(kf2->true_pose[2], kf1->true_pose[2]);
        if (abs(da) > 0.349) {
            ++nAngFilter;
            ld_free(kf1);
            ld_free(kf2);
            continue;
        }

        double dx = abs(kf2->true_pose[0] - kf1->true_pose[0]);
        double dy = abs(kf2->true_pose[1] - kf1->true_pose[1]);
        double d = dx * dx + dy * dy;
        if (d > 0.16) {
            ++nDisFilter;
            ld_free(kf1);
            ld_free(kf2);
            continue;
        }

        params->laser_ref = kf1;
        params->laser_sens = kf2;
        double odometry[3], ominus_laser[3], temp[3];
        pose_diff_d(kf2->odometry, kf1->odometry, odometry);
        ominus_d(params->laser, ominus_laser);
        oplus_d(ominus_laser, odometry, temp);
        oplus_d(temp, params->laser, params->first_guess);

        sm_icp(params, &result);
        if ( result.iterations > params->max_iterations/2 || !result.valid ) {
            ++nIterFilter;
            ld_free(kf1);
            ld_free(kf2);
            continue;
        }
        if ( !valid_transform(params, result.x) ){
            printf("[icp][optimization] drop keyframes: %d and %d\n", bi+i, bi+j);
            ++nDisFilter;
            ld_free(kf1);
            ld_free(kf2);
            continue;
        }

        /// finally a new edge born
        ++newEdges;
        edge.edge_id_from = i;
        edge.edge_id_to   = j;
        copy_d(result.x, 3, edge.transform);
        double w = compute_weight_new_edge(kf1->true_pose, kf2->true_pose,
                                           result.x, result.iterations);

        // use fixed weight
        edge.information[0] = 20 * w;
        edge.information[1] = 5 * w;
        edge.information[2] = 286.4789 * w; // 900/pi, 286.4789
        edges.push_back(edge);

        ld_free(kf2);
        ld_free(kf1);
    }
    printf("[icp][optimization] Add tatal new edges: %d / %d\n", newEdges, tatal_count);
//    printf("[icp] discard edge matching through distance: %d)\n", nDisFilter);
//    printf("[icp] discard edge matching through angle: %d\n", nAngFilter);
//    printf("[icp] discard edge matching through iterations: %d\n", nIterFilter);
}
