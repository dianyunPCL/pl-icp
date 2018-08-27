#include "csm/csm_all.h"
#include "csm/utils.h"
#include "csm/laser_data_yogo.h"
#include "icp/icp_yogo.h"
#include "icp/icp_optimization.h"

#include <iomanip>
#include <string.h>
#include <cstring>
#include <libgen.h>
#include <vector>
#include <sstream>
#include <dirent.h>
#include <ctime>

using namespace std;

ofstream fout;
string g2o_file = "080920-5-kf.g2o";
string data_path = "/data/laser_data_yogo/laser080920/laser5/";
string config_file = "../cfg/sm_config.txt";

int get_num_files(const char* path);
void print_result(const sm_params& params, const sm_result& result);
void write_trajectory(const vector<vector<double>> poses, const string& file_name);
void write_to_g2o_file(const vector<icp_vertex>& vertexs, const vector<icp_edge>& edges);


int main(int argc, const char*argv[]) {
    sm_set_program_name(argv[0]);

//    if (argc <= 2) {
//        cerr << "[icp] Usage: pl-icp data_path g2o_file" << endl;
//        return -1;
//    }
    data_path = argv[1];
//    g2o_file = argv[2];
    string trajectory_file = argv[2];

    struct sm_config config;
    struct sm_params params;
    struct sm_result result;

    vector<vector<double>> pose_raw;
    vector<double> pose_r(3);
    vector<laser_data> keyframes;
    vector<icp_vertex> vertexs;
    vector<icp_edge> edges;
    struct icp_edge edge;
    struct icp_vertex vertex;

    /** set params */
    config.readParams(config_file);
    set_params(&params, &config);

    /** Generate scan data from filename */
    int nfiles = get_num_files(data_path.c_str());
    vertexs.reserve(nfiles);
    edges.reserve(4*nfiles);
    printf("[icp][info] find tatal %d files in dataset %s\n", nfiles, data_path.c_str());

    std::vector<std::string> files;
    for (int i=0; i<nfiles; ++i) {
      std::stringstream ss;
      ss << data_path << i << ".txt";
      files.push_back(ss.str());
    }

    /** Read first scan */
    int nrays = 937;
    LDP laser_ref = ld_alloc_new(nrays);
    if ( !(laser_ref = ld_from_yogo_stream(files[0].c_str(), nrays, 0)) ) {
        sm_error("Could not read first scan.\n");
        return -1;
    }
    if(!ld_valid_fields(laser_ref))  {
        sm_error("Invalid laser data in first scan.\n");
        return -2;
    }
    double first_pose[] = {0., 0., 0.};
    copy_d(first_pose, 3, laser_ref->true_pose);
    copy_d(first_pose, 3, laser_ref->last_trans);
    copy_d(laser_ref->odometry, 3, laser_ref->estimate);
    // @Vance: 在这里计算首帧点的笛卡尔坐标，后面icp函数里就不用重复计算了
    ld_compute_cartesian(laser_ref);

    /** For the first scan, set estimate = odometry
     *  Add first vertex and first keyframe
     */
    vertex.vertex_id = 0;
    vertex.vertex_pose[0] = 0.0;
    vertex.vertex_pose[1] = 0.0;
    vertex.vertex_pose[2] = 0.0;
    vertexs.push_back(vertex);

    keyframes.reserve(nfiles);   // distribute space
    keyframes.push_back(*laser_ref);

    /** read data and run icp */
    cout << "[icp][info] dealting with laser_data ..." << endl;
    // vars for keyframe
    bool need_new_frame = true;
    int skip_frame = 0, drop_frames = 0;
    int kf_delta_frame = params.kf_delta_frame;
    int frame_id = 0;

    LDP laser_sens = ld_alloc_new(nrays);
    std::clock_t time_stt = clock();
    for ( int i=1; i<nfiles; ++i) {
        // run between delta frames
        if ( !need_new_frame && skip_frame < kf_delta_frame) {
            skip_frame++;
            continue;
        }
        frame_id++;
        skip_frame = 0;

        laser_sens = ld_from_yogo_stream(files[i].c_str(), nrays, frame_id);
        if(!ld_valid_fields(laser_sens))  {
            sm_error("Invalid laser data in (#%d in file).\n", i);
            return -3;
        }

        params.laser_ref  = laser_ref;
        params.laser_sens = laser_sens;

        /* Set first guess as the difference in odometry */
        double odometry[3];
        pose_diff_d(laser_sens->odometry, laser_ref->odometry, odometry);
        copy_d(odometry, 3, params.first_guess);


        /* Do the actual work */
        sm_icp(&params, &result);

        // check the result
        if(!result.valid){
            printf("[icp][error] ICP result unvalid! ICP stop!\n");
            ld_free(laser_ref);
            ld_free(laser_sens);
            return -2;
        }
        else {   // PL-ICP成功，获取结果
            /* Add the result to the previous estimate */
//            if(!valid_transform(result.x)) {
//                frame_id--;
//                continue;
//            }
            valid_transform(&params, result.x);
            oplus_d(laser_ref->estimate, result.x, laser_sens->estimate);
            copy_d(result.x, 3, laser_sens->last_trans);
            double gp[3];
            get_global_pose(gp, laser_ref->true_pose, result.x);
            copy_d(gp, 3, laser_sens->true_pose);
            keyframes.push_back(*laser_sens);
//            cout << std::setprecision(6);
//            cout << "[icp][info] #" << frame_id << "# valid transform: " << result.x[0]
//                 << ", " << result.x[1] << ", " << result.x[2]
//                 << ", iterations = " << result.iterations << endl;
//            cout << "[icp][info] #" << frame_id << "# global pose: " << gp[0]
//                 << ", " << gp[1] << ", " << gp[2] << endl;

//            if (frame_id > 497 && frame_id < 502) {
//                printf("-- #%d iteration: %d, nvalid: %d, valid_transform: %f, %f, %f, true pose: %f, %f, %f\n",
//                       frame_id, result.iterations, result.nvalid,
//                       result.x[0], result.x[1], result.x[2], gp[0], gp[1], gp[2]);
//            }

            // add vertex and edge
            vertex.vertex_id = frame_id;   // i
            copy_d(laser_sens->true_pose, 3, vertex.vertex_pose);
            vertexs.push_back(vertex);

            edge.edge_id_from = frame_id - 1;
            edge.edge_id_to = frame_id;
            copy_d(result.x, 3, edge.transform);
            edge.information[0] = 25;
            edge.information[1] = 5;
            edge.information[2] = 286.4789; // 900/pi
            edges.push_back(edge);

            // free useless data
            ld_free(laser_ref);
            laser_ref = laser_sens;

            // update need_new_frame flag
            need_new_frame = newKeyframeNeeded(&params, result.x);
        }
    }
    ld_free(laser_ref);

    double t = (clock() - time_stt) / (float)CLOCKS_PER_SEC;
    cout << "[icp][info] done. Tatal cost time: " << t << "s, average: "
         << 1000*t/frame_id << "ms" << endl;

    // add edges
//    time_stt = clock();
//    find_new_edges(&params, edges, keyframes);
//    t = (clock() - time_stt) / (float)CLOCKS_PER_SEC;
//    printf("[icp] find edges done. Cost time: %fs\n", t);

    // write SE2
//    write_to_g2o_file(vertexs, edges);

    for (int i=0; i<nfiles; ++i) {
        double gp[3];
        copy_d(keyframes[i].true_pose, 3, gp);
        pose_r = {gp[0], gp[1], gp[2]};
        pose_raw.push_back(pose_r);
    }
    write_trajectory(pose_raw, trajectory_file);


    return 0;
}


/** print icp results */
void print_result(const sm_params& params, const sm_result& result)
{    
    printf("2) Matching Result: [%f, %f, %f]\n", result.x[0], result.x[1], result.x[2]);
    printf("3) Total Iterations: %d\n", result.iterations);
    printf("4) Valid Correspondences: %d (%f%%)\n", result.nvalid, result.nvalid/9.37);
    printf("5) Total Correspondence Error: %f\n", result.error);
    if (params.do_compute_covariance) {
        printf("6) Covariance Matrix:\n");
        std::cout << *result.cov_x_m << std::endl;
    }
}

/** write the icp trajectory to a file */

void write_trajectory(const vector<vector<double>> poses, const string& file_name)
{
    try {
        fout.open(file_name.c_str());
        fout.close();
        fout.open(file_name.c_str(), ios::app);
        for (int i=0; i<poses.size(); ++i) {
            fout << poses[i][0] << " "
                 << poses[i][1]  << " "
                 << poses[i][2]  << "\n";
        }
        fout.close();
    } catch(...) {
        printf("[icp][error] write to result file failed!\n");
    }
}


/** write the vertexs and edges to a g2o file */
void write_to_g2o_file(const vector<icp_vertex>& vertexs, const vector<icp_edge>& edges)
{
    try {
        fout.open(g2o_file.c_str());
        fout.close();
        fout.open(g2o_file.c_str(), ios::app);
        for (auto& v : vertexs) {
            fout << "VERTEX_SE2 " << v.vertex_id << " " << v.vertex_pose[0] << " "
                 << v.vertex_pose[1] << " " << v.vertex_pose[2] << "\n";
        }
        for (auto& e : edges) {
            fout << "EDGE_SE2 " << e.edge_id_from << " " << e.edge_id_to << " "
                 << e.transform[0] << " " << e.transform[1] << " " << e.transform[2] << " "
                 << e.information[0] << " 0 0 " << e.information[1] << " 0 "
                 << e.information[2] << "\n";
        }
        fout.close();
    } catch(...) {
        printf("[icp][error] write to g2o file failed!\n");
    }
    printf("[icp][info] write to g2o file sucessed: %s\n", g2o_file.c_str());
}


/** for the file number calculation in the dataset */
static int filterDot(const struct dirent * dir) {
  if (strcmp(dir->d_name,".") == 0 || strcmp(dir->d_name, "..") == 0) {
    // 过滤掉 "."和".."
    return 0;
  } else {
    return 1;
  }
}
int get_num_files(const char* path)
{
    struct dirent **namelist;

    int file_nums = scandir(path, &namelist, filterDot, alphasort);

    free(namelist);
    return file_nums;
}

