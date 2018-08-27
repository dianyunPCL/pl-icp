#include "csm/csm_all.h"
#include "csm/utils.h"
#include "csm/laser_data_yogo.h"
#include "icp/icp_yogo.h"
#include "icp/icp_optimization.h"

#include <g2o/core/block_solver.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/core/factory.h>
#include <g2o/core/optimization_algorithm_factory.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <g2o/types/slam2d/vertex_se2.h>
#include <g2o/types/slam2d/edge_se2.h>

#include <iomanip>
#include <string.h>
#include <cstring>
#include <libgen.h>
#include <vector>
#include <sstream>
#include <dirent.h>
#include <ctime>
#include <thread>
#include <mutex>

using namespace std;
using namespace g2o;

typedef BlockSolver<BlockSolverTraits<-1,-1> > SlamBlockSolver;
typedef LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;

mutex mu;
ofstream fout;

string data_path = "/data/laser_data_yogo/laser080920/laser5/";
string config_file = "../cfg/sm_config.txt";


int get_num_files(const char* path);
void write_trajectory(const vector<vector<double>>& poses, const string& file_name);
void pose_graph(sm_params* params, vector<laser_data>& key_frames,
                vector<icp_edge> edges, int begin_index);
void thread_pose_graph(sm_params params, vector<laser_data>& keyframes,
                       vector<icp_edge> edges, int begin_index, int last_index);

int main(int argc, const char*argv[]) {
    sm_set_program_name(argv[0]);

//    if (argc <= 2) {
//        cerr << "[icp] Usage: pl-icp data_path g2o_file" << endl;
//        return -1;
//    }
    data_path = argv[1];
    string opt_file = argv[2];

    vector<laser_data> keyframes;
    vector<vector<double>> pose_output;
    vector<double> pose_out(3);
    vector<icp_edge> edges;
    struct icp_edge edge;

    struct sm_params params;
    struct sm_result result;
    struct sm_config config;

    /** set params */
    config.readParams(config_file);
    set_params(&params, &config);

    /** Generate scan data from filename */
    int nfiles = get_num_files(data_path.c_str());
    cout << "[icp] find tatal files: " << nfiles << endl;
    vector<string> files;
    for (int i=0; i<nfiles; ++i) {
        stringstream ss;
        ss << data_path << i << ".txt";
        files.push_back(ss.str());
    }

    /** Read first scan */
    int nrays = 937;
    LDP laser_ref = ld_alloc_new(nrays);
    LDP laser_sens = ld_alloc_new(nrays);
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

    /// vars for keyframe  
    pose_output.reserve(nfiles);
    keyframes.reserve(nfiles);
    keyframes.push_back(*laser_ref);    // add first frame
    bool need_new_frame = true;
    int frame_id = 0;
    int vertex_id = 0;
    int skip_frame = 0;
    int kf_delta_frame = params.kf_delta_frame;

    /// vars for optimization
    int optimize_time = 0;
    double information[] = {20, 5, 286.4789};
    double tatal_opt_time = 0.0;
    bool need_optimization = false;
    int opt_begin_index = 0, opt_last_index = 0;

    /** read data and run icp */
    cout << "[icp] dealting with laser_data ..." << endl;
    std::clock_t time_stt = clock();
    for ( int i=1; i<nfiles; ++i) {
        /** run between delta frames */
        if ( !need_new_frame && skip_frame < kf_delta_frame) {
            skip_frame++;
            continue;
        }
        frame_id++;
        skip_frame = 0;

        laser_sens = ld_from_yogo_stream(files[i].c_str(), nrays, frame_id);
        if(!ld_valid_fields(laser_sens))  {
            sm_error("[icp][] Invalid laser data in (#%d in file).\n", i);
            return -3;
        }

        vertex_id++;
        params.laser_ref = laser_ref;
        params.laser_sens = laser_sens;

        /** Set first guess as the difference in odometry */
        double odometry[3];
        pose_diff_d(laser_sens->odometry, laser_ref->odometry, odometry);
        copy_d(odometry, 3, params.first_guess);

        /** Do the actual work */
        sm_icp(&params, &result);

        /** check the result */
        if(!result.valid){
            printf("[icp] ICP result unvalid! ICP stop!\n");
            ld_free(laser_ref);
            ld_free(laser_sens);
            return -2;
        }
        else {   // PL-ICP成功，获取结果
            /* Add the result to the previous estimate */
            valid_transform(&params, result.x);
            oplus_d(laser_ref->true_pose, result.x, laser_sens->estimate);
            copy_d(result.x, 3, laser_sens->last_trans);
            double gp[3];
            get_global_pose(gp, laser_ref->true_pose, result.x);
            copy_d(gp, 3, laser_sens->true_pose);
            keyframes.push_back(*laser_sens);
//            cout << setprecision(6);
//            cout << "[icp][info] #" << frame_id << "# valid transform: " << result.x[0]
//                 << ", " << result.x[1] << ", " << result.x[2]
//                 << ", iterations = " << result.iterations << endl;
//            cout << "[icp][info] #" << frame_id << "# global pose: " << gp[0]
//                 << ", " << gp[1] << ", " << gp[2] << endl;

            // add vertexs and edges
            edge.edge_id_from = vertex_id - 1;
            edge.edge_id_to = vertex_id;
            copy_d(result.x, 3, edge.transform);
            double w = compute_weight_neighbor_edge(result.iterations);
            double inf[3];
            inf[3] = information[0] * w;
            inf[3] = information[1] * w;
            inf[3] = information[2] * w;
            copy_d(inf, 3, edge.information);
            edges.push_back(edge);

            // free useless data
            ld_free(laser_ref);
            laser_ref = laser_sens;

            // update need_new_frame flag
            need_new_frame = newKeyframeNeeded(&params, result.x);
            opt_last_index = frame_id;
        }

        /** begin optimization */
        need_optimization = needOptimization(&params, keyframes, opt_begin_index, opt_last_index);
        if (need_optimization) {
//            int opt_last_index = keyframes.size() - 1;
            if (opt_last_index != frame_id) { cout << "last_index != frame_id!" << endl; }

            std::clock_t op_start = clock();
//            pose_graph(&params, keyframes, edges, begin_index);
            std::thread thr(thread_pose_graph, params, std::ref(keyframes),
                            std::move(edges), opt_begin_index, opt_last_index);
//            thr.detach();
            thr.join();
            double t = (clock() - op_start) / (double)CLOCKS_PER_SEC;
            cout << "[icp][optimization] optimization # " << optimize_time
                 << " cost time: " << 1000*t << "ms" << endl;
            tatal_opt_time += 1000*t;

            // update laser_ref true pose after optimizatio
            copy_d(keyframes[opt_last_index].true_pose, 3, laser_ref->true_pose);
            vertex_id = 0;
            optimize_time++;
            opt_begin_index = opt_last_index;
            cout << "[icp][info] #" << frame_id << "# global pose up to: " << laser_ref->true_pose[0]
                 << ", " << laser_ref->true_pose[1] << ", " << laser_ref->true_pose[2] << endl;
        }
    }
    ld_free(laser_sens);

    double tt = (clock() - time_stt) / (double)CLOCKS_PER_SEC;
    cout << "[icp][info] pl-icp done. Tatal cost time: " << tt << "s, average: "
         << 1000*tt/frame_id << "ms" << endl;
    cout << "[icp][[optimization]] Optimization tatal cost time: " << tatal_opt_time
         << "ms, average: " << tatal_opt_time/optimize_time << "ms" << endl;


    for (int i=0; i<nfiles; ++i) {
        double gp[3];
        copy_d(keyframes[i].true_pose, 3, gp);
        pose_out = {gp[0], gp[1], gp[2]};
        pose_output.push_back(pose_out);
    }
    printf("[icp][info] write trajectories..\n");
    write_trajectory(pose_output, opt_file);

    return 0;
}


/** for optimization update **/
void pose_graph(sm_params* params, vector<laser_data>& keyframes,
                vector<icp_edge> edges, int begin_index)
{
//    int nframs = params->pg_max_frames;
    int nframs = keyframes.size() - begin_index;
    int li = nframs;

    SlamLinearSolver* linearSolver = new SlamLinearSolver();
    linearSolver->setBlockOrdering(false);
    SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
//    OptimizationAlgorithmGaussNewton* optimizationAlgorithm =
//         new OptimizationAlgorithmGaussNewton(blockSolver);
    OptimizationAlgorithmLevenberg* optimizationAlgorithm =
            new OptimizationAlgorithmLevenberg(blockSolver);

    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(optimizationAlgorithm);
    optimizer.setVerbose(false);

    /// add vertexs
    int nrays = keyframes[0].nrays;
    LDP laser = ld_alloc_new(nrays);
    laser = ld_from_keyframe(keyframes[begin_index]);
    g2o::VertexSE2* pose = new g2o::VertexSE2;
    pose->setId(0);
    pose->setEstimate(Eigen::Vector3d(laser->true_pose[0], laser->true_pose[1], laser->true_pose[2]));
    pose->setFixed(true);
    optimizer.addVertex(std::move(pose));
    ld_free(laser);

    for (size_t i=1; i<nframs; ++i) {
        laser = ld_from_keyframe(keyframes[begin_index+i]);
        g2o::VertexSE2* pose = new g2o::VertexSE2;
        pose->setId(i);
        pose->setEstimate(Eigen::Vector3d(laser->true_pose[0],
                          laser->true_pose[1], laser->true_pose[2]));
        optimizer.addVertex(std::move(pose));

        ld_free(laser);
    }

    /// generate new edges
    add_extra_edges(params, edges, keyframes, begin_index, li);

    /// add edges
    Eigen::Matrix3d martrixEdgetestbai = Eigen::Matrix3d::Identity();
    for (auto& e : edges) {
        martrixEdgetestbai(0,0) = e.information[0];
        martrixEdgetestbai(1,1) = e.information[1];
        martrixEdgetestbai(2,2) = e.information[2];

        g2o::EdgeSE2* icp_edge = new g2o::EdgeSE2;
        icp_edge->vertices()[0] = optimizer.vertex(e.edge_id_from);
        icp_edge->vertices()[1] = optimizer.vertex(e.edge_id_to);
        icp_edge->setMeasurement(Eigen::Vector3d(e.transform[0],
                                 e.transform[1], e.transform[2]));
        icp_edge->setInformation(martrixEdgetestbai);

        optimizer.addEdge(std::move(icp_edge));
    }

    optimizer.initializeOptimization(0);
    optimizer.optimize(params->pg_max_iterations);

    // update pose
    Eigen::Vector3d es;
    for(size_t i = 0; i < nframs; ++i) {
        g2o::VertexSE2* vSE2 = static_cast<g2o::VertexSE2*>(optimizer.vertex(i));
        es = vSE2->estimate().toVector();
        keyframes[i+begin_index].true_pose[0] = es[0];
        keyframes[i+begin_index].true_pose[1] = es[1];
        keyframes[i+begin_index].true_pose[2] = es[2];
    }

    optimizer.clear();
    Factory::destroy();
    OptimizationAlgorithmFactory::destroy();
    HyperGraphActionLibrary::destroy();
}


void write_trajectory(const vector<vector<double>>& poses, const string& file_name)
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
        printf("[icp][info] write to result file sucessed!\n");
    } catch(...) {
        printf("[icp][error] write to result file failed!\n");
    }
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

void thread_pose_graph(sm_params params, vector<laser_data>& keyframes,
                       vector<icp_edge> edges, int bi, int li)
{
//    mu.lock();
    int nframes = params.pg_max_frames;
//    int li = keyframes.size() - 1;
//    cout << "b and e:" << bi << ", " << li << endl;

    SlamLinearSolver* linearSolver = new SlamLinearSolver();
    linearSolver->setBlockOrdering(false);
    SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
    OptimizationAlgorithmLevenberg* optimizationAlgorithm =
            new OptimizationAlgorithmLevenberg(blockSolver);

    g2o::SparseOptimizer optimizer;
    optimizer.setAlgorithm(optimizationAlgorithm);
    optimizer.setVerbose(false);

    /// add vertexs
    int nrays = keyframes[0].nrays;
    LDP laser = ld_alloc_new(nrays);
    laser = ld_from_keyframe(keyframes[bi]);
    g2o::VertexSE2* pose = new g2o::VertexSE2;
    pose->setId(0);
    pose->setEstimate(Eigen::Vector3d(laser->true_pose[0], laser->true_pose[1], laser->true_pose[2]));
    pose->setFixed(true);
    optimizer.addVertex(std::move(pose));
    ld_free(laser);
    for (size_t i=1; i<nframes; ++i) {
        laser = ld_from_keyframe(keyframes[bi+i]);
        g2o::VertexSE2* pose = new g2o::VertexSE2;
        pose->setId(i);
        pose->setEstimate(Eigen::Vector3d(laser->true_pose[0],
                          laser->true_pose[1], laser->true_pose[2]));
        optimizer.addVertex(std::move(pose));

        ld_free(laser);
    }

    /// generate new edges
    add_extra_edges(&params, edges, keyframes, bi, li);

    /// add edges
    Eigen::Matrix3d martrixEdgetestbai = Eigen::Matrix3d::Identity();
    for (auto& e : edges) {
        martrixEdgetestbai(0,0) = e.information[0];
        martrixEdgetestbai(1,1) = e.information[1];
        martrixEdgetestbai(2,2) = e.information[2];

        g2o::EdgeSE2* icp_edge = new g2o::EdgeSE2;
        icp_edge->vertices()[0] = optimizer.vertex(e.edge_id_from);
        icp_edge->vertices()[1] = optimizer.vertex(e.edge_id_to);
        icp_edge->setMeasurement(Eigen::Vector3d(e.transform[0],
                                 e.transform[1], e.transform[2]));
        icp_edge->setInformation(martrixEdgetestbai);

        optimizer.addEdge(std::move(icp_edge));
    }

    optimizer.initializeOptimization();
    optimizer.optimize(params.pg_max_iterations);

    // update pose
    Eigen::Vector3d es;
    for(size_t i = 0; i < nframes; ++i) {
        g2o::VertexSE2* vSE2 = static_cast<g2o::VertexSE2*>(optimizer.vertex(i));
        es = vSE2->estimate().toVector();
        keyframes[i+bi].true_pose[0] = es[0];
        keyframes[i+bi].true_pose[1] = es[1];
        keyframes[i+bi].true_pose[2] = es[2];
    }

    optimizer.clear();
    Factory::destroy();
    OptimizationAlgorithmFactory::destroy();
    HyperGraphActionLibrary::destroy();

//    mu.unlock();

}


