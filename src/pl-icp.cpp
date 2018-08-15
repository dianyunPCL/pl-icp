#include <time.h>
#include <string.h>
#include <cstring>
#include <libgen.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <dirent.h>
#include <ctime>
#include <sstream>

#include "options/options.h"
#include "csm/csm_all.h"
#include "csm/utils.h"
#include "csm/laser_data_yogo.h"
#include "icp/pl_icp.h"
#include "icp/icp_optimization.h"

using namespace std;

struct {
    const char * file_in;
    const char * file_out;
    int format;

    /* which algorithm to run */
    int algo;

    int recover_from_error;
    int delta_frame;

    int debug;
} p;

extern void sm_options(struct sm_params*p, struct option*ops);

ofstream fout;
string g2o_file = "080920-5-kf.g2o";

int get_num_files(const char* path);
void print_result(const sm_params& params, const sm_result& result);
void write_trajectory(const sm_params& params, const sm_result& result, const string& file_name);
void write_to_g2o_file(const vector<icp_vertex>& vertexs, const vector<icp_edge>& edges);
bool newKeyframeNeeded(sm_params* param, double* trans);

int main(int argc, const char*argv[]) {
    sm_set_program_name(argv[0]);

    struct sm_params params;
    struct sm_result result;

    vector<struct icp_vertex> vertexs;
    vector<struct icp_edge> edges;
    vector<laser_data> key_frames;

    struct icp_edge edge;
    struct icp_vertex vertex;

    /** set params */
    struct option* ops = options_allocate(100);
    options_string(ops, "in", &p.file_in, "stdin", "Input file ");
    options_string(ops, "out", &p.file_out, "stdout", "Output file ");
    options_int(ops, "algo", &p.algo, 0, "Which algorithm to use (0:(pl)ICP 1:gpm-stripped 2:HSM) ");
    options_int(ops, "debug", &p.debug, 0, "Shows debug information");
    options_int(ops, "recover_from_error", &p.recover_from_error, 0, "If true, tries to recover from an ICP matching error");
    options_int(ops, "format", &p.format, 0,"Output format (0: log in JSON format, 1: log in Carmen format (not implemented))");
    sm_options(&params, ops);
    if(!options_parse_args(ops, argc, argv)) {
        fprintf(stderr, "\n\nUsage:\n");
        options_print_help(ops, stderr);
        return -1;
    }

    /** clear output file when exist */
    if (p.file_out != NULL) {
        fout.open(p.file_out);
        fout.close();
    }

    /** Generate scan data from filename */
    int file_num = get_num_files(p.file_in);
    vertexs.reserve(file_num);
    edges.reserve(4*file_num);
    printf("[pl-icp] find tatal %d files in dataset %s\n", file_num, p.file_in);

    std::vector<std::string> files;
    for (int i=0; i<file_num; ++i) {
      std::stringstream ss;
      ss << p.file_in << i << ".txt";
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

    /** For the first scan, set estimate = odometry
     *  Add first vertex and first keyframe
     */
    vertex.vertex_id = 0;
    vertex.global_pose[0] = 0.;
    vertex.global_pose[1] = 0.;
    vertex.global_pose[2] = 0.;
    vertexs.push_back(vertex);

    copy_d(vertex.global_pose, 3, laser_ref->global_pose);
    copy_d(laser_ref->odometry, 3, laser_ref->estimate);

    key_frames.reserve(file_num);   // distribute space
    key_frames.push_back(*laser_ref);

    /** read data and run icp */
    printf("[pl-icp] dealting with laser_data ...\n");
    // vars for keyframe
    bool need_new_frame = true;
    int skip_frame = 0;
    int count = 0;

    LDP laser_sens = ld_alloc_new(nrays);
    clock_t time_stt = clock();
    for ( int i=1; i<file_num; ++i) {
        // run between delta frames
        if ( !need_new_frame && skip_frame<5) {
            ++skip_frame;
            continue;
        }
        ++count;
        skip_frame = 0;

        laser_sens = ld_from_yogo_stream(files[i].c_str(), nrays, count);
        if(!ld_valid_fields(laser_sens))  {
            sm_error("Invalid laser data in (#%d in file).\n", i);
            return -3;
        }

        params.laser_ref  = laser_ref;
        params.laser_sens = laser_sens;

        /* Set first guess as the difference in odometry */
        double odometry[3], ominus_laser[3], temp[3];
        pose_diff_d(laser_sens->odometry, laser_ref->odometry, odometry);
        ominus_d(params.laser, ominus_laser);
        oplus_d(ominus_laser, odometry, temp);
        oplus_d(temp, params.laser, params.first_guess);

        /* Do the actual work */
        sm_icp(&params, &result);


        // check the result
        if(!result.valid){
            if(p.recover_from_error) {
                sm_info("One ICP matching failed. Because you passed  -recover_from_error, I will try to recover."
                " Note, however, that this might not be good in some cases. \n");
                sm_info("The recover is that the displacement is set to 0. No result stats is output. \n");

                /* For the first scan, set estimate = odometry */
                copy_d(laser_ref->estimate, 3, laser_sens->estimate);
                ld_free(laser_ref);
                laser_ref = laser_sens;
            } else {
                sm_error("One ICP matching failed. Because I process recursively, I will stop here.\n");
                sm_error("Use the option -recover_from_error if you want to try to recover.\n");
                ld_free(laser_ref);
                ld_free(laser_sens);
                return 2;
            }
        }
        else {   // PL-ICP成功，获取结果
            /* Add the result to the previous estimate */
            oplus_d(laser_ref->estimate, result.x, laser_sens->estimate);

            double transform[3] = {result.x[0], result.x[1], result.x[2]};
            vector<double> gp = get_global_pose(laser_ref->global_pose, result.x);
            laser_sens->global_pose[0] = gp[0];
            laser_sens->global_pose[1] = gp[1];
            laser_sens->global_pose[2] = gp[2];
            key_frames.push_back(*laser_sens);

            // add vertex and edge
            vertex.vertex_id = count;   // i
            copy_d(laser_sens->global_pose, 3, vertex.global_pose);
            vertexs.push_back(vertex);

            edge.edge_id_behind = count;
            edge.edge_id_ahead = count - 1;
            copy_d(transform, 3, edge.transform);
            edge.information[0] = 25;
            edge.information[1] = 50;
            edge.information[2] = 286.4789; // 900/pi
            edges.push_back(edge);

            if (p.file_out != NULL) {
                write_trajectory(params, result, p.file_out);
            }

            // free useless data
            ld_free(laser_ref);
            laser_ref = laser_sens;

            // update need_new_frame flag
            if (params.kf_dist_linear ==0 || params.kf_dist_angular == 0)
                need_new_frame = true;
            else
                need_new_frame = newKeyframeNeeded(&params, transform);
        }
    }
    ld_free(laser_sens);
    double t = (clock() - time_stt) / (float)CLOCKS_PER_SEC;
    printf("[pl-icp] done. Tatal cost time: %fs, average: %fms\n", t, 1000*t/count);

    // add edges
    printf("[pl-icp] start finding new edges...\n");
    time_stt = clock();
    find_new_edges(params, edges, key_frames);
    t = (clock() - time_stt) / (float)CLOCKS_PER_SEC;
    printf("[pl-icp] find edges done. Cost time: %fs\n", t);

    // write SE2
    write_to_g2o_file(vertexs, edges);

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
void write_trajectory(const sm_params& params, const sm_result& result, const string& file_name)
{
    try {
        fout.open(file_name.c_str(), ios::app);
        fout << params.laser_sens->global_pose[0] << " "
             << params.laser_sens->global_pose[1] << " "
             << params.laser_sens->global_pose[2] << " ";

        if (params.do_compute_covariance) {
            for (int i=0; i<3; ++i)
                for (int j=0; j<3; ++j)
                    fout << *gsl_matrix_ptr(result.cov_x_m, i, j) << " ";
        }
        fout << "\n";
        fout.close();
    } catch(...) {
        printf("[pl-icp] write to result file failed!\n");
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
            fout << "VERTEX_SE2 " << v.vertex_id << " " << v.global_pose[0] << " "
                 << v.global_pose[1] << " " << v.global_pose[2] << "\n";
        }
        for (auto& e : edges) {
            fout << "EDGE_SE2 " << e.edge_id_behind << " " << e.edge_id_ahead << " "
                 << e.transform[0] << " " << e.transform[1] << " " << e.transform[2] << " "
                 << e.information[0] << " 0 0 " << e.information[1] << " 0 "
                 << e.information[2] << "\n";
        }
        fout.close();
    } catch(...) {
        printf("[pl-icp] write to g2o file failed!\n");
    }
    printf("[pl-icp] write to g2o file sucessed: %s\n", g2o_file.c_str());
}

/** for the need of new keyframe judgment */
bool newKeyframeNeeded(sm_params* param, double* trans)
{
  if (fabs(trans[2] > param->kf_dist_angular))
      return true;

  double x = trans[0];
  double y = trans[1];
  double kf_dist_linear_sq = param->kf_dist_linear * param->kf_dist_linear;
  if (x*x + y*y > kf_dist_linear_sq)
      return true;

  return false;
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

