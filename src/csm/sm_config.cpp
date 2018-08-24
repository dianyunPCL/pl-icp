#include "csm/sm_config.h"

void sm_config::readParams(std::string filename)
{
    std::ifstream fin(filename.c_str());
    if (!fin) {
        std::cerr << "[icp][error] Config file does not exist." << std::endl;
        return;
    }
    while(!fin.eof()) {
        std::string line;
        std::getline(fin, line);
        if (line[0] == '#') {
            // 以‘＃’开头的是注释
            continue;
        }
        int pos = line.find('#');
        if (pos != -1) {
            // 从井号到末尾的都是注释
            line = line.substr(0, pos);
        }

        // 去掉空格和回车
        pos = line.find(' ');
        while (pos != std::string::npos) {
            line.erase(pos, 1);
            pos = line.find(' ');
        }
        pos = line.find('\r');
        if (pos != std::string::npos)
            line.erase(pos, 1);

        // 参数值读取与存储
        pos = line.find("=");
        if (pos == -1)
            continue;
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos+1, line.length());

        params[key] = value;
    }

    std::cout << "[icp][info] Read config file done." << std::endl;
}


void set_params(sm_params* params, const sm_config* config)
{
    params->max_angular_correction_deg   = config->getParams<double>("max_angular_correction_deg");
    params->max_linear_correction        = config->getParams<double>("max_linear_correction");
    params->max_iterations               = config->getParams<int>("max_iterations");
    params->epsilon_xy                   = config->getParams<double>("epsilon_xy");
    params->epsilon_theta                = config->getParams<double>("epsilon_theta");
    params->max_correspondence_dist      = config->getParams<double>("max_correspondence_dist");
    params->sigma                        = config->getParams<double>("sigma");
    params->use_corr_tricks              = config->getParams<int>("use_corr_tricks");
    params->restart                      = config->getParams<int>("restart");
    params->restart_threshold_mean_error = config->getParams<double>("restart_threshold_mean_error");
    params->restart_dt                   = config->getParams<double>("restart_dt");
    params->restart_dtheta               = config->getParams<double>("restart_dtheta");
    params->clustering_threshold         = config->getParams<double>("clustering_threshold");
    params->orientation_neighbourhood    = config->getParams<int>("orientation_neighbourhood");
    params->use_point_to_line_distance   = config->getParams<int>("use_point_to_line_distance");
    params->do_alpha_test                = config->getParams<int>("do_alpha_test");
    params->do_alpha_test_thresholdDeg   = config->getParams<double>("do_alpha_test_thresholdDeg");
    params->outliers_maxPerc             = config->getParams<double>("outliers_maxPerc");
    params->outliers_adaptive_mult       = config->getParams<double>("outliers_adaptive_mult");
    params->outliers_adaptive_order      = config->getParams<double>("outliers_adaptive_order");
    params->do_visibility_test           = config->getParams<int>("do_visibility_test");
    params->outliers_remove_doubles      = config->getParams<int>("outliers_remove_doubles");
    params->do_compute_covariance        = config->getParams<int>("do_compute_covariance");
    params->debug_verify_tricks          = config->getParams<double>("debug_verify_tricks");
    params->min_reading                  = config->getParams<double>("min_reading");
    params->max_reading                  = config->getParams<double>("max_reading");
    params->use_ml_weights               = config->getParams<int>("use_ml_weights");
    params->use_sigma_weights            = config->getParams<int>("use_sigma_weights");
    params->laser[0]                     = config->getParams<double>("laser[0]");
    params->laser[1]                     = config->getParams<double>("laser[1]");
    params->laser[2]                     = config->getParams<double>("laser[2]");
    params->kf_delta_frame               = config->getParams<int>("kf_delta_frame");
    params->kf_dist_linear               = config->getParams<double>("kf_dist_linear");
    params->kf_dist_angular              = config->getParams<double>("kf_dist_angular");
    params->pg_max_iterations            = config->getParams<int>("pg_max_iterations");
    params->pg_max_frames                = config->getParams<int>("pg_max_frames");
}
