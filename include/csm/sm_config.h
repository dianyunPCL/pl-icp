#ifndef CONFIG_H
#define CONFIG_H

#include "csm/algos.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <boost/lexical_cast.hpp>


struct sm_config
{
    std::unordered_map<std::string, std::string> params;

    void readParams(std::string filename="../../cfg/sm_config.txt");

    template <class T>
    T getParams(const std::string& key) const
    {
        auto iter = params.find(key);
        if (iter == params.end()) {
            std::cerr << "[icp][error] Parameter name " << key << " not found!" << std::endl;
        }
        return boost::lexical_cast<T>(iter->second);
    }
};


void set_params(struct sm_params* params, const struct sm_config* config);


#endif // PARAMETER_READER_H
