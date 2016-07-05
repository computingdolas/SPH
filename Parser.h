/*
 * Parser.h
 *
 *  Created on: 03-Jul-2016
 *      Author: sagardolas
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <string>
#include <map>
#include <iostream>
#include <vector>
#include <fstream>

#include "Type.h"
#include "cudaDeviceBuffer.h"

class Parser{
public:

	int num_params;
    int num_particles;
    std::string filename;
    std::map<std::string,std::string> params;

    std::vector<real_t> mass;
    std::vector<real_t> pos;
    std::vector<real_t> vel;

    //Constructor
    Parser(std::string filename){
        this->num_params = 0;
        this->filename = filename;
    }

    //Parse the parameters
    void readParameters();

    //Read input configuration
    void readInputConfiguration();

    // Fill the cudaDeviceBuffers
    void fillBuffers(cudaDeviceBuffer<real_t> &mass, cudaDeviceBuffer<real_t> &velocity, cudaDeviceBuffer<real_t> &position ) ;

};

#endif /* PARSER_H_ */
