
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <iostream>

// User Defined Header Files
#include "cudaDeviceBuffer.h"
#include "Parser.h"
#include "kernels.cuh"

// Constants
# define pi 3.14159265358979323846

// User Defined Functions
real_t RestDistance(const u_int nump,const real_t volumne) {

	real_t n = nump / volumne ;
	real_t a = std::cbrt(3.0 / (4 * pi * n )) ;
	return 0.893 * a ;
}

int main(int argc, char **argv){

	// Read from the file.
    Parser p(argv[1]);
    p.readParameters();
    p.readInputConfiguration();

    // Simulation Parameters
    real_t time_end = std::stod(p.params["time_end"]) ;
    real_t timestep_length = std::stod(p.params["timestep_length"]) ;
    u_int vtk_out_freq = std::stol(p.params["vtk_out_freq"]) ;
    u_int threads_per_blocks = std::stol(p.params["cl_workgroup_1dsize"]) ;
    std::string vtk_name = p.params["vtk_out_name_base"] ;
    real_t xmin = std::stod(p.params["x_min"]); // We will assume cubic domain
    real_t xmax = std::stod(p.params["x_max"]);
    real_t ymin = std::stod(p.params["y_min"]);
    real_t ymax = std::stod(p.params["y_max"]);
    real_t zmin = std::stod(p.params["z_min"]);
    real_t zmax = std::stod(p.params["z_max"]);
    real_t re 	= std::stod(p.params["re"]) ;

    // Computing the Cell length
    const real_t celllength = 2.0 * re / 3.0 ;
    const u_int numcellx = (xmax - xmin) / celllength ;

    // Number of particles
    const u_int numparticles = p.num_particles ;

    // Number of Cell
    const u_int numcells = numcellx * numcellx * numcellx ;

    // Buffers
    cudaDeviceBuffer<real_t> mass(numparticles,"Scalar") ;
    cudaDeviceBuffer<real_t> position(numparticles,"Vector") ;
    cudaDeviceBuffer<real_t> density(numparticles,"Scalar") ;
    cudaDeviceBuffer<real_t> pressure(numparticles,"Scalar") ;
    cudaDeviceBuffer<real_t> velocity(numparticles,"Vector") ;
    cudaDeviceBuffer<real_t> forceold(numparticles,"Vector") ;
    cudaDeviceBuffer<real_t> forcenew(numparticles,"Vector") ;
    cudaDeviceBuffer<int> cell_list(numcells,"Scalar");
    cudaDeviceBuffer<int> particle_list(numparticles,"Scalar");
    cudaDeviceBuffer<real_t> const_args(9,"Scalar");
    cudaDeviceBuffer<u_int> num_cells(3,"Scalar");

    // Fill the buffers from Initial list
    p.fillBuffers(mass,velocity,position) ;
	// Velocity according to Maxwell Boltzmann distribution

    //Filling in the host data for the constant arguments
    const_args[0] = xmin;
    const_args[1] = xmax;
    const_args[2] = ymin;
    const_args[3] = ymax;
    const_args[4] = zmin;
    const_args[5] = zmax;
    const_args[6] = celllength;
    const_args[7] = celllength;
    const_args[8] = celllength;

    //Number of cells per dimension assuming cubic domain
    num_cells[0] = numcellx;
    num_cells[1] = numcellx;
    num_cells[2] = numcellx;

    // Allocating memory on Device
    mass.allocateOnDevice();
    position.allocateOnDevice();
    velocity.allocateOnDevice();
    density.allocateOnDevice() ;
    forceold.allocateOnDevice();
    forcenew.allocateOnDevice();
    cell_list.allocateOnDevice();
    particle_list.allocateOnDevice();
    const_args.allocateOnDevice();
    num_cells.allocateOnDevice();

    //Copy to Device
    mass.copyToDevice();
    position.copyToDevice();
    velocity.copyToDevice();
    density.copyToDevice() ;
    forceold.copyToDevice();
    forcenew.copyToDevice();
    cell_list.copyToDevice();
    particle_list.copyToDevice();
    const_args.copyToDevice();
    num_cells.copyToDevice();

	// Algorithm to launch
    // Calculate the number of blocks to launch
    u_int numblocks  ;

	{
    	// Cell Discretization
		// Density Updation
		// Pressure Calculation
		// Force due to pressure
		// Force to viscosity
		// Force due to gravity

		// Boundary sweep for pressure and viscosity

	}
}


