
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <iostream>
#include <random>

// User Defined Header Files
#include "cudaDeviceBuffer.h"
#include "Parser.h"
#include "kernels.cuh"
#include "VTKWriter.h"

// Constants
# define pi 3.14159265358979323846

// User Defined Functions
real_t RestDistance(const u_int nump,const real_t volumne) {

	real_t n = nump / volumne ;
	real_t a = std::cbrt(3.0 / (4 * pi * n )) ;
	return 0.893 * a ;
}

void MaxWellBoltzmannVel(cudaDeviceBuffer<real_t> &vel){

	u_int num_particles = vel.size()/ 3  ;
	real_t a1,a2,a3,r,s;

	for (u_int p = 0; p < num_particles; ++p) {
		u_int pindex = p * 3 ;

	    do {

	        a1 = 2.0 * std::rand() / ((double)RAND_MAX + 1.0 ) - 1.0 ;
	        a2 = 2.0 * std::rand() / ((double)RAND_MAX + 1.0 ) - 1.0 ;
	        a3 = 2.0 * std::rand() / ((double)RAND_MAX + 1.0 ) - 1.0 ;

	        r = (a1 * a1) + (a2 * a2) + (a3 * a3) ;

	    } while (r >= 1.0);

	    s = std::sqrt(-2.0 * log(r)/r) ;

	    vel[pindex] = a1 *s ;
	    vel[pindex+1] = a2 *s ;
	    vel[pindex+2] = a3 *s ;
	}
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
    real_t nu = std::stod(p.params["nu"]);
    real_t k = std::stod(p.params["k"]);
    real_t restpressure = std::stod(p.params["p0"]);
    real_t restdensity = std::stod(p.params["rho"]);

    // Computing the Cell length
    const real_t celllength = re  ;
    const u_int numcellx = (xmax - xmin) / celllength ;

    // Number of particles
    const u_int numparticles = p.num_particles ;
    real_t d = RestDistance(numparticles,xmax*ymax*zmax);

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

    // Algorithm to launch real_t d = RestDistance(numparticles,xmax*ymax*zmax);
    // Calculate the number of blocks to launch
    u_int blocks_p,blocks_c,threads_p,threads_c;
    threads_p = threads_per_blocks;
    threads_c = threads_per_blocks;

    if(numparticles%threads_per_blocks == 0){
        blocks_p = numparticles/threads_p;
    }
    else{
        blocks_p  = numparticles/threads_p+1;
    }


    if(numcells%threads_per_blocks){
        blocks_c = numcells/threads_c;
    }
    else{
        blocks_c  = numparticles/threads_c+1;
    }

    VTKWriter writer(vtk_name);
	{

        InitializeCellList<<<blocks_c,threads_c>>>(cell_list.devicePtr,numcells) ;

		InitializePartList<<<blocks_p,threads_p>>>(particle_list.devicePtr,numparticles) ;

		UpdateList<<<blocks_p,threads_p>>>(cell_list.devicePtr,particle_list.devicePtr,
							position.devicePtr,celllength,numparticles,numcellx) ;

        //CalculateDensity<<<blocks_p,threads_p >>>(mass.devicePtr,cell_list.devicePtr,particle_list.devicePtr,
                        //density.devicePtr,position.devicePtr,re,numparticles,celllength,
                        //numcellx) ;

        //CalculatePressure<<<blocks_p,threads_p>>>(pressure.devicePtr,density.devicePtr,
            //			restpressure,restdensity,k,numparticles) ;

        //CalculateForce<<<blocks_p,threads_p>>>(velocity.devicePtr,forcenew.devicePtr,cell_list.devicePtr,
                //	particle_list.devicePtr,mass.devicePtr,pressure.devicePtr,
                //	density.devicePtr,position.devicePtr,numparticles,celllength,numcellx,re,nu) ;

		BoundarySweep<<<blocks_p,threads_p>>>   (forcenew.devicePtr,density.devicePtr,mass.devicePtr,timestep_length,position.devicePtr,d,numparticles,re,const_args[0],1);

        int iter=0;
        for (real_t t = 0.0;t<=time_end && iter==0; t+= timestep_length) {
            if(iter % vtk_out_freq == 0){
                // copy to host back
                forcenew.copyToHost();
                forceold.copyToHost();
                position.copyToHost();
                velocity.copyToHost();
                writer.writeVTKOutput(mass,position,velocity,numparticles);
            }
			positionUpdate<<<blocks_p,threads_p>>>(forcenew.devicePtr,position.devicePtr,velocity.devicePtr,
					mass.devicePtr,numparticles,timestep_length) ;

            copyForces<<<blocks_p,threads_p>>>(forceold.devicePtr,forcenew.devicePtr,numparticles);

            InitializeCellList<<<blocks_c,threads_c>>>(cell_list.devicePtr,numcells) ;

			InitializePartList<<<blocks_p,threads_p>>>(particle_list.devicePtr,numparticles) ;

			UpdateList<<<blocks_p,threads_p>>>(cell_list.devicePtr,particle_list.devicePtr,
								position.devicePtr,celllength,numparticles,numcellx) ;

            //CalculateDensity<<<blocks_p,threads_p >>>(mass.devicePtr,cell_list.devicePtr,particle_list.devicePtr,
                            //density.devicePtr,position.devicePtr,re,numparticles,celllength,
                        //	numcellx) ;

            //CalculatePressure<<<blocks_p,threads_p>>>(pressure.devicePtr,density.devicePtr,
                //			restpressure,restdensity,k,numparticles) ;

            //CalculateForce<<<blocks_p,threads_p>>>(velocity.devicePtr,forcenew.devicePtr,cell_list.devicePtr,
                    //	particle_list.devicePtr,mass.devicePtr,pressure.devicePtr,
                        //density.devicePtr,position.devicePtr,numparticles,celllength,numcellx,re,nu) ;
			BoundarySweep<<<blocks_p,threads_p>>>   (forcenew.devicePtr,density.devicePtr,mass.devicePtr,timestep_length,
						position.devicePtr,d,numparticles,re,const_args[0],1);

			velocityUpdate<<<blocks_p,threads_p>>>(forceold.devicePtr,forcenew.devicePtr,mass.devicePtr,
					velocity.devicePtr,numparticles,timestep_length);

            iter++;

		}
	}
}


