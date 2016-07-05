/*
 * kernels.cuh
 *
 *  Created on: 03-Jul-2016
 *      Author: sagardolas
 */

#ifndef KERNELS_CUH_
#define KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include "Type.h"

// Constants
# define pi 3.14159265358979323846

// Function to find out the thread id
__device__ u_int GetGlobalId(){

	return blockIdx.x *blockDim.x + threadIdx.x;
}

// Function to find out the global cell index
__device__ u_int GetGlobalCellId(const int i,const int j, const int k, const u_int cellx){

	return i + j * cellx + cellx * cellx * k ;

}

// Initialize the cell list .. cell parallel
__global__ void InitializeCellList(int *cell_list,const u_int num_cell){

	u_int cdx = GetGlobalId() ;
	cell_list[cdx] = -1 ;
}

// Calculate the norm
__device__ real_t norm(const real_t *r){

    real_t sum = 0.0 ;
    for (int i =0 ; i < 3 ; ++i){
        sum  += r[i] * r[i] ;
    }
    return sqrt(sum) ;

}

__device__ real_t W(const real_t re, const real_t r ) {

	const real_t c = ( 315.0 /(64.0 * pi * pow(re,9))) ;
	const real_t rij = pow(pow(re,2) - pow(r,2),3) ;
	return  (c * rij) ;

}

// DeltaW for pressure
__device__ real_t deltaW(const real_t re, const real_t r){

	const real_t c = ( 45.0/(pi * pow(re,6) * r)) ;
	const real_t rij = pow(re - r,3) ;
	return (c * rij) ;

}

// Delta function for viscosity
__device__ real_t deltaWvis(const real_t re, const real_t r){

	const real_t c = ( 45.0/(pi * pow(re,6))) ;
	const real_t rij = (re - r) ;
	return c * rij ;
}

// Initialize the Particle List
__global__ void InitializePartList(int *Particle_list, const u_int numparticles){

	u_int pdx = GetGlobalId();
	Particle_list[pdx] = pdx ;
}

// Update the list
__global__ void UpdateList(int * cell_list,
							int * particle_list ,
							const real_t *position,
							const real_t cell_length,
							const u_int num_particles,
							const u_int num_cell){

	u_int id = GetGlobalId() ;

	if (id < num_particles){

		// Find the index of the particles
		u_int pid = id * 3;

		// Find the cordinates of the particles
		real_t pos[3] = {position[pid],position[pid+1],position[pid+2]} ;

		// Find the 3D index of the cell in which it lies in
		u_int i = pos[0] / cell_length ;
		u_int j = pos[1] / cell_length ;
		u_int k = pos[2] / cell_length ;

		// Find the global cell id
		u_int cell_index  = GetGlobalCellId(i,j,k,num_cell) ;

		//Register the particle in linked list
		int old = atomicExch(&cell_list[cell_index],id) ;
		particle_list[id] = old ;
	}
}

// Find Neighbour list
__device__ void findNeighbour(const u_int xindex,
							const u_int yindex ,
							const u_int zindex,
							const u_int numcellx,
							u_int *neighbour_list){

	int tempx,tempy, tempz ;

	u_int count = 0 ;
	for (int i = -1	; i <=1; ++i) {
		for(int j =-1 ; j <=2 ; ++j ){
			for (int k = -1 ; k <=2 ; ++k){
				if (i != 0 || j != 0 || k !=0) {
					tempx 	=  (int)xindex + i ;
					tempy 	=  (int)yindex + j ;
					tempz 	=  (int)zindex + k ;

					 if(!(((tempx <0) || (tempx>static_cast<int>(numcellx-1))  ) || ( ( tempy<0 ) ||
							 ( tempy>static_cast<int>(numcellx-1)) ) ||
							 ( (tempz<0) || (tempz> static_cast<int>(numcellx-1) ) )) ){

						 neighbour_list[count] = GetGlobalCellId(tempx,tempy,tempz,numcellx) ;
						 ++count ;
					 }
				}
			}
		}
	}
	neighbour_list[count] = GetGlobalCellId(xindex, yindex, zindex,numcellx);
}

// Calculate the density
__global__ void CalculateDensity(const real_t *mass,
								const int *cell_list,
								const int *particle_list,
								real_t *density,
								const real_t *position,
								const real_t re,
								const u_int num_particles,
								const real_t cell_length,
								const u_int numcellx ){


	u_int pid = GetGlobalId() ;

	if (pid  < num_particles) {

		density[pid] = 0.0 ;
		u_int vpid = pid * 3 ;
		real_t pos[3] = {position[vpid],position[vpid+1],position[vpid+2]} ;

		// Find out in which cell it lies
		int cell_index[3] = {0} ;
		cell_index[0] =	pos[0]/cell_length ;
		cell_index[1] =	pos[1]/cell_length ;
		cell_index[2] = pos[2]/cell_length ;
		int cell_id = GetGlobalCellId(cell_index[0],cell_index[1],cell_index[2],numcellx) ;

		u_int nl[26] = {0} ;
		findNeighbour(cell_index[0],cell_index[1],cell_index[2],numcellx,nl) ;

		// Traversing the own cell
		for (int head_id = cell_list[cell_id]; head_id !=-1	; head_id = particle_list[head_id]) {
			if (head_id != pid){

				// find out the index of the particle
				u_int nvpid = head_id * 3 ;

				// find out the position of the neighbour particle
				real_t npos[3] = {position[nvpid],position[nvpid+1],position[nvpid+2]} ;

				// find out the relative vector
				real_t relvec[3] = {0.0,0.0,0.0} ;
				for (int var = 0; var < 3; ++var) {
					relvec[var] = npos[var]- pos[var] ;
				}

				// Find out the norm of relative vector
				real_t rnorm = norm(relvec) ;

				// Add contribution to the density
				density[pid] += W(re,rnorm) * mass[pid] ;

			}
		}

		// Traversing the NeighbourHood cell
		u_int temp = 0 ;
		for (u_int n = nl[temp]; n!= cell_id;++temp){
			for (int head_id = cell_list[n]; head_id !=-1; head_id = particle_list[head_id]) {
				if (head_id != pid){

					// find out the index of the particle
					u_int nvpid = head_id * 3 ;

					// find out the position of the neighbour particle
					real_t npos[3] = {position[nvpid],position[nvpid+1],position[nvpid+2]} ;

					// find out the relative vector
					real_t relvec[3] = {0.0} ;
					for (int var = 0; var < 3; ++var) {
						relvec[var] = npos[var]- pos[var] ;
					}

					// Find out the norm of relative vector
					const real_t r = norm(relvec) ;

					// Add contribution to the density
					density[pid] += mass[pid]*W(re,r) ;
				}
			}
		}
	}
}

//Calculate the pressure for each particle
__global__ void CalculatePressure(real_t *pressure ,
								const real_t *density,
								const real_t restpressure,
								const real_t restdensity,
								const real_t k,
								const u_int num_particles){

	u_int pid = GetGlobalId() ;

	if (pid < num_particles ){
		pressure[pid] = restpressure + k*(density[pid] - restdensity) ;
	}
}

// Calculate the Force due to pressure and viscosity
__global__ void CalculateForce(const real_t *velocity,
								real_t * force,
								const int *cell_list,
								const int *particle_list ,
								const real_t *mass,
								const real_t *pressure ,
								const real_t *density,
								const real_t *position,
								const u_int num_particles,
								const real_t cell_length,
								const u_int numcellx,
								const real_t re,
								const real_t nu){

	u_int pid = GetGlobalId() ;

	if (pid  < num_particles) {

		u_int vpid = pid * 3 ;
		force[vpid] = 0.0 ;
		force[vpid+1] = 0.0 ;
		force[vpid+2] = 0.0 ;

		real_t pos[3] = {0} ;
		pos[0] =  position[vpid] ;
		pos[1] =  position[vpid+1] ;
		pos[2] =  position[vpid+2] ;

		// Find out in which cell it lies
		int cell_index[3] = {pos[0]/cell_length,pos[1]/cell_length,pos[2]/cell_length };
		int cell_id = GetGlobalCellId(cell_index[0],cell_index[1],cell_index[2],numcellx) ;

		u_int nl[26] = {0} ;
		findNeighbour(cell_index[0],cell_index[1],cell_index[2],numcellx,nl) ;

		// Traversing the own cell
		for (int head_id = cell_list[cell_id]; head_id !=-1	; head_id = particle_list[head_id]) {
			if (head_id != pid){

				// find out the index of the particle
				u_int nvpid = head_id * 3 ;

				// find out the position of the neighbour particle
				real_t npos[3] = {position[nvpid],position[nvpid+1],position[nvpid+2]} ;

				// find out the relative vector
				real_t relvec[3] = {0.0,0.0,0.0} ;
				for (int var = 0; var < 3; ++var) {
					relvec[var] = npos[var]- pos[var] ;
				}

				// Find out the norm of relative vector
				real_t rnorm = norm(relvec) ;

				// Add contribution to the forces due to pressure
				real_t constant =  mass[head_id] * ( (pressure[pid] + pressure[head_id])/ (2 * density[head_id])) * deltaW(re,rnorm) ;
				force[vpid] += -constant * relvec[0] ;
				force[vpid+1] += -constant * relvec[1] ;
				force[vpid+2] += -constant * relvec[2] ;

				// Add contribution due to viscosity
				real_t constantvis = nu * mass[head_id] * deltaWvis(re,rnorm) / density[head_id] ;
				force[vpid] += constantvis * (velocity[nvpid] - velocity[vpid]) ;
				force[vpid+1] += constantvis * (velocity[nvpid+1] - velocity[vpid+1]) ;
				force[vpid+2] += constantvis * (velocity[nvpid+2] - velocity[vpid+2]) ;

			}
		}

		// Traversing the NeighbourHood cell
		u_int temp = 0 ;
		for (u_int n = nl[temp]; n!= cell_id;++temp){
			for (int head_id = cell_list[n]; head_id !=-1; head_id = particle_list[head_id]) {
				if (head_id != pid){

					// find out the index of the particle
					u_int nvpid = head_id * 3 ;

					// find out the position of the neighbour particle
					real_t npos[3] = {position[nvpid],position[nvpid+1],position[nvpid+2]} ;

					// find out the relative vector
					real_t relvec[3] = {0.0} ;
					for (int var = 0; var < 3; ++var) {
						relvec[var] = npos[var]- pos[var] ;
					}

					// Find out the norm of relative vector
					const real_t r = norm(relvec) ;

					// Add contribution to the forces
					real_t constant =  mass[head_id] * ( (pressure[pid] + pressure[head_id]) / (2 * density[head_id])) * deltaW(re,rnorm) ;
					force[vpid] += -constant * relvec[0] ;
					force[vpid+1] += -constant * relvec[1] ;
					force[vpid+2] += -constant * relvec[2] ;

					// Add contribution due to viscosity
					real_t constantvis = nu * mass[head_id] * deltaWvis(re,rnorm) / density[head_id] ;
					force[vpid] += constantvis * (velocity[nvpid] - velocity[vpid]) ;
					force[vpid+1] += constantvis * (velocity[nvpid+1] - velocity[vpid+1]) ;
					force[vpid+2] += constantvis * (velocity[nvpid+2] - velocity[vpid+2]) ;

				}
			}
		}
	}
}

// Boundary Sweep for Pressure and Density
__global__ void BoundarySweep() {

	// Traverse all the cell and add the effects of pressure near
	// boundary and density
	// due to boundary particles


}

// Updating velocity ...
__global__ void velocityUpdate(const real_t *forceOld,
								const real_t *forceNew,
								const real_t *mass,
								real_t *velocity,
								const u_int numparticles,
								const real_t deltat ){

	u_int pid = GetGlobalId() ;

	if (pid < numparticles){
		u_int vpid = pid * 3 ;
	    velocity[vpid] += ( (forceNew[vpid] + forceOld[vpid]) * deltat ) / (2.0 * mass[pid] ) ;
	    velocity[vpid+1] += ( (forceNew[vpid+1] + forceOld[vpid+1]) * deltat ) / (2.0 * mass[pid] ) ;
	    velocity[vpid+2] += ( (forceNew[vpid+2] + forceOld[vpid+2]) * deltat ) / (2.0 * mass[pid] ) ;

	}

}

// Updating the position
__global__ void positionUpdate(const real_t *force,
								real_t *position,
								const real_t *velocity,
								const real_t *mass,
								const u_int num_particles,
								const real_t deltat
								) {


	u_int pid = GetGlobalId() ;

	if (pid < num_particles){

		u_int vpid = pid * 3 ;
        position[vpid]   += (deltat * velocity[vpid] ) + ( (force[vpid] *  deltat * deltat) / ( 2.0 * mass[vpid]) ) ;
        position[vpid+1] += (deltat * velocity[vpid+1] ) + ( (force[vpid+1] * deltat * deltat) / ( 2.0 * mass[vpid]) ) ;
        position[vpid+2] += (deltat * velocity[vpid+2] ) + ( (force[vpid+2] * deltat * deltat) / ( 2.0 * mass[vpid]) ) ;
	}
}

// Finding the boundary cells and the index of it
#endif /* KERNELS_CUH_ */
