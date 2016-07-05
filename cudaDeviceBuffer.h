/*
 * cudaDeviceBuffer.h
 *
 *  Created on: 30-Jun-2016
 *      Author: sagardolas
 */

#ifndef CUDADEVICEBUFFER_H_
#define CUDADEVICEBUFFER_H_

#include <stdio.h>
#include <vector>
#include <cuda_runtime.h>
#include <iostream>
#include <string>

#include "Type.h"



#define DIM 3  //..dimension of the system we ar working on ..//
template <typename type>
class cudaDeviceBuffer {

private:

    //PhysicalQuantity phyVar ;
    std::vector<type> data ; // Data on the host ...//
    u_int actualSize ; // Actual size to be allocated ...//
    u_int numBytes_ ;

public:


    type *devicePtr ; // Pointer on the Device ..///

    // Constructing the Buffer......//
    cudaDeviceBuffer(u_int numParticles_,const std::string physicalquantity) ;
    ~cudaDeviceBuffer() ;

    // Accessing the data.....//
    const type& operator[](u_int index_) const ;
    type& operator[](u_int index_) ;

    // Memory Operations....//
    void copyToHost() ;
    void copyToDevice() ;

    // Allocate and Deallocate Memory on device....//
    void allocateOnDevice() ;
    void freeMemoryOnDevice() ;

    //Calculate Memory...//
    void bytesToAllocate() ;

    //Reseting the data array....//
    void reset() ;

    // Check Error....//
    void checkError(const cudaError_t err) ;
};

template<typename type>
cudaDeviceBuffer<type>::cudaDeviceBuffer(u_int _numParticles,const std::string physicalquantity ){

    if (physicalquantity == "Scalar" ) {
    	actualSize = _numParticles ;
    }
    else
        actualSize = _numParticles * DIM;

    // Allocating the data
    data.resize(actualSize,0.0) ; // May be problem here

    // Calculating the Number of Bytes
    bytesToAllocate() ;

}
template <typename type>
cudaDeviceBuffer<type>::~cudaDeviceBuffer<type>(){

    freeMemoryOnDevice() ; // Do not call this function untill and unless you are sure that you have allocated memory in the device

}
template<typename type>
const type& cudaDeviceBuffer<type>::operator[](u_int index_) const{

    return data[index_] ;
}

template<typename type>
type& cudaDeviceBuffer<type>::operator[](u_int index_){

    return data[index_] ;
}

template<typename type>
void cudaDeviceBuffer<type>::copyToDevice(){

    checkError(cudaMemcpy(devicePtr, &data[0], numBytes_, cudaMemcpyHostToDevice)) ;
}
template<typename type >
void cudaDeviceBuffer<type>::copyToHost(){

    checkError(cudaMemcpy(&data[0], devicePtr, numBytes_, cudaMemcpyDeviceToHost)) ;
}

template<typename type>
void cudaDeviceBuffer<type>::allocateOnDevice(){

    checkError(cudaMalloc(&devicePtr, numBytes_));
}

template<typename type>
void cudaDeviceBuffer<type>::freeMemoryOnDevice(){

    checkError(cudaFree(devicePtr)) ;
}

template<typename type>
void cudaDeviceBuffer<type>::bytesToAllocate(){

    //calculate number of bytes
    numBytes_ = actualSize * sizeof(type) ;
}

template<typename type>
void cudaDeviceBuffer<type>::reset(){

    // Not defined yet
}

template<typename type>
void cudaDeviceBuffer<type>::checkError(const cudaError_t err){

    if(err!= cudaSuccess){
        std::cout<<cudaGetErrorString(err)<<std::endl ;
        exit(-1) ;
    }
}
#endif /* CUDADEVICEBUFFER_H_ */
