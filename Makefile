#Compiler FLags 
CFLAGS = -c -std=c++11 -O3 
CCg = g++ 

#Compiler flags for Cuda 
CFLAGScuda = -std=c++11 -arch sm_20 -O3 
CC = nvcc

all	:sph

sph	:Parser.o VTKWriter.o
		$(CC) $(CFLAGScuda) Parser.o VTKWriter.o Simulation.cu -o sph

Parser.o	:Parser.cpp
		$(CCg) $(CFLAGS) Parser.cpp

VTKWriter.o	:VTKWriter.cpp
		$(CCg) $(CFLAGS) VTKWriter.cpp

clean		:
		rm -rf *.o *.out  
