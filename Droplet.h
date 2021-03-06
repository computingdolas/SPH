/*
 * Droplet.h
 *
 *  Created on: 06-Jul-2016
 *      Author: sagardolas
 */

#ifndef DROPLET_H_
#define DROPLET_H_


#include "Type.h"

#include <string>
#include <fstream>
#include <iostream>
#include <random>

class Droplet {

private:

	real_t radius ;
	real_t x_center ;
	real_t y_center ;
	real_t z_center ;
	u_int bubble_particles ;
	u_int num_bubbles;
	u_int tank_height ;
	u_int tank_particles ;

public:

    Droplet(const real_t _radius , const real_t _xcenter,const real_t _ycenter,
            const real_t _zcenter, const u_int _bubbleparticles,
			const u_int _numbubbles,const u_int _tankheight,const u_int _tankparticles ) ;

	~Droplet() ;

    void MakeBubbleAndTank() ;

};

Droplet::Droplet(const real_t _radius , const real_t _xcenter,
			const real_t _ycenter,const real_t _zcenter ,const u_int _bubbleparticles,
			const u_int _numbubbles,const u_int _tankheight,const u_int _tankparticles ) : radius(_radius),
					x_center(_xcenter),y_center(_ycenter),z_center(_zcenter),bubble_particles(_bubbleparticles),num_bubbles(_numbubbles),
					tank_height(_tankheight),tank_particles(_tankparticles) {

}

Droplet::~Droplet(){

}

void Droplet::MakeBubbleAndTank(){

	std::ofstream out ;
	out.open("bubbleWithTank.dat") ;

    real_t mass  = 1.0;
    real_t vel = -10.0 ;
    real_t offset1 = 0.5 ;

    u_int numpxy = tank_particles / tank_height ;

    real_t dist = 0.5;

    std::cout<<"The distance between particles"<<dist<<std::endl  ;
    u_int numpx = std::sqrt(numpxy) ;
    std::cout<<"tank particles "<<tank_particles<<std::endl;
    out<<12500<<std::endl;

    u_int iter = 0 ; 
	// For the tank
    /*for (u_int varz = 0; varz < 5; ++varz) {
        for (u_int vary = 0; vary < 50; ++vary) {
            for (u_int varx = 0; varx < 50; ++varx) {
                out<<mass<<" "<<varx * dist + offset1<<" "<<vary * dist + offset1<<" "<<varz * dist + offset1<<" "<<"0"<<" "<<"0"<<" "<<"0"<<std::endl ;
                ++iter ; 
			}
		}
    }*/

    iter = 0;
    for (u_int varz = 0; varz < 5; ++varz) {
        for (u_int vary = 0; vary < 50; ++vary) {
            for (u_int varx = 0; varx < 50; ++varx) {
                out<<mass<<" "<<varx * dist + 10.0<<" "<<vary * dist + 10.0<<" "<<varz * dist + 10.0<<" "<<"-5"<<" "<<"0"<<" "<<"0"<<std::endl ;
                ++iter ;
            }
        }
    }


    /*
	iter = 0 ; 
	// for the bubble
	real_t a1,a2,a3,r,s;
    for (u_int drop = 0; drop < bubble_particles; ++drop) {
		do {

	        a1 = 2.0 * std::rand() / ((double)RAND_MAX + 1.0 ) - 1.0 ;
	        a2 = 2.0 * std::rand() / ((double)RAND_MAX + 1.0 ) - 1.0 ;
	        a3 = 2.0 * std::rand() / ((double)RAND_MAX + 1.0 ) - 1.0 ;

	        r = (a1 * a1) + (a2 * a2) + (a3 * a3) ;

	    } while (r >= 1.0);

		s = std::sqrt(-2.0 * log(r)/r) ;

		a1 = a1 * s * radius ;
		a2 = a2 * s * radius ;            
		a3 = a3 * s * radius ;
        out<<mass<<" "<<a1+x_center<<" "<<a2+ y_center<<" "<<a3+z_center<<" "<<"0"<<" "<<"0"<<" "<<-10<<std::endl ;
        ++iter ; 
        a1 = a2 = a3 = r = 0.0 ; 
	}
        */
	out.close();
}

#endif /* DROPLET_H_ */
