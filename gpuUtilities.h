/* 
 * File:   gpuUtilities.h
 * Author: saifmulla
 *
 * Created on 02 October 2013, 12:19
 */

#ifndef GPUUTILITIES_H
#define	GPUUTILITIES_H

//#include "gpuPolyFoam.h"

extern std::vector<Vec3> ommpositions, ommforces, ofpositions;

//declare utilities functions to be used 
/**
 * extract OF positions and generate and OMM equivalent 
 * array to be passed to OMM system
 */

int extractOFPostoOMM(std::vector<Vec3>& posinnm,struct poly_solver_t* sol);
/**
 * setOMMBox
 * extract information from OF to generate relevant data structures
 * such as reference units and box information 
 */
void setOMMBox(struct poly_solver_t* solver,const boundBox& bBoxOF,const double dt); 



#endif	/* GPUUTILITIES_H */

