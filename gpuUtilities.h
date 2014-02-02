/* 
 * File:   gpuUtilities.h
 * Author: saifmulla
 *
 * Created on 02 October 2013, 12:19
 */

#ifndef GPUUTILITIES_H
#define	GPUUTILITIES_H

//#include "gpuPolyFoam.h"
#include <../cuda-5.0/include/CL/cl_platform.h>

extern std::vector<Vec3> ommpositions, ommforces, ofpositions;

/**
 * initialize the openmm system
 * this function creates list of number of particles in the syste
 * by retriving relevant information from OF and assigns to OMM
 * system, additionally initialises
 * context
 * system
 * integrator
 */
void initialiseOMM(struct poly_solver_t* solver);

/**
 * extract sigma and eps values
 * to create a string of lennardJones parameters
 * which could be used as argument to custom non bonded force function
 */
void extractCoeffParameters(const MOLECULE& mol,
				const polyIdPairs& plid,
				std::vector<std::string>& coeffStr);

/**
 * extract OF particles
 * and correspondingly add to the OMM system
 */
void extractOFParticles(struct poly_solver_t* solver, 
                        CustomNonbondedForce* const nonbonded);
/**
 * extract OF positions and generate and OMM equivalent 
 * array to be passed to OMM system
 */

int extractOFPostoOMM(std::vector<Vec3>& posinnm,
                      const struct poly_solver_t* sol,
                      const boundBox& bb);

/**
 * extract OF velocities and generate OMM equivalent array
 * to be passed to OMM system
 */
int extractOFVeltoOMM(std::vector<Vec3>& velinnm, const struct poly_solver_t* sol, int particles);
/**
 * setOMMBox
 * extract information from OF to generate relevant data structures
 * such as reference units and box information 
 */
void setOMMBox(struct poly_solver_t* solver,const boundBox& bBoxOF,const double dt);

/**
 * add particle to custom non bonded force
 * by retriving from potential dict list
 */
void addParticlesToNonBonded(CustomNonbondedForce* const nonbonded,
                                const struct poly_solver_t* solver);

/**
 * set obtained and calculate forces from omm to 
 * openfoam
 */

int setOFforce(struct poly_solver_t* solver, const std::vector<Vec3>& ommForces);

/**
 * getOMMState
 * get information from openmm
 */
void getOMMState(const Context* context,std::vector<Vec3>& statearray);

/**
 * extractOFQ
 * this function extracts Q tensor for each molecule and assigns it to 
 * OpenMM array with conversion
 */
int extractOFQ(const struct poly_solver_t* solver, std::vector<OpenMM::Tensor>& moleculeQ);
/**
 * extractOFSiteRefPositions
 * this function extracts site reference positions from each molecule of OF system and pushes 
 * to vector of Vec3 type which will be equal to number of molecules
 */
int extractOFSiteRefPositions(const struct poly_solver_t* solver, std::vector<OpenMM::Vec3>& siteRefPositions);
/**
 * extractMoleculePositions
 * this function extracts Positions of a molecule from each molecule of OF system and pushes
 * to vector of Vec3 type which will be equal to number of molecules
 */
int extractMoleculePositions(const struct poly_solver_t* solver, std::vector<OpenMM::Vec3>& moleculePositions);

/**
 * extractMoleculePI
 * this function extracts PI values from each molecule of OF system and pushes to OMM
 * vector of Vec3 type which will be equal to number of molecules
 */
int extractMoleculePI(const struct poly_solver_t* solver, std::vector<OpenMM::Vec3>& moleculePI);

/**
 * setOFpositions
 */
void setOFPositions(struct poly_solver_t* solver, const std::vector<Vec3>& posInNm);

/**
 * setOFVelocities
 */
void setOFVelocities(struct poly_solver_t* solver, const std::vector<Vec3>& velInNm);
//                enum STATES st);
#endif	/* GPUUTILITIES_H */

