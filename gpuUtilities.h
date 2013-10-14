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

/**
 * initialize the openmm system
 * this function creates list of number of particles in the syste
 * by retriving relevant information from OF and assigns to OMM
 * system, additionally initialises
 * context
 * system
 * integrator
 */
int initialiseOMM(struct poly_solver_t* solver);

/**
 * extract sigma and eps values
 * to create a string of lennardJones parameters
 * which could be used as argument to custom non bonded force function
 */
void extractLennardJonesParameters(const MOLECULE& mol,
				const polyIdPairs& plid,
				std::string& epsString,
				std::string& sigmaString);

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

int setOFforce(struct poly_solver_t* solver, const std::vector<Vec3>& ommForce);
#endif	/* GPUUTILITIES_H */

