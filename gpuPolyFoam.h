/**
 * //identity: $Id$
 * filename: gpuPolyFoam.h
 * @developer: Saif Mulla
 * @created: 21/09/2013
 * Application: gpuPolyFoam
 * Description:
 */
#ifndef GPUPOLYFOAM_H
#define	GPUPOLYFOAM_H

#include <stdio.h>
#include <iostream>


#include "fvCFD.H"
#include "clockTimer.H"
#include "selectIdPairs.H"
#ifdef USE_OMM
#include "OpenMM.h"
using namespace OpenMM;
#endif

//#ifdef WATER
//#include "mdWater.H"
//#include "molecularField.H"
#ifdef MONO
#include "mdAtomistic.H"
#else
#include "mdPoly.H"
#endif

//CONSTANT VARIABLES 
#define FEM2SEC 1e-15
#define NANSEC 1e-9

using namespace std;


//- use typedef to declare clouds in common term 
#ifdef MONO
	typedef atomisticMoleculeCloud MOLECULE;
#else
	typedef polyMoleculeCloud MOLECULE;
#endif
/**
 * structure to represent the classes into single block
 * so that it becomes easier to transfer
 */

struct poly_solver_t
{
	reducedUnits* redUnits;
	potential* pot;

#ifdef MONO
	atomisticMoleculeCloud* molecules; //if atomistic molecule
#else
	polyMoleculeCloud* molecules;//if polyatomic molecule
#endif
	clockTimer* evolveTimer;
//declare OMM data structure
#ifdef USE_OMM
	System* system;
	Context* context;
	Integrator* integrator;
	Vec3 bBoxOMMinNm;
	double refTime, refMass, refLength, refForce, refCharge, deltaT;
        selectIdPairs* sid; //open foam
#endif
	poly_solver_t() : 
	redUnits(0), pot(0), molecules(0),
#ifdef USE_OMM
	evolveTimer(0),system(0), context(0), integrator(0){}
#else
	evolveTimer(0) {}
#endif
	~poly_solver_t() {
		delete redUnits; delete pot; 
		delete molecules; delete evolveTimer;
#ifdef USE_OMM
		printf("Solver deleted \n");
		delete context; delete integrator; delete system;
#endif
	}
};

#endif
