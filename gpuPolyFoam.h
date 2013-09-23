/**
 * //identity: $Id$
 * filename: gpuPolyFoam.h
 * @developer: Saif Mulla
 * @created: 21/09/2013
 * Application: gpuPolyFoam
 * Description:
 */

#include <stdio.h>
#include <iostream>


#include "fvCFD.H"
#include "clockTimer.H"
#include "mdAtomistic.H"
#include "selectIdPairs.H"
#ifdef USE_OMM
#include "OpenMM.h"
using namespace OpenMM;
#endif

#ifdef WATER
#include "mdWater.H"
#include "molecularField.H"
#elif defined POLY
#include "mdPoly.H"
#endif

//CONSTANT VARIABLES 
#define FEM2SEC 1e-15

using namespace std;


//- use typedef to declare clouds in common term 
#ifdef WATER
	typedef molecularField MOLECULE;
#elif defined MONO
	typedef atomisticMoleculeCloud MOLECULE;
#elif defined POLY
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
#ifdef WATER
	molecularField* molecules; //if water molecule
#elif defined MONO
	atomisticMoleculeCloud* molecules; //if atomistic molecule
#elif defined POLY
	polyMoleculeCloud* molecules;//if polyatomic molecule
#endif
	clockTimer* evolveTimer;
	selectIdPairs* sid;
//declare OMM data structure
#ifdef USE_OMM
	System* system;
	Context* context;
	Integrator* integrator;
	Vec3 bBoxOMMinNm;
	double refTime, refMass, refLength, refForce, deltaT;
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
		delete context; delete integrator; delete system;
#endif
	}
};

extern std::vector<Vec3>& ommpositions, ommforces, ofpositions;
//declare utilities functions to be used 
/**
 * extract OF positions and generate and OMM equivalent 
 * array to be passed to OMM system
 */
int extractOFPostoOMM(std::vector<Vec3>& posinnm,struct poly_solver_t* sol);

