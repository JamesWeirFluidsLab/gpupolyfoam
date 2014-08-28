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
#include <exception>


#include "fvCFD.H"
#include "clockTimer.H"

#include "polyIdPairs.H"
#include "selectIdPairs.H"
#include "OpenMM.h"
using namespace OpenMM;

#ifdef WATER
#include "mdWater.H"
#include "molecularField.H"
#elif defined MONO
#include "mdAtomistic.H"
#else
#include "mdPoly.H"
#endif

//CONSTANT VARIABLES 
#define FEM2SEC 1e-15
#define E18 1e18
#define E12 1e-12
#define NM 1e9
#define NANO 1e-9
#define NM2RUF 1.66053886e-12 
#define DALTON 1.66053886e-27 //from KONSTANTINOS water code
#define CHARGE 1.602176487e-19 

using namespace std;


//- use typedef to declare clouds in common term 
#ifdef MONO
	typedef atomisticMoleculeCloud MOLECULE;
	typedef atomisticMolecule TypeMolecule;
	typedef selectIdPairs Pairs;
#elif defined WATER
	typedef molecularField MOLECULE;
	typedef waterMolecule TypeMolecule;
	typedef selectIdPairs Pairs;
#else
	typedef polyMoleculeCloud MOLECULE;
	typedef polyMolecule TypeMolecule;
	typedef polyIdPairs Pairs;
#endif
/**
 * structure to represent the classes into single block
 * so that it becomes easier to transfer
 */

enum STATES {
    Forces = 1,
    Positions = 2,
    Velocities = 3
};

struct poly_solver_t
{
    reducedUnits* redUnits;
    potential* pot;
#ifdef MONO
    atomisticMoleculeCloud* molecules; //if atomistic molecule
#elif defined WATER
    molecularField* molecules;
#else 
    polyMoleculeCloud* molecules;//if polyatomic molecule
#endif
    clockTimer* evolveTimer;
//declare OMM data structure

    clockTimer* ommTimer;
    clockTimer* openFoamTimer;
    System* system;
    Context* context;
    VelocityVerletIntegrator* integrator;
    Vec3 bBoxOMMinNm;
    double refTime, refMass, refLength, refForce, refCharge, refVelocity, deltaT,rCutInNM;
    Pairs* plid; //open foam polyIdPairs class for gpu solver
    int uniqueMolecules;
    bool isMolecular;
    poly_solver_t() : 
    redUnits(0), pot(0), molecules(0),
    evolveTimer(0),system(0), context(0), integrator(0), isMolecular(false){}
    ~poly_solver_t() {
        delete redUnits; delete pot; 
        delete molecules; delete evolveTimer;
        delete context; delete integrator; delete system;
    }
};

/**
 * get the scalling function for the current 
 * simulation
 */
const Foam::word getScallingFunction(const MOLECULE& mol);

#endif
