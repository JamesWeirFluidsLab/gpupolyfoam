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
#include "mdPoly.H"
#include "clockTimer.H"


/**
 * structure to represent the classes into single block
 * so that it becomes easier to transfer
 */

struct poly_solver_t
{
	reducedUnits* redUnits;
	potential* pot;
	polyMoleculeCloud* molecules;
	clockTimer* evolveTimer;
	poly_solver_t() : redUnits(0), pot(0), molecules(0), evolveTimer(0) {}
	~poly_solver_t() {delete redUnits; delete pot; delete molecules; delete evolveTimer;}
};

