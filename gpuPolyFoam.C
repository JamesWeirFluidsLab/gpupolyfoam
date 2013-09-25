/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
    \\/      M anipulation   |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    mdPolyFoam

Description
    Molecular dynamics solver for fluid dynamics -- polyatomics only

\*---------------------------------------------------------------------------*/

#include "gpuPolyFoam.h"

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


	poly_solver_t* solver = new poly_solver_t();

    solver->redUnits = new reducedUnits(runTime, mesh);

	solver->pot = new potential(mesh,*solver->redUnits);

	solver->molecules = new MOLECULE(runTime,mesh,*solver->pot,*solver->redUnits);
	
	solver->sid = new selectIdPairs(mesh, *solver->pot);
    
	solver->evolveTimer = new clockTimer(runTime,"evolve",true);
    
    double dt = mesh.time().deltaT().value();
    const boundBox bBoxOpenFOAM = mesh.bounds();

    
    // Bounding box in OpenMM in nanometres [nm]:
    Vec3  bBoxOpenMmInNm = Vec3(Lx, Ly, Lz);
    //obtain reference properties
#ifdef USE_OMM
    solver->refLenght = solver->redUnits.refLenght();
    solver->refMass = solver->redUnits.refMass();
    solver->refForce = solver->redUnits.refForce();
    solver->refCharge = solver->redUnits.refCharge();
    solver->refTime = solver->redUnits.refTime();
    // Convert to femtoseconds [fs] for OpenMM:
    solver->deltaT = = deltaT*redUnits.refTime()/FEM2SEC;
    // Convert individual bounding-box length-scales to nanometres [nm]:
    double Lx = bBoxOpenFOAM.span().x()*solver->refLength/NANSEC;
    double Ly = bBoxOpenFOAM.span().y()*solver->refLength/NANSEC;
    double Lz = bBoxOpenFOAM.span().z()*solver->refLength/NANSEC;
    //create boxsize for OMM
    solver->bBoxOMMinNm = Vec3(Lx,Ly,Lz);

#endif

    
    Info << "\nStarting time loop\n" << endl;

	while (runTime.loop()){
        
        Info << "Time = " << runTime.timeName() << endl;
		
        solver->evolveTimer->startClock();
		
        solver->molecules->evolve();
        
        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
        
	}
    
    Info << "End\n" << endl;
    
    delete solver;
    
    return 0;
}
