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
#ifdef USE_OMM
#include "gpuUtilities.h"
#include "gpuUtilities.C"
#endif

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


        poly_solver_t* solver = new poly_solver_t();

        solver->redUnits = new reducedUnits(runTime, mesh);

	solver->pot = new potential(mesh,*solver->redUnits);

	solver->molecules = new MOLECULE(runTime,mesh,*solver->pot,*solver->redUnits);

//	solver->sid = new selectIdPairs(mesh, *solver->pot);
    
        double dt = mesh.time().deltaT().value();
        const boundBox bBoxOF = mesh.bounds();

    
    //obtain reference properties
#ifdef USE_OMM
        setOMMBox(solver,bBoxOF,dt);
	//extract OF positions to OMM
	std::vector<Vec3> ommpos;
	extractOFPostoOMM(ommpos,solver);
#endif

    
    Info << "\nStarting time loop\n" << endl;
	solver->evolveTimer = new clockTimer(runTime,"evolve",true);
      while (runTime.loop()){
        
        Info << "Time = " << runTime.timeName() << endl;
		
        solver->evolveTimer->startClock();
#ifdef USE_OMM
        printf("OMM");
#else
        solver->molecules->evolve();
#endif
        
        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
        
	}
  
    Info << "End\n" << endl;
    
    delete solver;
    
    return 0;
}
