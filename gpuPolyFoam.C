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


	poly_solver_t* poly = new poly_solver_t();

    	poly->redUnits = new reducedUnits(runTime, mesh);

	poly->pot = new potential(mesh,*poly->redUnits);

	poly->molecules = new polyMoleculeCloud(runTime,mesh,*poly->pot,*poly->redUnits);
	poly->evolveTimer = new clockTimer(runTime,"evolve",true);

    Info << "\nStarting time loop\n" << endl;

	while (runTime.loop()){
        	Info << "Time = " << runTime.timeName() << endl;
		poly->evolveTimer->startClock();
		poly->molecules->evolve();
        	runTime.write();

        	Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	}
    Info << "End\n" << endl;
/*
    potential pot(mesh, redUnits);

    polyMoleculeCloud molecules(runTime, mesh, pot, redUnits);

    Info << "\nStarting time loop\n" << endl;

    clockTimer evolveTimer(runTime, "evolve", true);

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << endl;

        evolveTimer.startClock();

        molecules.evolve();

        evolveTimer.stopClock();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info << "End\n" << endl;
*/
delete poly;
    return 0;
}
