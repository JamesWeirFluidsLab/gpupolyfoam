/*---------------------------------------------------------------------------*\
 * =========                   |
 * \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
 *  \\    /    O peration      |
 *   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
 *    \\/      M anipulation   |
 * -------------------------------------------------------------------------------
 * License
 *    This file is part of OpenFOAM.
 * 
 *    OpenFOAM is free software; you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published by the
 *    Free Software Foundation; either version 2 of the License, or (at your
 *    option) any later version.
 * 
 *    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 *    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *    for more details.
 * 
 *    You should have received a copy of the GNU General Public License
 *    along with OpenFOAM; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * Application
 *    mdPolyFoam
 * 
 * Description
 *    Molecular dynamics solver for fluid dynamics -- polyatomics only
 * 
 * \*---------------------------------------------------------------------------*/

#include "gpuPolyFoam.h"
#include "gpuUtilities.h"
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
  
  solver->evolveTimer = new clockTimer(runTime,"evolve",true);
  
  
  //obtain reference properties
  #ifdef USE_OMM
  int num = 0;
  int nummols = 0;
  
  solver->plid = new polyIdPairs(mesh, *solver->pot);
  std::vector<int> molecularInfo;
  extractMolecularInfo(solver,molecularInfo);
  
  double dt = mesh.time().deltaT().value();
  const boundBox bBoxOF = mesh.bounds();
  setOMMBox(solver,bBoxOF,dt);
  initialiseOMM(solver);
  
  std::vector<OpenMM::Vec3> posInNm,velInNm,
    molPositions, moleculePI, sitePositions, momentOfInertia,siteForces;

  std::vector<OpenMM::Tensor> moleculeQ;
  std::vector<std::vector<unsigned int> > moleculeStatus;
  std::vector<OpenMM::Vec3> siteRefPositions;
  
  num = extractOFPostoOMM(posInNm,solver,bBoxOF);
  nummols = extractOFVeltoOMM(velInNm,solver,num);
  int t = 0;
  extractOFQ(solver,moleculeQ);
  t = extractOFSiteRefPositions(solver, siteRefPositions);
  t = extractMoleculePositions(solver, molPositions);
  t = extractMoleculePI(solver, moleculePI);
  extractMomentOfInertia(solver, momentOfInertia, moleculeStatus);

  Info << "extracted " << num 
  << " particles and " << nummols
  << " from OF" << nl;
  OpenMM::State state;    
  solver->context->setPositions(posInNm);
  solver->context->setMoleculeVelocities(velInNm);
  solver->context->setMoleculeQ(moleculeQ);
  solver->context->setMoleculePositions(molPositions);
  solver->context->setSiteRefPositions(siteRefPositions);
  solver->context->setMoleculePI(moleculePI);
  solver->context->setMomentOfInertia(momentOfInertia);
  solver->context->setMoleculeStatus(moleculeStatus);

  
  solver->ommTimer = new clockTimer(runTime, "openMMTimer", true);
  solver->openFoamTimer = new clockTimer(runTime, "openFoamTimer", true);
  solver->openFoamTimer->startClock();
  #endif

  Info << "\nStarting time loop\n" << endl;
  
  while (runTime.loop()){
    
    Info << "Time = " << runTime.timeName() << endl;
    
    solver->evolveTimer->startClock();
    
    #ifdef USE_OMM
    solver->molecules->preliminaries();
    solver->openFoamTimer->stopClock();
    solver->ommTimer->startClock();
    solver->integrator->step(1);
    velInNm.clear();
    state = solver->context->getState(State::MoleculeVel|
    State::MoleculePos|State::Forces|State::Positions,true);
    velInNm = state.getMoleculeVel();
    setOFVelocities(solver,velInNm);
    siteForces = state.getForces();
    setOFforce(solver,siteForces);
    molPositions = state.getMoleculePos();
    setOFPositions(solver,molPositions);
    sitePositions = state.getPositions();
    setOFSitePositions(solver,sitePositions);
    solver->ommTimer->startClock();
    solver->openFoamTimer->startClock();
    solver->molecules->postPreliminaries();

    #endif
    solver->evolveTimer->stopClock();
    
    runTime.write();
    
    Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
    
  }
  
  Info << "End\n" << endl;
  
  delete solver;
  
  return 0;
}


const Foam::word getScallingFunction(const MOLECULE& mol)
{
  Foam::word scallingFunction; //explicity namespace declaration 
  
  IOdictionary potentialDict
  (
    IOobject
    ("potentialDict",mol.mesh().time().system(),
     mol.mesh(),IOobject::MUST_READ,IOobject::NO_WRITE)
  );
  
  const dictionary& pairdict = potentialDict.subDict("pair");
  
  try
  {
    if(pairdict.found("electrostatic")){
      const dictionary& pairpotentialdict = pairdict.subDict("electrostatic");
      Foam::word temp = pairpotentialdict.lookup("energyScalingFunction");
      scallingFunction = temp;
  }
  else{
    throw 20;
  }
}
catch(int msg){
  Info << "getScallingFunction():: " << msg << nl;
}


return scallingFunction;
}
