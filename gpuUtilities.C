#include "gpuUtilities.h"

int extractOFPostoOMM(std::vector<Vec3>& posinnm,struct poly_solver_t* sol)
{
	posinnm.clear();
	IDLList<polyMolecule>::iterator mol(sol->molecules->begin());
	for (mol = molecules->begin(); mol != molecules->end(); ++mol)
        {
            numOfAtoms++;

            Foam::vector rI = (mol().position() - globalBBox.min())*sol->refLength/NANSEC;

            posInNm.push_back(Vec3(rI.x(), rI.y(), rI.z()));
        } 
}


//set OMM box
void setOMMBox(struct poly_solver_t* solver, const boundBox& bBoxOF,const double dt)
{
	solver->refLength = solver->redUnits->refLength();
	solver->refMass = solver->redUnits->refMass();
	solver->refForce = solver->redUnits->refForce();
	solver->refCharge = solver->redUnits->refCharge();
	solver->refTime = solver->redUnits->refTime();
	// Convert to femtoseconds [fs] for OpenMM: 
	solver->deltaT = dt*solver->refTime/FEM2SEC;
	// Convert individual bounding-box length-scales to nanometres [nm]
	double Lx = bBoxOF.span().x()*solver->refLength/NANSEC;
	double Ly = bBoxOF.span().y()*solver->refLength/NANSEC;
	double Lz = bBoxOF.span().z()*solver->refLength/NANSEC;
	//create boxsize for OMM 
	solver->bBoxOMMinNm = Vec3(Lx,Ly,Lz);
	
}


