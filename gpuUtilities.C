#include "gpuUtilities.h"


int initialiseOMM(struct poly_solver_t* solver)
{
	int status = 0;
	// Load OpenMM platform-specific plugins:
     	Platform::loadPluginsFromDirectory(Platform::getDefaultPluginsDirectory());
	int numPlatforms = Platform::getNumPlatforms();
	cout<<numPlatforms<<" OMM platforms found on this system"<<std::endl;
	for(int i=0;i<numPlatforms;i++)
   	{	
        	Platform& tempPlatform = OpenMM::Platform::getPlatform(i);
		std::string tempPlatformName = tempPlatform.getName();
		cout << "Platform " << i << 
		" is " << tempPlatformName.c_str() << std::endl;
	}
	cout<<"==================================\n"<<std::endl;
	solver->system = new System();
	MOLECULE& tempmol = *solver->molecules;

        //get the scalling function required for this simulation
        Foam::word scallingfunction = getScallingFunction(tempmol);
        
	//string variables for eps and sigma string
	std::string epsstring,sigmastring;
	extractLennardJonesParameters(tempmol,*solver->plid,epsstring,sigmastring);	
	std::cout<<"Eps string = "<<epsstring<<"\n"<<"Sigma string = "<<sigmastring<<"\n";
	
    CustomNonbondedForce* nonbonded;
    
    if(scallingfunction == "noScaling"){
        nonbonded = new CustomNonbondedForce(
        "4*eps*((sigma/r)^12-(sigma/r)^6)+(138.9354561469*q/r); "
        "eps= "+epsstring+"; sigma=" + sigmastring +"; q=q1*q2"
        );
    }
    else if(scallingfunction == "shifted"){
        nonbonded = new CustomNonbondedForce(
                "4*eps*((sigma/r)^12-(sigma/r)^6)+(138.9354561469*q)*(1/r-1/rCut); "
                "eps= "+epsstring+"; sigma=" + sigmastring +"; q=q1*q2"
                );
        nonbonded->addGlobalParameter("rCut",solver->rCutInNM);
    }
    else if(scallingfunction == "shiftedForce"){
        // ================================================
        // Coulomb potential with shifted force correction:
        // ================================================
        nonbonded = new CustomNonbondedForce(
                "4*eps*((sigma/r)^12-(sigma/r)^6)+(138.9354561469*q)*(1/r-1/rCut+(r-rCut)/rCut^2); "
                "eps= "+epsstring+"; sigma=" + sigmastring +"; q=q1*q2"
                );
        nonbonded->addGlobalParameter("rCut",solver->rCutInNM);
    }
    else{
        Info << "Error: (Solver) Cannot find electrostatic potential" <<
                "specification in potentialDict" << nl;
        exit(-1);
    }
    
    //add per particles to nonbonded force
    addParticlesToNonBonded(nonbonded,solver);
    //set OMM box size
    solver->system->setDefaultPeriodicBoxVectors(
        Vec3(solver->bBoxOMMinNm[0],0,0),
        Vec3(0,solver->bBoxOMMinNm[0],0),
        Vec3(0,0,solver->bBoxOMMinNm[0])
        );
        
	extractOFParticles(solver,nonbonded);
	solver->system->addForce(nonbonded);
    Platform& platform = Platform::getPlatformByName("OpenCL");
    solver->integrator =  new VerletIntegrator(solver->deltaT*OpenMM::PsPerFs);
    solver->context = new Context(*omm->system,*omm->integrator,platform);
    
    Info << "Initialised sytem on OMM with " <<
        solver->system.getNumParticles() <<
        " particles " << nl;
    Info << "Using OMM "<<
        omm->context->getPlatform().getName().c_str() <<
    << " platform " << nl;
    
	return status;
}


void extractLennardJonesParameters(const MOLECULE& mol,
			const polyIdPairs& p,
			std::string& epsString,
			std::string& sigmaString)
{
	epsString = "";
	sigmaString = "";

	const List<word>& idList = mol.pot().siteIdList();
	if(idList.size()<=0){
		FatalErrorIn("ExtractLennardJonesParameters::gpuPolyFoam") << nl
			<< "potential siteidlist cannot be zero "
			<< nl << abort(FatalError);
	}

	int counter = 0;
	forAll(p.sigma(), i)
    	{
      		forAll(p.sigma()[i], j)
      		{
           		const word& idAStr = idList[i];
           		const word& idBStr = idList[j];

			scalar epsAB = p.epsilon()[i][j]*6.02214129e23/1e3f;
			scalar sigmaAB = p.sigma()[i][j]/NANSEC;

			std::string epsABStr;
			std::stringstream outEpsAB;
           		outEpsAB << epsAB;
           		epsABStr = outEpsAB.str();

           		std::string sigmaABStr;
           		std::stringstream outSigmaAB;
           		outSigmaAB << sigmaAB;
          		sigmaABStr = outSigmaAB.str();

           		std::string epsStr2 = epsABStr+"*("+idAStr+"1*"+idBStr+"2)";
           		std::string sigmaStr2 = sigmaABStr+"*("+idAStr+"1*"+idBStr+"2)";

           		if(counter == 0)
           		{
              			epsString = epsString+""+epsStr2;
              			sigmaString = sigmaString+""+sigmaStr2;
           		}
           		else
           		{
             			 epsString = epsString+" + "+epsStr2;
              			 sigmaString = sigmaString+" + "+sigmaStr2;
           		}
           		counter++;
		}
	}	
		
}

void extractOFParticles(struct poly_solver_t* solver,
                        CustomNonbondedForce* const nonbonded)
{
    const int species = solver->plid->sigma().size();
    Info << "Number of Special found " << species << nl;
    
    MOLECULE& mol = *solver->molecules;
    IDLList<polyMolecule>::iterator m(mol.begin());
    for (m = mol.begin(); m != mol.end(); ++m)
    {
            const polyMolecule::constantProperties& constprop = 
                    mol.constProps(m().id());
            int molsize = constprop.sites().size();

            //now traverse throught the N sites of each molecule size obtained
            //from previous retrive
            int itr = 0;
            while(itr<molsize)
            {
                scalar tempmolmass = constprop.sites()[itr].siteMass()
                        *solver->refMass/DALTON;
                scalar tempmolcharge = constprop.sites()[itr].siteCharge()
                        *solver->refCharge/CHARGE;
                int midx = solver->system->addParticle(tempmolmass1);
                int sid = constprop.sites()[itr].siteId();
                
                std::vector<double> params(species+1);
                for(int k=0;k<species+1;k++)
                    params[k]=0;
                params[sid] = 1;
                // charge:
                params[species] = (double) tempmolcharge;
                nonbonded->addParticle(params);
                
                itr++;
            }
    }//first for loop ends
}


int extractOFPostoOMM(std::vector<Vec3>& posinnm,struct poly_solver_t* sol,
                      const boundBox& bb)
{
	posinnm.clear();
    int numMols = 0;
    int numParticles = 0;
    
	IDLList<polyMolecule>::iterator mol(sol->molecules->begin());
	for (mol = molecules->begin(); mol != molecules->end(); ++mol)
    {
        int molsize = mol().sitePositions().size();
        int m = 0;
        
        while (m<molsize) {
            Foam::vector rI = (mol().sitePositions()[m] - bb.min())*sol->refLength*NANSEC;
            posInNm.push_back(Vec3(rI.x(), rI.y(), rI.z()));
            m++;
        }
        numParticles += molsize;
        numMols++;
    }
    
    return numParticles;
 
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
	
	//convert to reduced omm rcut size
	solver->rCutInNM = solver->molecules->pot().pairPotentials().rCutMax()*solver->refLength/NANSEC;
	// Convert individual bounding-box length-scales to nanometres [nm]
	double Lx = bBoxOF.span().x()*solver->refLength/NANSEC;
	double Ly = bBoxOF.span().y()*solver->refLength/NANSEC;
	double Lz = bBoxOF.span().z()*solver->refLength/NANSEC;
	//create boxsize for OMM 
	solver->bBoxOMMinNm = Vec3(Lx,Ly,Lz);
	
}


void addParticlesToNonBonded(CustomNonbondedForce* const nonbonded,
                                const struct poly_solver_t* solver)
{
    const List<word>& idlist = solver->molecules->pot().siteIdList();
    
    forAll(solver->plid->sigma(),i){
        const word& idAstr = idlist[i];
        nonbonded->addPerParticleParameter(idAstr);
    }
    nonbonded->addPerParticleParameter("q");
    
    //set nonbondedmethod
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(solver->rCutInNM);
}
