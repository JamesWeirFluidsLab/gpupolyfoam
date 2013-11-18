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
    Info << "Scalling function " << scallingfunction << nl;
       
	//dynamic coeff string array
	std::vector<std::string> coeffStr;
        
	extractLennardJonesParameters(tempmol,*solver->plid,coeffStr);
        std::string coefftype = solver->plid->coeffType();
        Info << "CoeffType " << coefftype << nl;
        
        const label cs = solver->plid->coeffSize();
        const List<label>& coeffids = solver->plid->coeffNumIds();
        const List<word>& coeffnames = solver->plid->coeffNames();
        
        
        std::stringstream tempstr;;
        
        for(int i=0;i<cs;++i){
            tempstr << coeffnames[i] << "= " << coeffStr[i] << ";";
        }
        
        std::string formulastr;
        
        if(coefftype == "lennardJones"){
            formulastr = 
            "4*epsilon*((sigma/r)^12-(sigma/r)^6)+(138.9354561469*q)*(1/r-1/rCut+(r-rCut)/rCut^2); " 
            + tempstr.str() + " q=q1*q2";
        }
        else if(coefftype == "morse"){
            formulastr = "D(exp[-2alpha(r-r0)]-2exp[-alpha(r-r0)]);"
                    + tempstr.str();
        }
        else{
            Info << "no coeff type found, quitting ..." << nl;
            exit(-1);
        }
                       
        //form a custom equation string from the obtained
        //parameters
        
	
	
    CustomNonbondedForce* nonbonded;
    
    /*if(scallingfunction == "noScaling"){
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
    else*/ if(scallingfunction == "shiftedForce"){
        // ================================================
        // Coulomb potential with shifted force correction:
        // ================================================
//        nonbonded = new CustomNonbondedForce(
//                "4*epsilon*((sigma/r)^12-(sigma/r)^6)+(138.9354561469*q)*(1/r-1/rCut+(r-rCut)/rCut^2); "
//                "epsilon= "+epsstring+"; sigma=" + sigmastring +"; q=q1*q2"
//                );
        nonbonded = new CustomNonbondedForce(formulastr);
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
    
    Info << "DEBUG: CUSTOMNONBONDEDFORCE INFO"<<nl;
    Info << "Num particles "<< nonbonded->getNumParticles() << nl;
    Info << "Num exclusion "<< nonbonded->getNumExclusions() << nl;
    Info << "Num perparticle params "<< nonbonded->getNumPerParticleParameters() << nl;
    Info << "Num global params "<< nonbonded->getNumGlobalParameters() << nl;
    Info << "Nonbonded method "<< nonbonded->getNonbondedMethod() << nl;
    Info << "Cutoff distance "<< nonbonded->getCutoffDistance() << nl;
    Info << "==============================="<<nl;    

    solver->system->addForce(nonbonded);
    Platform& platform = Platform::getPlatformByName("OpenCL");
    solver->integrator =  new VerletIntegrator(solver->deltaT*OpenMM::PsPerFs);
    solver->context = new Context(*solver->system,*solver->integrator,platform);
    
    Info << "Initialised sytem on OMM with " <<
        solver->system->getNumParticles() <<
        " particles " << nl;
    Info << "Using OMM "<<
        solver->context->getPlatform().getName().c_str() <<
    	" platform " << nl;
    
	return status;
}


void extractLennardJonesParameters(const MOLECULE& mol,
			const polyIdPairs& p,
			std::vector<std::string>& coeffStr)
{

	const List<word>& idList = mol.pot().siteIdList();
	const List<word>& coefflist = p.coeffNames();
	label coeffsize = coefflist.size();
	std::vector<std::string> coeffstr(coeffsize);

	int counter = 0;
	for(label c = 0; c < coeffsize; ++c)
	{
		counter = 0;
		forAll(p.coeffVals()[c],i){
			forAll(p.coeffVals()[c][i],j){
				const word& idAStr = idList[i];
				const word& idBStr = idList[j];
//TODO: delete later if all is well				
				scalar coeffAB;
                                
				if(coefflist[c] == "sigma")
                                    coeffAB = p.coeffVals()[c][i][j]/1e-9;
				else if(coefflist[c] == "epsilon")
                                    coeffAB = p.coeffVals()[c][i][j]*6.02214129e23/1e3f;
                                else if(coefflist[c] == "a")
                                    coeffAB = 0;//TODO: write conversion for a
                                else if(coefflist[c] == "alpha")
                                    coeffAB = 0;
                                else if(coefflist[c] == "r0")
                                    coeffAB = 0;
                                
				std::string ABstr;
				std::stringstream outAB;		
				outAB << coeffAB;
				ABstr = outAB.str();
				std::string tempstr = ABstr+"*("+idAStr+"1*"+idBStr+"2)";

				if(counter == 0)
					coeffstr[c] = coeffstr[c]+""+tempstr;
				else
					coeffstr[c] = coeffstr[c]+" + "+tempstr;			
				counter++;
			}
		}
	}

	for(counter = 0; counter < coeffsize; ++counter)
		coeffStr.push_back(coeffstr[counter]);
}

void extractOFParticles(struct poly_solver_t* solver,
                        CustomNonbondedForce* const nonbonded)
{
    const int species = solver->plid->nIds();
    
    int* midx=NULL;
	 
    MOLECULE& mol = *solver->molecules;
    IDLList<polyMolecule>::iterator m(mol.begin());
    for (m = mol.begin(); m != mol.end(); ++m)
    {
            const polyMolecule::constantProperties& constprop = 
                    mol.constProps(m().id());
            int molsize = constprop.sites().size();


            //allocate memory for array holding molecule index based on molecule size
            midx = (int*) malloc(sizeof(int)*molsize);

            //now traverse throught the N sites of each molecule size obtained
            //from previous retrive
            int itr = 0;
            while(itr<molsize)
            {
                double tempmolmass = constprop.sites()[itr].siteMass()
                        *solver->refMass/DALTON;
                double tempmolcharge = constprop.sites()[itr].siteCharge()
                        *solver->refCharge/CHARGE;
                midx[itr]  = solver->system->addParticle(tempmolmass);
                int sid = constprop.sites()[itr].siteId();
         
                std::vector<double> params(species+1);
                for(int k=0;k<species+1;k++)
                    params[k]=0;
                params[sid] = 1;
                params[species] = tempmolcharge;
                
                nonbonded->addParticle(params);
                
                itr++;
            }
            //add exclusion for each molecule obtained
            if(molsize>1)
			for(int i=0;i<molsize-1;i++)
				for(int j=i+1;j<molsize;j++)
					nonbonded->addExclusion(midx[i],midx[j]);
        	    // clear the mol id list for further assignment
            free(midx);

    }//first for loop ends
}


int extractOFPostoOMM(std::vector<Vec3>& posInNm,struct poly_solver_t* sol,
                      const boundBox& bb)
{
    posInNm.clear();
    int numMols = 0;
    int numParticles = 0;
    
    IDLList<polyMolecule>::iterator mol(sol->molecules->begin());
	for (mol = sol->molecules->begin(); mol != sol->molecules->end(); ++mol)
    {
        int molsize = mol().sitePositions().size();
        int m = 0;

        while (m<molsize) {
            Foam::vector rI = (mol().sitePositions()[m] - bb.min())*sol->refLength*NM;
            posInNm.push_back(Vec3(rI.x(), rI.y(), rI.z()));
            m++;
        }
        numParticles += molsize;
        numMols++;
    }
    
    return numParticles;
 
}

int setOFforce(struct poly_solver_t* solver, const std::vector<Vec3>& ommForces)
{
	IDLList<polyMolecule>::iterator mol(solver->molecules->begin());
	int numparticle = 0;

    for (mol = solver->molecules->begin(); mol != solver->molecules->end(); ++mol){
        int molsize = mol().siteForces().size();
        int m = 0;
		while(m<molsize){
			mol().siteForces()[m] = Foam::vector                      
                        (
                            ommForces[numparticle][0],
                            ommForces[numparticle][1],
                            ommForces[numparticle][2]
                        )*NM2RUF/solver->refForce;
			numparticle++;
			m++;
		}
	}
	return numparticle;
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
	solver->rCutInNM = solver->molecules->pot().pairPotentials().rCutMax()*solver->refLength/NANO;
	// Convert individual bounding-box length-scales to nanometres [nm]
	double Lx = bBoxOF.span().x()*solver->refLength/NANO;
	double Ly = bBoxOF.span().y()*solver->refLength/NANO;
	double Lz = bBoxOF.span().z()*solver->refLength/NANO;
	//create boxsize for OMM 
	solver->bBoxOMMinNm = Vec3(Lx,Ly,Lz);
	
}


void addParticlesToNonBonded(CustomNonbondedForce* const nonbonded,
                                const struct poly_solver_t* solver)
{
    const List<word>& idlist = solver->molecules->pot().siteIdList();
    
    for(int i = 0; i < solver->plid->nIds(); ++i){
        const word& idAstr = idlist[i];
	Info << "List " << i << ": "<< idAstr << nl;
        nonbonded->addPerParticleParameter(idAstr);
    }
    nonbonded->addPerParticleParameter("q");
    
    //set nonbondedmethod
    nonbonded->setNonbondedMethod(CustomNonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(solver->rCutInNM);
}

void getOMMState(const Context* context,std::vector<Vec3>& statearray)
//        enum STATES st)
{
	statearray.clear();
	State state = context->getState(OpenMM::State::Forces|State::Energy,true);
	statearray = state.getForces();
        
}
