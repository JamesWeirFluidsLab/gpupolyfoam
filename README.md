**GPU PolyFOAM**
================

This OF solver is coupled with OpenMM GPU library to delegate functionality on
GPU. This solver is intended to run polymolecule simulations consisting molecule
sizes of N. It's a generic solver which enables running MonoAtomic and
polyAtomic simulations.



**Functionality available**
---------------------------

\* Force Calculation on GPU



**Download**
------------

Please follow steps to download Gpu PolyFoam solver.

-   *change directory to *`cd
    $FOAM_APP/solvers/discreteMethods/molecularDynamics/  `

-   on changing directory run this command `git clone
    git@bitbucket.org:saifmulla/gpupolyfoam.git`

-   if the comand executes successfully you must see a new directory created
    automantically with the name **gpupolyfoam**



**Installation** 
-----------------



This installation assumes you already have OpenMM installed in your system and
its path

is set in your bashrc. Moreover you must also create another environment
variable inside your bashrc file which essentially specifies the parent
directory of OpenMM installation.

inside you bashrc create a varible `export FOAM_OMM=<path-to-omm-install-dir>`



*Note: Do not forget to source your bashrc to reflect this new variable*



Essentially the *FOAM_OMM* varible is used as a relative path inside the

gpuPolyFoam options file, this is done to keep the installation standard and
error free therefore please do not change this variable inside the options file



Please execute the following step for installing **GPU PolyFOAM** solver

-   change directory to the newly created **gpupolyfoam** directory

-   inside the directory execute *wmake USE_OMM=POLY* for polyAtomic class



or alternatively for motoatomic you must execute

    *wmake USE_OMM=MONO*

*Note: please take care of upper case letters which are vital for appropriate
compilation*



if you do not get any compilation errors then there should be an executable
generated inside

*$FOAM_APPBIN* folder with the name **gpuPolyFoam**

you should use the same executable name to execute GPU based polyAtomistic
simulation.



                    **Thanks**

==================


