# RigidBodyIB
Rigid Body Immersed Boundary Method

This code provides implementations of the empirical fits to the marker-marker (blob-blob) mobility described in the Appendices in the paper:

**An immersed boundary method for rigid bodies**, B. Kallemov and A. Pal Singh Bhalla and B. E. Griffith and A. Donev, Communications in Applied Mathematics and Computational Science, **11(1)**, 79-141 (2016)
[DOI](http://dx.doi.org/10.2140/camcos.2016.11.79) [arXiv](http://arxiv.org/abs/1602.02170)

These fits are made available to other users of the IB method -- let us know if they are useful or if you find problems.
This code is now part of the [IBAMR library](https://github.com/ibamr/ibamr).
 
For steady Stokes flow (beta=infinity) fitting coefficients are available for the
3 (Roma+Peskin), 4 (Peskin) and 6 (Bao+Peskin) point kernels
WARNING: For unsteady Stokes, i.e., finite viscous CFL number beta, fitting is available *only* for IB6 kernel (Bao+Peskin)
 
Note that instead of taking beta as an argument these functions take viscosity, grid spacing and time step.
This is because there are some cancelations of zeros and infinities for beta=0 or infinity that are best avoided numerically.
Also note that the viscosity argument to these functions is the physical viscosity multiplied by a constant kappa,
where kappa depends on the temporal discretization used
(e.g., kappa=1 for Backward Euler, kappa=1/2 for Crank-Nicolson, and kappa may change from stage to stage for RK integrators)

Dimension: NDIM identifier must be defined before any call of the functions
 
The following function are available 

getHydroRadius
  - Returns a hydrodynamic radius of a single "blob" for the corresponding Peskin kernel (in units of grid spacing h)
  input parameters 
  IBKernelName  name of the kernel {"IB3","IB4","IB6"}  

getEmpiricalMobilityComponents
  - Generates empirical f(r) and g(r) functions values for a the NDIM X NDIM block of Mobility Matrix for two markers 
  input parameters
  IBKernelName  		name of the kernel {"IB3","IB4","IB6"}  
  MU			fluid viscosity  
  rho			fluid density 
  Dt			time step used in the fluid solver
  r			distance between markers
  DX			grid spacing
  resetAllConstants 	whether all constant must be reset if beta(viscous CFL number) is changed (otherwise will use the previous beta for fitting formula)
  L_domain		length of the domain (used only only for 2 dimensional steady stokes)
  F_MobilityValue		a pointer to return f(r) value
  G_Mobilityvalue		a pointer to return g(r) value


getEmpiricalMobilityMatrix
  - Generates an empirical Mobility Matrix
  input parameters
  IBKernelName  		name of the kernel {"IB3","IB4","IB6"}  
  MU			fluid viscosity  
  rho			fluid density 
  Dt			time step used in the fluid solver
  DX			grid spacing
  X			a pointer to array of size NDIM*N that containts coordinates of markers 
  N			numbers of markers
  resetAllConstants 	whether all constant must be reset if beta(viscous CFL number) is changed (otherwise will use the previous beta for fitting formula)
  PERIODIC_CORRECTION	periodic corrections for f(r) function, set to 0.0 for all other cases or if it's not known.
  L_domain		length of the domain (used only only for 2 dimensional steady stokes)
  MM			a pointer to return a moblity matrix stored in column-major format. The allocated space size for MM must be non less than sizeof(double)*(NDIM*N)^2 

    
    
getRPYMobilityMatrix
  - Generates an Mobility Matrix based on Ronte-Prager-Yamakawa (RPY) approximation
  input parameters
  IBKernelName  		name of the kernel {"IB3","IB4","IB6"}  
  MU			fluid viscosity  
  DX			grid spacing
  X			a pointer to array of size NDIM*N that containts coordinates of markers 
  N			numbers of markers
  PERIODIC_CORRECTION	periodic corrections for f(r) function, set to 0.0 for all other cases or if it's not known.
  MM			a pointer to return a moblity matrix stored in column-major format. The allocated space size for MM must be non less than sizeof(double)*(NDIM*N)^2 

