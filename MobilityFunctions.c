/*

 \brief This is library to generate a mobility matrix based on Peskin kernel functions and staggered grid solver
 \date $April 15, 2015$
 \author Bakytzhan Kallemov, Courant Institute, NYU
 
This code accompanies the paper:
"An immersed boundary method for rigid bodies", B. Kallemov and A. Pal Singh Bhalla and B. E. Griffith and A. Donev, submitted to Comput. Methods Appl. Mech. Engrg., 2015
where the notation and ideas are explained.
These fits are made available to other users of the IB method -- let us know if they are useful or if you find problems.
We expect to release the complete code as part of the IBAMR library in the near future.
 
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

*/
#ifdef __cplusplus
extern "C"{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef KRON
#define KRON(i, j) ((i==j)?1:0) //Kronecker symbol
#endif

static const double MOB_FIT_FG_TOL = 1.0e-5;//min distance between blobs to apply empirical fitting 
static const double ZERO_TOL = 1.0e-10;//tolerance for zero value  
static double MOB_FIT_FACTOR; //constant for normalization
static int reUse=0; //flag for reuse data 
typedef enum _KERNEL_TYPES{IB3, IB4, IB6, UNKNOWN_TYPE = -1}KERNEL_TYPES;

//Hydro radii for each kernel IB3,IB4,IB6
const double MOB_FIT_HydroRadius[] = {0.91, 1.255, 1.4685};

static double HRad; //current hydrodynamic radius

//coefficients of empirical fit for steady stokes 
static double F_s[7];
static double G_s[4];
//coefficients of empirical fit for finite beta 
static double F_b[10];
static double G_b[6];
//coefficients of empirical fit for f(0)
static double Z_b[5];
#if (NDIM==2)
static double Z_s[3];
#endif


static KERNEL_TYPES GetKernelType(
    const char* IBKernelName)
{
    if (!strcmp(IBKernelName,"IB_3")) return IB3;
    if (!strcmp(IBKernelName,"IB_4")) return IB4;
    if (!strcmp(IBKernelName,"IB_6")) return IB6;
    
    fprintf (stderr, "MobilityFunctions: Unknown interpolation kernel type. Available options are IB_3, IB_4, IB_6. \n");
    exit (EXIT_FAILURE);
}

/*!
 * returns Hydrodynamic radius value
 */
double
getHydroRadius(
    const char* IBKernelName)
{
    KERNEL_TYPES CurrentKernelType=GetKernelType(IBKernelName);
    return MOB_FIT_HydroRadius[CurrentKernelType];
}


/*!
 * returns squared norm of the vector
 */
static double
get_sqnorm(
    const double *a_vec)
{
#if (NDIM==3)
    return a_vec[0]*a_vec[0]+a_vec[1]*a_vec[1]+a_vec[2]*a_vec[2];
#elif(NDIM==2)
  return a_vec[0]*a_vec[0]+a_vec[1]*a_vec[1];
#endif
}

/*!
 * returns a value based on linear interpolation of array
 */
static double InterpolateLinear(
    const double*Xin, 
    const double*Yin, 
    const int N, 
    double X0)
{
    if (fabs(X0-Xin[0])<ZERO_TOL) return Yin[0];
    if (X0<Xin[0]) return Yin[0]+(Yin[1]-Yin[0])/(Xin[1]-Xin[0])*(X0-Xin[0]);
    if (X0>Xin[N-1]) return Yin[N-1]+(Yin[N-1]-Yin[N-2])/(Xin[N-1]-Xin[N-2])*(X0-Xin[N-1]);
    //find nearest neighbour index
    int indx;
    for (indx=0;indx<N-1;indx++)
    {
	if (X0<=Xin[indx+1]) break;	
    }
    return Yin[indx]+(Yin[indx+1]-Yin[indx])/(Xin[indx+1]-Xin[indx])*(X0-Xin[indx]);
}

static void InterpolateConstants(
    KERNEL_TYPES MOB_FIT_current, 
    const double beta)
{
    //data for initialization

    //setting hydrodynamic radius
    HRad=MOB_FIT_HydroRadius[MOB_FIT_current];
    
#if (NDIM==3)
    //coefficients for steady stokes 3D
    const double F_stokes_coeff[][3] = {{1.769,  1.263,   1.097}, {1.637,  1.535,   1.60533}, {1.362,  0.8134,  0.566}, {0.6825, 0.1812,  0.0886}, {1.074,  0.6436,  0.4432}, {0.6814, 0.181,   0.08848}, {0.2609, 0.1774,  0.1567}};
    const double G_stokes_coeff[][3] = {{4.314,  11.08,  17.54}, {0.1106, 0.1759, 0.2365}};
    
    //coefficients for fitting f_beta(0) 3D
    //***********now exist only for IB6 kernel****************
    const int num_cases=8;
    //Number cases is 8, affects array sizes 
    double M_betas[] = {0., 0.1, 0.25, 0.5, 1., 10., 100., 1000.0};
    
    //coefficient to empirical fitting selfmobility f_beta(0)
    double F_beta_zero[]={1.0/0.0230502, -0.131445787110707,0.287509131179951,61.267737026060075};
    
    double F_beta_coeff[][8] = {{2.68030156113612,2.29464251115180,2.09297672079185,1.91583152825747,1.86066266441860,1.08621320400389,1.81206666469132,0.310021323508742}, 
					{1.39009017417989,1.40943006156447,1.35006023047952,1.94612277585268,3.86995647537320,11.5138612720591,1.11416207665514,-0.0372181749138520},
					{-0.884827389490450,-0.756475504578064,-0.608404711842046,-0.564784303490654,-0.651175884725500,-0.427005600208834,0.00756892008838000,0.104913106939363}, 
					{1.0,1.06262915172038,1.49440529353864,2.03210062812245,2.18657011795234,3.37689526048115,3.67794210129942,3.67794210129942}, 
					{2.92685539141722,2.67245966301464,2.42487239022736,2.70289015143191,3.34961045584841,3.07595994173387,0.174756918359867,0.129697932247353}, 
					{29.0233576593274,26.6542400026965,24.6035062102207,24.4821402826182,28.4998841491606,32.3675245644131,17.3333226235307,2.20851998706445}, 
					{1.67431017814765,1.55636776784166,1.39866915583103,1.25760476142372,1.26454914224257,1.54512892334494,1.04695565772442,0.470157376261480}, 
					{1.18077832670187,0.866186289480926,0.456603306920829,0.157941170393203,0.133374232562159,2.93244835654164,1.07958559850185,0.221651722533943}, 
					{0.0314600987678820,0.0296450899554450,0.0322750684446850,0.0296570052703910,0.0140708796915960,0.0426413220771180,0.0232409154706480,0.00157997335725900}, 
					{0.0320036255275870,0.0210726756307760,0.0110662759921420,0.00380530089032800,0.00265521991315500,0.00541862366430600,5.09744623703744e-05,0}};
    
    double G_beta_coeff[][8] = {{2.56932494798493,0.884354880187973,0.113504063191842,-0.221849887772746,0.542698003099015,0.276258520317914,0.330745197269599,-7.88722383366528e-05}, 
					{-0.171344575435234,0.431984145806452,0.482467326122353,0.461388580082004,-0.0463322271962110,-0.131688108575531,-0.114177160773643,0.0119064417188160}, 
					{-0.0846371631980100,-0.0746824104260310,-0.0822639676502030,-0.0689510022688410,0.108601073154871,0.105609456589694,0.0941995195204770,0.0575579089490820}, 
					{1.0,0.356391856987450,0.495958169969290,0.651495025037754,0.608684335056712,9.63355010424247,2.60454636502040,3.17882803024967}, 
					{2.21329860848447,0.894411380099355,0.698529765627973,0.552224660456298,0.723993778966245,0.423143094120794,0.359827255706337,0.267507608759639}, 
					{0., 0.00283198188374900,0.00380833037508200,0.00266380633550800,-0.00766798833782700,-0.00440457501794600,-0.00100982221976800,1.90168352243139e-05}};
    

#elif(NDIM==2)
//coefficients for fitting f_beta(0) 2D
    const int num_cases=7;
    const double F_beta_zero[3][5] ={{0.1245, 0.05237, 1.873, 0.5579, -0.0159}, {0.07016, 0.01957, 1.229, 0.2613, -0.00793}, {0.05307, 0.01246, 1.015, 0.1882,-0.005932}};
    const double F_stokes_zero[3][3]={{0., 0., 0.}, {0., 0., 0.}, {0.1383, 0.07947, 0.7942}};
    
    double M_betas[7] = {0., 0.1, 0.25, 0.5, 1.0, 5.0, 10};
    
    //coefficients for steady stokes in 2d
    const double F_stokes_coeff[3][5] =  {{0., 0., 0., 0., 0.}, {0., 0., 0., 0., 0.}, {0.011515464372641,0.009893387219365,0.049266457632092,-1.919065487969394,0.340774530921814}};
    const double G_stokes_coeff[3][4] = {{0., 0., 0., 0.}, {0., 0., 0., 0.}, {0.1658, 0.007702, -0.1309, 0.1752}};

    //coefficients for finite beta in 2d
    double F_beta_coeff[8][7] ={
	//IB6
	{0.401420183520072, 0.401420183520072, 0.556142821029073, 0.664888194861374, 0.685894454477942, 0.718841089908404, 0.779082390669803},
	{14.804050962552742, 14.804050962552742, 10.904841649143368, 18.767257210239800, 84.161092829863650, 88.264694040072570, 77.033837028058230},
	{-0.279298012100466, 0.138980245204639, -0.009865807538010, -0.105903623979614, -0.136643790173548, -0.108325904537003, -0.070939669895227},
	{0.143313203071711, -0.831842069290423, -0.502458949732120, -0.269698049364202, -0.127451784740875, 0.028550283976351,0.014242399241400},
	{-8.452399536533306e-04, 0.308011998536034, 0.198315636699142, 0.122821649591300, 0.069589352321454, -1.213280410389923e-05, -1.547596093318852e-05},
	{-0.346117229077536, 1.824423513673857, 1.290756405391252, 0.893777292699120, 0.601816678029090, 0.075988811036186, 0.120200851694187},
	{0.149778804192806, -1.104621035097404, -0.701248785909740, -0.410203830606991, -0.214828358133586, 0.017288913634831,0.004378615023975},
	{-9.654473927525071e-04, 0.318896526093837, 0.205898404796718, 0.127863548863031, 0.072415937247846, 1.721308517520611e-04, 1.217042463586860e-04}};

    double F0_coeff[8] ={0.401420183520072,14.804050962552742,-0.279298012100466,0.143313203071711,-8.452399536533306e-04,-0.346117229077536,0.149778804192806,-9.654473927525071e-04};

    
    double G_beta_coeff[6][7] = { 
	//IB6
	{41.140687788555795, 41.140687788555795, 41.140769660186210, 41.139661947715510, 41.139222653458140, 41.105326347582200, 47.160007792647610},
	{14.217072827929760, 14.217072827929760, 14.217087881280210, 14.216779422829347, 14.216641036782928, 14.204778284605013, 16.468575623609084},
	{423.3, 3.485776377241552e+02,2.012387396385607e+02, 1.229443140730820e+02, 1.378753563734097e+02, 2.730874202363090e+02, 3.571369756322808e+02},
	{-620.4, -4.416320420444578e+02, -1.892987328094158e+02, -59.978486591337024, -55.000065547079785, -93.500092223773600, -89.959564014124340},
	{353.1, 2.427725412415102e+02, 1.140535742418950e+02, 42.114254817059020, 27.352127041224016, 30.082757196477147, 31.289709513774948},
	{1.918, 1.669052214523884, 1.336813769864806, 0.959267692774396, 0.666701713724319, 0.291457599992040, 0.196141486054449}};
#endif

    // The current version supports only IB_6 kernel for the time dependent case
    if  ((MOB_FIT_current!=IB6) && beta)
    {    
	fprintf (stderr,"MobilityFunctions: WARNING! Current version does not support IB_3 and IB_4 for time dependent case.\n");
	//exit (EXIT_FAILURE);
    }

    int cnt;
#if (NDIM==3)
    //********3D case
    //setting coeficients for steady stokes fitting
    for (cnt=0;cnt<7;cnt++) F_s[cnt] = F_stokes_coeff[cnt][MOB_FIT_current];
    for (cnt=0;cnt<2;cnt++) G_s[cnt] = G_stokes_coeff[cnt][MOB_FIT_current];
    
    //setting coeficients for time dependent fitting (curently only for IB6 kernel)
    for (cnt=0;cnt<4;cnt++)  Z_b[cnt] = F_beta_zero[cnt];
    for (cnt=0;cnt<10;cnt++) F_b[cnt] = InterpolateLinear(M_betas, F_beta_coeff[cnt], num_cases, beta);
    for (cnt=0;cnt<6;cnt++)  G_b[cnt] = InterpolateLinear(M_betas, G_beta_coeff[cnt], num_cases, beta);
   
#elif(NDIM==2)
    for (cnt=0;cnt<5;cnt++) Z_b[cnt] = F_beta_zero[MOB_FIT_current][cnt];
    for (cnt=0;cnt<3;cnt++) Z_s[cnt] = F_stokes_zero[MOB_FIT_current][cnt];
    for (cnt=0;cnt<5;cnt++) F_s[cnt] = F_stokes_coeff[MOB_FIT_current][cnt];
    for (cnt=0;cnt<4;cnt++) G_s[cnt] = G_stokes_coeff[MOB_FIT_current][cnt];

    if (beta>=ZERO_TOL)
    	for (cnt=0;cnt<8;cnt++) F_b[cnt] = InterpolateLinear(M_betas, F_beta_coeff[cnt], num_cases, beta);
    else
	for (cnt=0;cnt<8;cnt++) F_b[cnt] = F0_coeff[cnt];
    
    for (cnt=0;cnt<6;cnt++) G_b[cnt] = InterpolateLinear(M_betas, G_beta_coeff[cnt], num_cases, beta);
    printf("beta is %g\n",beta);
#endif
    return;
};

/*!
 * initialization of all constants. 
 */
static void InitializeAllConstants(
    const char* IBKernelName,
    const double MU,
    const double rho,
    const double Dt,
    const double DX)
{

    KERNEL_TYPES CurrentKernelType=GetKernelType(IBKernelName);
    double beta;
    //finding beta
    if (MU <= ZERO_TOL) 
	beta=0.0; //invisid case
    else 
	beta=MU*Dt/(rho*DX*DX);
    
#if (NDIM==3)
//******3D case
    if ((rho < ZERO_TOL)||(beta>=1000.1)) 
	MOB_FIT_FACTOR = 1./MU/DX;  //3D steady stokes
    else  
	MOB_FIT_FACTOR = Dt/(rho*DX*DX*DX);
#elif(NDIM==2)
//*******2D case
    if ((rho<ZERO_TOL)||(beta>10.1)) //2D steady stokes
	MOB_FIT_FACTOR = 1./MU;  
    else
	MOB_FIT_FACTOR = Dt/(rho*DX*DX);
#endif
    
    InterpolateConstants(CurrentKernelType, beta);
}    

/* !
**  _F_R_INF
*/
static double 
_F_R_INF(
    const double rr,
    const double Dx,
    const double L_domain)
{
    const double r=rr/Dx;
#if (NDIM==3)
    const double factor=1.0/(8.0*M_PI);
    if (r<0.8) 
	return factor/(3.0/4.0*HRad+F_s[6]*r*r);
    else
	return factor*(exp(-F_s[0]*r)*F_s[1]/HRad +(F_s[2]*r+F_s[3]*r*r*r)/(1.0 + F_s[4]*r*r + F_s[5]*pow(r,4)));
#elif(NDIM==2)
    if (L_domain<ZERO_TOL) 
    {
	fprintf(stderr, "IBEMpiricalMobility:_F_R_INF()  L_domain must be non specified in 2D!. Abort.\n");
	exit(EXIT_FAILURE); 
    }
    double f_0=(Z_s[0]+Z_s[1]*log(L_domain/Z_s[2]));
    if (r<0.1)
	return f_0;
    else
	return (f_0+(F_s[0]*r*r+F_s[1]*r*r*r+F_s[2]*r*r*r*log(r))/(1.0+F_s[3]*r+F_s[4]*r*r-F_s[2]*4*M_PI*r*r*r));
#endif
}// _F_R_INF

/* !
**  _G_R_INF
*/
static double 
_G_R_INF(
    const double rr,
    const double Dx)
{
    const double r=rr/Dx;
#if (NDIM==3)
    const double factor=1.0/(8.0*M_PI);
    if (r <  MOB_FIT_FG_TOL) return 0.0;
    else  return (factor*r*r)/(G_s[0] + G_s[1]*r*r + r*r*r);
#elif(NDIM==2)
    return  (G_s[0]*r*r+G_s[1]*r*r*r)/(1.0+G_s[2]*r+G_s[3]*r*r+G_s[1]*r*r*r)/4/M_PI;
#endif
}// _G_R_INF


/* !
**  _F_R_BETA
*/
static double 
_F_R_BETA( 
    const double rr,
    const double Dx,
    const double beta,
    const double L_domain)
{
    const double r=rr/Dx;

#if (NDIM==3)
    double f_0 = (1.0+Z_b[1]*sqrt(beta)+Z_b[2]*beta)/(Z_b[0] + Z_b[3]*beta+Z_b[2]*6.0*M_PI*HRad*beta*beta);

    //invisid case
    if (beta < ZERO_TOL) 
	return f_0/4.0/M_PI*((4.0*M_PI-F_b[4]*r*r)/(1.0+F_b[0]*r+F_b[1]*r*r+F_b[2]*r*r*r+F_b[4]*f_0*pow(r,5)) + (F_b[5]*r*exp(-F_b[6]*r)+F_b[7]*r)/(1.0+F_b[8]*r*r*r+F_b[9]*pow(r,5)));
    else if (beta<1000.1)
	return f_0/4.0/M_PI*((4.0*M_PI+F_b[4]*(-r*r+pow(r,4)*exp(-F_b[3]*r/sqrt(beta))/(2.0*beta)))/(1.0+F_b[0]*r+F_b[1]*r*r+F_b[2]*r*r*r+F_b[4]*f_0*pow(r,5)) + (F_b[5]*r*exp(-F_b[6]*r)+F_b[7]*r)/(1.0+F_b[8]*r*r*r+F_b[9]*pow(r,5)));
    else
	return 	_F_R_INF(rr,Dx,0.0);

#elif(NDIM==2)
    const double factor = -0.5/M_PI;
    if (beta<=10.1)
    {
	if (r< 0.8) //if distance is less self mobility is used 
	    return (Z_b[0]+Z_b[1]*beta*beta)/(Z_b[4]*pow(beta,4) + Z_b[3]*beta*beta*beta + Z_b[2]*beta + 1.0)/(1.0+r*r/(1.7+1.5*log(1.0+beta)));
	else 
	    return factor*( exp(-F_b[0]*r/sqrt(beta?beta:1.0))*r*log(r)/(F_b[1]+2.0*r)/(beta?beta:1.0) + (F_b[2]+F_b[3]*r+F_b[4]*r*r)/(1.0+F_b[5]*r*r+F_b[6]*r*r*r+F_b[7]*pow(r,4)));
    }
    else
    {
	return 	_F_R_INF(rr,Dx,L_domain);
    }

#endif
}

/* !
**  _G_R_BETA
*/
static double 
_G_R_BETA( 
    const double rr,
    const double Dx,
    const double beta)
{
    const double r=rr/Dx;

#if (NDIM==3)
    double f_0 = (1.0+Z_b[1]*sqrt(beta)+Z_b[2]*beta)/(Z_b[0] + Z_b[3]*beta+Z_b[2]*6.0*M_PI*HRad*beta*beta);
    if (beta < ZERO_TOL) 
	return f_0*3.0/(4.0*M_PI)*G_b[4]*r*r/(1.0+G_b[0]*r+G_b[1]*r*r+G_b[2]*r*r*r+G_b[4]*f_0*pow(r,5));
    else if (beta<1001.0)
	return f_0*3.0/(4.0*M_PI)*G_b[4]*(r*r+pow(r,4)*exp(-G_b[3]*r/sqrt(beta))/(6.0*beta))/(1.0+G_b[0]*r+G_b[1]*r*r+G_b[2]*r*r*r+G_b[5]*pow(r,4)+G_b[4]*f_0*pow(r,5));
    else 
	return 	_G_R_INF(rr,Dx);
#elif(NDIM==2)
    if (beta<=10.1)
    {
	const double factor =1.0/M_PI;
	if (r< MOB_FIT_FG_TOL) 	return 0.0;
	else 
	{
	    double g_value = factor*r/(exp(-G_b[5]*r)*(G_b[2]+G_b[3]*r+G_b[4]*r*r)+r*r*r);
	    if (beta > ZERO_TOL) g_value +=factor/beta*exp(-G_b[0]*r/sqrt(beta))*r/(G_b[1]+4.0*r);
	    return g_value;
	}
    }
    else
    {
	return _G_R_INF(rr,Dx);
    }
#endif
}

/*!
 * computes Empirical Mobility components f(r) and g(r)
 */
void getEmpiricalMobilityComponents(
    const char* IBKernelName,
    const double MU,
    const double rho,
    const double Dt,
    const double r,
    const double DX,
    const int resetAllConstants,
    const double L_domain,
    double *F_MobilityValue,
    double *G_Mobilityvalue)
{
    //Reuse same static constants for efficiency
    if (resetAllConstants) reUse=0;
    if (!reUse)
    {
	InitializeAllConstants(IBKernelName, MU, rho, Dt,DX);
	reUse=1;
    }
    double beta;
    //finding beta
    if (MU <= ZERO_TOL) 
	beta=0.0; //invisid case
    else 
	beta=MU*Dt/(rho*DX*DX);

    if (rho < ZERO_TOL)
    {
	*F_MobilityValue = MOB_FIT_FACTOR*_F_R_INF(r,DX, L_domain); //steady stokes term for f(r)
	*G_Mobilityvalue = MOB_FIT_FACTOR*_G_R_INF(r,DX);//steady stokes term for g(r)
    }
    else
    {
	*F_MobilityValue = MOB_FIT_FACTOR*_F_R_BETA(r,DX,beta, L_domain); //time-dependent f(r)
	*G_Mobilityvalue = MOB_FIT_FACTOR*_G_R_BETA(r,DX,beta);//time-dependent g(r)
    }
    return;
}

/*!
 * returns Self-Mobility value
 */
double getEmpiricalSelfMobility(
    const char* IBKernelName,
    const double MU,
    const double rho,
    const double Dt,
    const double DX,
    const int resetAllConstants,
    const double L_domain)
{
    if (resetAllConstants) reUse=0;
    if (!reUse)
    {
	InitializeAllConstants(IBKernelName, MU, rho, Dt, DX);
	reUse=1;
    }
    double beta;
    //finding beta
    if (MU <= ZERO_TOL) 
	beta=0.0; //invisid case
    else 
	beta=MU*Dt/(rho*DX*DX);
    
    if (rho < ZERO_TOL)
	return MOB_FIT_FACTOR*_F_R_INF(0.,DX, L_domain);
    else
	return MOB_FIT_FACTOR*_F_R_BETA(0., DX,beta, L_domain); 
}


/*!
 * generates and returns Empirical Mobility Matrix
 */
void
getEmpiricalMobilityMatrix(
    const char* IBKernelName,
    const double MU,
    const double rho,
    const double Dt,
    const double DX,
    const double *X, 
    const int N,
    const int resetAllConstants,
    const double PERIODIC_CORRECTION,
    const double L_domain,
    double *MM)
{

    int row,col;
    for (row = 0; row < N; row++)
	for (col = 0; col <=row; col++)    
	{
	    double r_vec[NDIM];
	    const int size = N*NDIM;
	    int cdir;
	    for (cdir = 0; cdir < NDIM; cdir++) 
	    {
		r_vec[cdir] = X[row*NDIM+cdir] - X[col*NDIM+cdir]; //r(i) - r(j)
	    }
	    
	    const double rsq = get_sqnorm(r_vec); 
	    const double r   = sqrt(rsq);
	    double F_R, G_R;
	    
	    getEmpiricalMobilityComponents(IBKernelName, MU, rho, Dt, r, DX, resetAllConstants, L_domain, &F_R, &G_R);

	    int idir, jdir;
	    for (idir = 0; idir < NDIM; idir++)
		for (jdir = 0; jdir <=idir; jdir++) 
		{
		    const int index = (col*NDIM+jdir)*size +row*NDIM+idir; // column-major for LAPACK
		    MM[index] = F_R*KRON(idir,jdir);
		    if (row != col) 
		    {
			MM[index] += G_R*r_vec[idir]*r_vec[jdir]/rsq;
			MM[(row*NDIM+idir)*size +col*NDIM+jdir]=MM[index];
			if  (idir!=jdir)
			    MM[(row*NDIM+jdir)*size +col*NDIM+idir]=MM[index];
		    }
		    if  (idir!=jdir)
			MM[(col*NDIM+idir)*size +row*NDIM+jdir]=MM[index];
		}
	}
    return;
}


/*!
 * generates and returns RPY Mobility Matrix
 */
void
getRPYMobilityMatrix(
    const char* IBKernelName,
    const double MU,
    const double DX,
    const double *X, 
    const int N,
    const double PERIODIC_CORRECTION,
    double *MM)
{
    HRad=getHydroRadius(IBKernelName)*DX;
    const double mu_tt = 1./(6.0*M_PI*MU*HRad);
    
    double r_vec[NDIM];
    int size = N*NDIM;
    int row,col;
    for (row = 0; row < N; row++)
	for (col = 0; col <=row; col++)    
	{
	    if (row==col)
	    {
		int idir,jdir;
		for (idir = 0; idir < NDIM; idir++)
		    for (jdir = 0; jdir < NDIM; jdir++) 
		    {
			const int index = (col*NDIM+jdir)*size +row*NDIM+idir; // column-major for LAPACK
			MM[index] = (mu_tt- PERIODIC_CORRECTION)*KRON(idir,jdir);
		    }
	    }
	    else
	    {
		int cdir;
		for (cdir = 0; cdir < NDIM; cdir++) 
		{
		    r_vec[cdir] = X[row*NDIM+cdir] - X[col*NDIM+cdir]; //r(i) - r(j)
		}
		
		const double rsq = get_sqnorm(r_vec); 
		const double r   = sqrt(rsq);
		int idir,jdir;
		for (idir = 0; idir < NDIM; idir++)
		    for (jdir = 0; jdir <=idir; jdir++) 
		    {
			const int index = (col*NDIM+jdir)*size +row*NDIM+idir; // column-major for LAPACK
			if (r<=2.0*HRad)
			{
			    MM[index] = (mu_tt*(1-9.0/32.0*r/HRad)- PERIODIC_CORRECTION)*KRON(idir,jdir)
				+mu_tt*r_vec[idir]*r_vec[jdir]/rsq*3.0*r/32./HRad;
			    MM[(row*NDIM+idir)*size +col*NDIM+jdir] = MM[index];
			    if (idir!=jdir)
			    {
				MM[(col*NDIM+idir)*size +row*NDIM+jdir] = MM[index];
				MM[(row*NDIM+jdir)*size +col*NDIM+idir] = MM[index];
			    }
			}
			else
			{
			    double cube=HRad*HRad*HRad/r/r/r;
			    MM[index] = (mu_tt*(3.0/4.0*HRad/r+1.0/2.0*cube) - PERIODIC_CORRECTION)*KRON(idir,jdir)
				+ mu_tt*r_vec[idir]*r_vec[jdir]/rsq*(3.0/4.0*HRad/r-3.0/2.0*cube);
			    MM[(row*NDIM+idir)*size +col*NDIM+jdir] = MM[index];
			    if (idir!=jdir)
			    {
				MM[(col*NDIM+idir)*size +row*NDIM+jdir] = MM[index];
				MM[(row*NDIM+jdir)*size +col*NDIM+idir] = MM[index];
			    }

			}
		    }//jdir
	    }
	}//column loop
    return;
}

#ifdef __cplusplus
}
#endif //extern "C"
