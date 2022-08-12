/*
**********************************************************
Enviroment settings for IMPACTA code
Including boundary conditions and differential info

1
AGRT

22/2/08 - Ion motion added

12/2/07

12/3/07 - The boundary condition stuff in here doesn't 
really do anything. The code is in the operators .h file
The other namespaces do though.

200 - MPI messages communicating ParVecs.

17/4/07 - more changes - diagnostic variables
3/5/07  - vector checking code inserted at end

27/7/07 Many small additions
15/10/07 - MPIChunkSize added as constant
feb 08 - iterations for matrix added
1/7/08 - Fixed boundaries changed to simpler algorithm
**********************************************************

*/
/*
Later need to read some of these namespaces from a file.

These are namespaces for different differential types
The name is fairly self explanatory
the codes are
diff_type
20 - centred - (as  in Roberts paper 12/2/07)
10 - upwind (won't work yet 23/1/07)
coords
0 - x-y cartesian
1 - r-theta cylindrical
2 - r-z cylindrical

In the transformation the coordinates
x -> r
y -> theta/z
*/

/*There will be a lot of these - one for each boundary (4-6)
  for each component (f0,f1....B - 6) yielding something like
  32 separate boundary variables.
*/
//Definitions
#define ULINE  "----------------------------------------------------------"<<std::endl
#define T_ULINE  "**********************************************************"<<std::endl
#define TAB '\t'
//Color definitions - using ANSI codes - may not be portable

#ifndef IMPACTAWITHCOLOUR

#define RED "" // makes text red
#define GREEN "" // makes text green
#define YELLOW "" // makes text yellow
#define BLUE "" // makes text blue
#define PURPLE "" // makes text purple
#define CYAN "" // makes text cyan
#define WHITE "" // makes text white
#define BRED "" // makes text red
#define BGREEN "" // makes text green
#define BYELLOW "" // makes text yellow
#define BBLUE "" // makes text blue
#define BPURPLE "" // makes text purple
#define BCYAN "" // makes text cyan
#define BWHITE "" // makes text white
#define ENDFORMAT "" // ends text formatting

#else

#define RED "\033[0;"<<31<<"m" // makes text red
#define GREEN "\033[0;"<<32<<"m" // makes text green
#define YELLOW "\033[0;"<<33<<"m" // makes text yellow
#define BLUE "\033[0;"<<34<<"m" // makes text blue
#define PURPLE "\033[0;"<<35<<"m" // makes text purple
#define CYAN "\033[0;"<<36<<"m" // makes text cyan
#define WHITE "\033[0;"<<37<<"m" // makes text white
#define BRED "\033[1;"<<31<<"m" // makes text red
#define BGREEN "\033[1;"<<32<<"m" // makes text green
#define BYELLOW "\033[1;"<<33<<"m" // makes text yellow
#define BBLUE "\033[1;"<<34<<"m" // makes text blue
#define BPURPLE "\033[1;"<<35<<"m" // makes text purple
#define BCYAN "\033[1;"<<36<<"m" // makes text cyan
#define BWHITE "\033[1;"<<37<<"m" // makes text white
#define ENDFORMAT "\033[0m" // ends text formatting

#endif
/*
  jan 2011 - added switch-array, this can be switched on and off
*/

class IMPACTA_Switch_Array
{
 private:
  double *myarray;
  double *otherarray;
  int Nel;
 public:
  IMPACTA_Switch_Array()
    {
      myarray=new double[1];
      otherarray=new double[1];
      myarray[0]=0.0;
      otherarray[0]=0.0;
      Nel=1;
    }
    IMPACTA_Switch_Array(int numberofelements)
    {
      myarray=new double[numberofelements];
      otherarray=new double[numberofelements];
      Nel=numberofelements;

      for (int i=0;i<numberofelements;++i)
	{
	  myarray[i]=0.0;
	  otherarray[i]=0.0;
	}
    }
    IMPACTA_Switch_Array(int numberofelements,int value)
    {
      myarray=new double[numberofelements];
      otherarray=new double[numberofelements];
      Nel=numberofelements;

      for (int i=0;i<numberofelements;++i)
	{
	  myarray[i]=value;
	  otherarray[i]=0.0;
	}
    }
    //*************************************
    //Access methods
    //*************************************
    
    double get(int el) 
    {
      return myarray[el];
    }
    double operator[](int el)
    {
       return myarray[el];
    }
   void set(int el,int val) 
    {
      myarray[el]= (double) val;
    }
    void set(int el,double val) 
    {
      myarray[el]=val;
    }
   void inc(int el,double val) 
    {
      myarray[el]+=val;
    }
    void switch_off()
    {
      int summation=0;
      for (int i=0;i<Nel;++i)
	summation+=myarray[i];
      if (summation) {
	double *temp;
	temp=myarray;
	myarray=otherarray;
	otherarray=temp;
      }
    }
 void switch_on()
    {
      int summation=0;
      for (int i=0;i<Nel;++i)
	summation+=myarray[i];
      if (!summation) {
	double *temp;
	temp=myarray;
	myarray=otherarray;
	otherarray=temp;
      }
    }
};

namespace IMPACT_Boundaries
{
  // {xlower, xupper, ylower, yupper}

  IMPACTA_Switch_Array fix_f0(4);
  IMPACTA_Switch_Array fix_f1_E(4);
  IMPACTA_Switch_Array fix_f2(4);
  IMPACTA_Switch_Array fix_B(4);
  IMPACTA_Switch_Array fix_ni(4);
  IMPACTA_Switch_Array fix_Ci(4);
  IMPACTA_Switch_Array fixed_any(4);
  
  int ifboundBx=0;
  int ifboundBy=0;

  int boundary_type=0;
  
}
/* These are namespaces for warnings*/
namespace warnings
{
  const int IMPACT_nowarnings=0;
}

namespace IMPACTA_smoothing
{
  double smooth_kernel[5]={1.0,0.0,0.0,0.0,0.0};
}

namespace zerotolerance
{
  double zerothreshold=1.e-25;
  double laggedtolerance=1.e-12;
  
  //for the next three, see Petsc manual
  double KSP_atol= 1.e-50; // default 1.e-50
  double KSP_rtol= 1.e-4; //default 1.e-2
  double KSP_dtol= 1.e5; //default 1.e5
  int MatrixSolver_ItMax=5000; // max number of iterations default 10^5

  int linear_solution_PC=0;
  int explicit_PC=0;
  int on_first_lag_low_rtol=0;
  double low_rtol=1e-2;
  //  double dispj_off_frac=1e-12;
  
  // These are for converging the matrix efficiently
  double iterate_matrix=0.0;
  double init_matrix_it=1e-3; //this is the initial matrix iteration
  //value taken if matrix does not converge
  double matrix_it_multiplier=10.0; //what the iteration value changes by
  double iterate_matrix_orig_val=0.0;
  double minimum_it_matrix_val=1e-10; //to prevent accidental appearance of 0 on diagonal
  
  double picard_divergence=1.0;//whether picard iterations diverge or not.
  int max_picard_its = 100; // maximum number of picard iterations

  int adaptive_timesteps=1; // for adaptive timestep iteration
  int adaptive_multiplier=2;
  int adaptivetimes=5;// times the matrix solve should be good before
  int adaptivetimesoriginal=5;// the adaptive_timestep number reduces
   
  double equil_percent=1.0e-6;
  
  double RB_D_tolerance = 1e-6;
  double CD_tolerance = 1e-12;
  double CD_frac_tol = 1e-7;
  int RB_D_itmax = 50;
  int RB_iterate = 1;
}
namespace equation_switches
{
  //These switch whole equations and parts of equations on or off
  int f0_equation_on = 1;
  int evolvef0 = 0;
  int f1_equation_on = 1;
  int f2_equation_on = 0;
  int f3_equation_on = 0;
  int E_equation_on = 1;
  int B_equation_on = 1;
  //components of the f0eqn
  double df0_by_dt_on = 1.0;
  double inf0_grad_f1_on = 1.0;
  double inf0_Edf1dv_on = 1.0;
  double Cee0_on = 1.0;
  //f1 eqn
  double e_inert_on=1.0;
  double inf1_vgradf0_on = 1.0;
  double inf1_Edf0dv_on = 1.0;
  double inf1_f1xB_on = 1.0;
  double Cei_on = 1.0;
  int Cee1_on = 1;
  //E eqn
  double disp_j_on = 1.0;
  //
  double e_viscosity_on=1.0;
  double df3_by_dt_on=1.0; // This is actually now for Weibel damping f3
  //
  int relaxtoeq=1;
  //
  double Bimp_in_E_equation=1.0; //rememver - if zero affects dEbydtfix
  double jimp_in_E_equation=1.0; //remember - can't be zero
  double Eimp_in_B_equation=1.0;
  double dEbydtfix=0.5; //1e-8;
  double Gamma_A=2.0; // for initializing distribution

  // This helps fix matrix solve problems
  double NEW_T_NORM=1.0;



  // Ray tracing
  bool tr_bool;

}
namespace if_dump_switches
{
  int ifdumpall=1;
  int dump_ne=1;
  int dump_ni=0;
  int dump_Ci=0;
  int dump_Z=0;
  int dump_Te=1;
  int dump_Ue=1;
  int dump_je=1;
  int dump_q=1;
  int dump_P=1;
  int dump_Q=1;
  int dump_E=1;
  int dump_B=1;
  int dump_wt=1;
  int dump_VN=1;
  int dump_eta=1;
  int dump_f0=0;
  int dump_f1=0;
  int dump_f2=0;
  int f_ndump=1;
  int view_matrix=0;
}
/* Following are the values of global constants etc.*/
double Calculate_Lambda();
namespace globalconsts
{
  double Y=1.0;
  const double pi=3.141592653589793238;
  const double twopi=2.0*pi;
  const double fourpi=4.0*pi;
  const double fourthirdspi=4.0*pi/3.0;
  const double fourthirds=4.0/3.0;
  const double oneover3 = 1.0/3.0;
  const double oneover6 = 1.0/6.0;
  const double pibysix=pi/6.0;
  const double sqrt2 = sqrt(2.0);
  
  double oneover_atomic_Z = 0.1;
  double c_L = 10;//3e8; // speed of light over v electron thermal
  double epsilon0=1.0;//8.85e-12;
  double e_charge = -1.0; //1.6e-19 // electron charge.
  double e_mass = 1; //9.11e-31 //electron mass
  double me_over_mp = 0.000544617024;
  double omega_p_nuei_n = 10; //w0*ln
  const int LenName=20;
  const int Ncoords=2;
  char IMPACT_Coords[Ncoords]={'x','y'};
  double deltasquared=1.0;
  double Bimp=1.0;

  //  conversion of T in keV to normalized units
  const double keVnorm = 511.0; //to be multiplied by v^2/c^2


  // Actual value of density!
  // collision frequency in units of 1/s 
  double Lambda_ei = Calculate_Lambda();
  double logLambda = log( Lambda_ei); 
  double N_par_Debye =  Lambda_ei/(9.0*oneover_atomic_Z);
  double real_T_J = 9.11e-31/c_L/c_L*2.998e8*2.998e8; // real temperature in J
  double ne24 = 16.0/9.0*pi*pi*pow(8.85e-12/1.602e-19/1.602e-19,3)/N_par_Debye/N_par_Debye*real_T_J*real_T_J*real_T_J/1.0e6/1.0e24;
  double real_wp = sqrt(1.602e-19*1.602e-19*ne24*1.0e6*1.0e24/8.85e-12/9.11e-31);
  double nuei = sqrt(2.0/pi)*logLambda/ Lambda_ei*real_wp;
 
  // I don't take into account Z or T dependance as is logarithmic anyway, 
  // and calculation would really slow things down
  double ionize_coll_log = log(1e5/c_L/c_L*oneover_atomic_Z)/logLambda; // normalized collision logarithm from Bethe formula

  double laser_frequency = 1.78e15; // angular frequency of laser
  const double Z_minimum=1e-3; // minimum ionization state.
  
  int MPIChunkSize = 1;//10000;
  std::string initial_time_str="";

  // For input maths
const int nsymbols=5;
const int nfuncs=9;
std::string IMPACTA_symbols[nsymbols]={"+","-","*","/","^"};
std::string IMPACTA_funcs[nfuncs]={"(","COS(","SIN(","TAN(","EXP(",
				   "TANH(","SINH(","COSH(","LOG("};
//numbers as strings useful
const int nnumbers=16;
std::string IMPACTA_numbers[nnumbers]={"1","2","3","4","5","6","7","8","9","0"
				       ,".", "LX","LY","PI","X","Y"};
}
namespace IMPACTA_ions
{
  int ionization_on=0;
  bool quasin_on=false;
  std::string ionization_model[3] = {"off","Saha","Thomas-Fermi"};
  int ion_motion=0;
  // double a1bar=0.00544617024;
  //double a2bar=a1bar;
  double atomic_mass_number=10.0;
  double ion_temperature=1.0;
  double alpha_ion = (globalconsts::me_over_mp/atomic_mass_number
                      /(globalconsts::oneover_atomic_Z));

  double ionization_time_frequency = 4.11e16/globalconsts::nuei; // omega_a normalized to e-i collision time
  double ionization_v_osc = 0.168 *(1.78e15/globalconsts::laser_frequency)*globalconsts::c_L; // 0.168 is a0 for 1 micron laser.
  const int N_Z_smooth = 10; //number of steps between smoothing of Z
}

namespace VISIT_PARAMS
{
  const int data_id0=0;//NONE
  const int data_id1=1;//ne
  const int data_id2=2;//ni
  const int data_id3=3;//Te
  const int data_id4=4;//Bz
  const int data_id5=5;
  const int data_id6=6;
  const int data_id7=7;
  const int data_id8=8;

  const int vec_id1=0; //E field
  const int vec_id2=1; //B field
  const int vec_id3=2; // j
  const int vec_id4=3; // q
}
namespace IMPACT_Diagnostics
{
  double init_start_time=0.0,init_end_time=0.0;//to measure initialization time
  double n_start_time=0.0,n_end_time=0.0; //to measure timestep CPU time
  double start_time=0.0,end_time=0.0; //to measure total CPU time
  double total_nl_its=0.0,nl_times=0.0; //to measure ave nonlinear iterations
  double total_picard_its=0.0,picard_times=0.0; //to measure 
  const int graphbars=50;

  double timegraph[graphbars+1];
  double total_delta_its=0.0,delta_times=0.0; //to measure delta in RB coeffs
  int showmatrix=0;
  int showvector=0;
  int noxnbody=0;
  int divcheck=1; // Check for diveregnce in Picard iteration

  int output_precision=10; //number of digits precision for output
}
std::string IMPACT_GetTime(double time)
{
  std::ostringstream outtime;
  int days,hours,mins;
  double secs=time;
  days=int (secs)/86400;
  secs-=86400.0*days;
  hours=int (secs)/3600;
  secs-=3600.0*hours;
  mins=int (secs)/60;
  secs-=60.0*mins;
if (days==1) outtime<<days<<" day, ";
  if (hours==1) outtime<<hours<<" hour, ";
  if (mins==1) outtime<<mins<<" min, ";
  if (days>1) outtime<<days<<" days, ";
  if (hours>1) outtime<<hours<<" hours, ";
  if (mins>1) outtime<<mins<<" mins, ";
  outtime<<secs<<" s";
 return outtime.str();
}
//Using this namespace means that space will be saved in
//the matrix for the components present in Robs version
//of IMPACT.
namespace IMPACT_Original
{
  int f0_on=1,f1_on=1,f2_on=0,f3_on=0,E_on=1,num_B=3;
}
namespace IMPACT_Messages
{
  char Version[20]="Wolfhound Alpha";
  char Date[11]="01/06/2012";
  char DataHeader[2][200]= {"IMPACTA Moment Data file (.imd)\nAGRT07\nNotes:\nThe data can be returned into an IMPACT_Moment structure\nby using the << operator.","IMPACTA Distribution Data file (.imd)\nAGRT07\nNotes:\nThe data can be returned into an IMPACT_Dist class\nby using the << operator."};
  std::string Root_Dir = "";
  std::string Data_Directory=Root_Dir+"impacta_data/";
  std::string Input_Directory=Root_Dir+"impacta_data_in/";
  std::string Field_Dir="FLD/";
  std::string DistFunc_Dir="DST/";
  std::string Moment_Dir="MNT/";
  std::string Constants_Dir="ION/";
  std::string InDeck=Root_Dir+"imstdin";
  
  int if_show_function_input=0;
}

namespace IMPACT_Input_Deck
{
  const int Extra_cmdln=6;
  const int MAXVAR=136; //maximum number of variable names for input deck +1
  const int MAXLEN = 22; //maximum length of variable name
 
  int Var_check[MAXVAR];
  char Var_names[MAXVAR][MAXLEN]=
    {
      "OUT_DIR",        //1
      "IN_DIR",
      "OUT_DIGS",      
      "N_MAX",
      "DT",
      "NDUMP",
      "T0",
      "NX",
      "NY",
      "NV",            //10
      "COORDINATES",
      "XMIN",          
      "YMIN",
      "DXGRID",
      "DYGRID",    
      "DVGRID",
      "NUM_B",
      "E_ON",
      "F1_ON",
      "F2_ON",        // 20
      "F3_ON",      
      "ZERO_TOL",
      "LAG_TOL",
      "KSP_ATOL", "KSP_RTOL","KSP_DTOL",
      "KSP_MAX_IT","KSP_PC_TYPE","KSP_METHOD","KSP_MONITOR", //30
      "A","Z0","TI",
      "Z(X,Y)","Z(X)", "Z(Y)",  // 36
      "NE(X,Y)","NE(X)", "NE(Y)",
      "NI(X,Y)","NI(X)", "NI(Y)",   // 42
      "TE(X,Y)","TE(X)", "TE(Y)",
      "BX(X,Y)","BY(X,Y)","BZ(X,Y)", // 48
      "B0HAT","B0(X)", "B0(Y)",
      "W_PE_OVER_NU_EI", 
      "V_TE_OVER_C",
      "F0_EQUATION",
      "F1_EQUATION",        //55
      "F2_EQUATION",
      "F3_EQUATION",
      "E_EQUATION",
      "B_EQUATION",
      "DF0/DT", "INF0_GRAD_F1","INF0_EDF1DV","CEE0",  // 63
      "DF1/DT","INF1_VGRADF0","INF1_EDF0DV","INF1_F1XB","CEI","CEE1",  // 69
      "DF2/DT","DAMP_WEIBEL","DISP_J",                // 72
      "INITIAL_CONDITION","INITIAL_STEP_DT_FRAC",  // 74
      "RB_D_IT_TOL", "RB_D_IT_MAX", "RB_ITERATE",   // 77
      "HEAT_SOURCE","HEATING(X,Y)","HEATING(X)","HEATING(Y)",  //80
      "HEATING(T)","VOSC","POLARIZATION",  //83
      "GRADN_Z(X,Y)/N","GRADN_Z(X)/N","GRADN_Z(Y)/N","GRADN_Z(T)/N", //87
      "GRADT_Z(X,Y)/T","GRADT_Z(X)/T","GRADT_Z(Y)/T","GRADT_Z(T)/T",  //91
      "IF_DUMP_NE","IF_DUMP_NI","IF_DUMP_CI","IF_DUMP_Z", //95
      "IF_DUMP_TE","IF_DUMP_UE","IF_DUMP_JE","IF_DUMP_Q","IF_DUMP_P",    //100
      "IF_DUMP_E","IF_DUMP_B","IF_DUMP_WT","IF_DUMP_VN","IF_DUMP_ETA",   // 105
      "IF_DUMP_F0","IF_DUMP_F1","IF_DUMP_F2","F_NDUMP",  //109
      "X_BOUND","Y_BOUND", //111
      "FIX_F0","FIX_F1_E","FIX_F2","FIX_B","FIX_NI","FIX_CI",         // 117
      "SMOOTHING_KERNEL","IONIZATION_MODEL","QUASINEUTRALITY",      // 120
      "TR_BOOL","BEAM_WIDTH","BEAM_RES","BOUNDARY","SHAPE",
      "N_C","X_0","Y_0","THETA","DTHETA",  // 130
      "DS","INTENSITY_DIFF_STEPS","TEMP_MULTIPLIER","FRACTION",
      "END"}; // END goes at the end..
}

//________________________________________________________________
//   Useful functions
//
//This function compares the given string with the commandline arguments
// and returns true(1) or false(0) 
int IMPACT_Cmp_CmdLine(int argnum, char **argstr, std::string string)
{
  int result=0;
  char null[3]="";
  for (int i = 0; i<argnum;++i)
    if (!strcmp(argstr[i],string.c_str())) 
      {
	result = 1;
	argstr[i]=null;
      }
  
  return result;
}
std::string IMPACT_Get_CmdLine(int *argnum, char **argstr, std::string string)
{
  int dircheck=0;
  char null[3]="";
  std::string answer;
  for (int i=0;i<*argnum;++i)
    if (!strcmp(argstr[i], string.c_str())) {dircheck=i;argstr[i]=null;}
  if (dircheck>0)
    {answer= argstr[dircheck+1]; argstr[dircheck+1]=null;}
  else answer = "*nofield*";
  return answer;
}
// list of known functions
int IMPACT_isFunc(char * func)
{
  
  char func3[4]="";
  char func4[5]="";
  for (int i=0;i<5;++i)
    func[i]=toupper(func[i]);
  for (int i=0;i<3;++i)
    func3[i]=toupper(func[i+1]);
  for (int i=0;i<4;++i)
    func4[i]=toupper(func[i+1]);
  int result =0;
  //note the numbers correspond to dofunction
  if(!strcmp(func4,"TANH")) result = 6;
  if (!result){
  if(!strcmp(func3,"SIN")) result = 1;
  if(!strcmp(func3,"COS")) result = 2;
  if(!strcmp(func3,"TAN")) result = 3;
  if(!strcmp(func3,"EXP")) result = 4;
  if(!strcmp(func3,"LOG")) result = 5;
  }
  return result;
}
double IMPACT_function(double value,int result)
{
  double answer=0.0;
  switch (result)
    {
    case 1:
      answer = sin(value);
      break;
    case 2:
     answer = cos(value);
      break;
    case 3:
      answer = tan(value);
      break;
    case 4:
      answer = exp(value);
      break;
    case 5:
      answer = log(value);
      break;
    case 6:
      answer = tanh(value);
      break;
    }
  return answer;
}
//executes the given function
void IMPACT_dofunction(double *gridtemp,int N,int result)
{

  for (int i=0;i<N;++i) gridtemp[i] = IMPACT_function(gridtemp[i],result);
}

double IMPACT_isConst(char* number, double max,int N)
{
  int length=strlen(number);
  std::istringstream numtemp;
  numtemp.str(number);			    
  for (int i=0;i<length;++i)
    number[i]=toupper(number[i]);
  double result =0;
  if(!strcmp(number,"PI")) result = globalconsts::pi;
  if (!strcmp(number,"MAX")) 
    {
      if (max==0.0) {std::cout<< "IMPACT: ERROR - in input deck, max requested but max - min = 0.0\n"; exit(0);} 
      result = max*(double)N/(double)(N-1);
    }
  if (!result) numtemp>>result;
  return result;
}
double IMPACT_isConst2D(char* number, double imax,double jmax,int Nx,int Ny)
{
  int length=strlen(number);
  std::istringstream numtemp;
  numtemp.str(number);			    
  for (int i=0;i<length;++i)
    number[i]=toupper(number[i]);
  double result =0;
  if(!strcmp(number,"PI")) result = globalconsts::pi;
  if (!strcmp(number,"ILIM")) 
    {
      if (imax==0.0) {std::cout<< "IMPACT: ERROR - in input deck, max requested but max - min = 0.0\n"; exit(0);} 
      result = imax*(double)Nx/(double)(Nx-1);
    }
if (!strcmp(number,"JLIM")) 
    {
      if (jmax==0.0) {std::cout<< "IMPACT: ERROR - in input deck, max requested but max - min = 0.0\n"; exit(0);} 
      result = jmax*(double)Ny/(double)(Ny-1);
    }
  if (!result) numtemp>>result;
  return result;
}
		
void chk()
{
  std::cout<<"EXIT OK"<<std::endl;exit(0);}
void chk(int check)
{
  std::cout<<"EXIT OK - "<< check<<std::endl;exit(0);}
void chk(double check)
{
  std::cout<<"EXIT OK - "<< check<<std::endl;exit(0);}
void chkns(int check)
{
  std::cout<<"OK - "<< check<<std::endl;}
void chkns(double check)
{
  std::cout<<"OK - "<< check<<std::endl;}

void chkMPI()
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (rank==0) std::cout<<"EXIT OK"<<std::endl; 
  MPI_Finalize(); 
  exit(0); 
}

/* for (int i = 31; i <= 37; i++)
    {
      std::cout << "\033[0;" << i << "mHello!\033[0m" << std::endl;
      std::cout << "\033[1;" << i << "mHello!\033[0m" << std::endl;
      }*/
  //_____________________________________________
 //Error handling - exit if alocation is outside range
 namespace no_vec_err_checking
{
  int vecerrchkon=0;
  inline void Vec_check_err(int *element, int *size)
  {
  }
}
  namespace with_vec_err_checking
{
  int vecerrchkon=1;
  inline void Vec_check_err(int *element,int *size)
  {
    
    try {if (*element>*size || *element<1)  throw 10;}
    catch (int err)
      {if (err==10)
	  {std::cout << std::endl<<"Exiting - Tried to access non-existent vector element";
	    std::cout<<std::endl;
	    std::cout<<"element="<<*element<<std::endl;
	    exit(0);}
      }
  }
}
