/*
**********************************************************
Implicit Magnetized Plasma and Collisional Transport with Anisotropy
IMPACTA
Version Greyhound Beta
AGRT University of Michigan

7/4/07
The main program for the numerical model

27/7/07 - The main program has become quite untidy - I am 
going to clean it up by diverting the loop into functions

3/08 - Added a few steps so that the solver can equilibriate
a bit.
3/08 - Adaptive timestep griding added
14/7 xnbody funcs removed
**********************************************************
*/
#include "impacta_headers.h" //this contains all the necessary headers
int main(int argnum, char **argstr_in)
{

  try
    {

  //**********************************************************
  //                   INITIALIZATION
  //**********************************************************
  IMPACT_Diagnostics::init_start_time=clock();
  IMPACT_Help(argnum, argstr_in); // help for command line flags
 
  MPI_Init(&argnum, &argstr_in); // initialize MPI

/*
    Here, we sort out command line input so that more can be added
    This is so that the PETSC library options can easily
    be called from the command deck.
  */
  char **argstr(IMPACT_Extend_Args(argstr_in,&argnum));

  /*make new config object - this will contain all necessary 
  objects such as x,y,v grids, dt, etc.
  This is input from a file named in the environment.h file*/
  IMPACT_Config config1(IMPACT_Input(&argnum, argstr));



  //set up boundary conditions
  IMPACTA_initiate_bounds();
  /*
    This bit is to deal with the case of reflect Bz boundary conditions
    The boundary conditions Bz and f1/E are flipped
    So we need to reverse the fixed boundaries also
    AGRT 2011
  */ 
  IMPACT_Boundaries::ifboundBx = !strcmp(boundaries[0],"reflect_Bz");
  IMPACT_Boundaries::ifboundBy = !strcmp(boundaries[1],"reflect_Bz");
  
  // calculate real units 
  update_real_quantities();
  
  //Headers
  IMPACT_Start_Messages();

  double fixedtempBzf1E=0.0;
  if ( IMPACT_Boundaries::ifboundBx) {
    for (int i=0;i<2;++i) {
      fixedtempBzf1E=IMPACT_Boundaries::fix_f1_E[i];
      IMPACT_Boundaries::fix_f1_E.set(i,IMPACT_Boundaries::fix_B[i]);
      IMPACT_Boundaries::fix_B.set(i,fixedtempBzf1E);
    }
  }
  if ( IMPACT_Boundaries::ifboundBy) {
    for (int i=2;i<4;++i) {
      fixedtempBzf1E=IMPACT_Boundaries::fix_f1_E[i];
      IMPACT_Boundaries::fix_f1_E.set(i,IMPACT_Boundaries::fix_B[i]);
      IMPACT_Boundaries::fix_B.set(i,fixedtempBzf1E);
    }
  }

  // Now switch off elements for Sparse count + initial condition
  IMPACT_Boundaries::fix_f0.switch_off();
  IMPACT_Boundaries::fix_f1_E.switch_off();
  IMPACT_Boundaries::fix_f2.switch_off();
  IMPACT_Boundaries::fix_B.switch_off();
  IMPACT_Boundaries::fix_Ci.switch_off();
  IMPACT_Boundaries::fix_ni.switch_off();
  IMPACT_Boundaries::fixed_any.switch_off();

  //Configure MPI and create object containing MPI information
  IMPACT_MPI_Config MPIconf(IMPACT_MPI_Setup(&argnum, argstr,&config1));
 
   //Create new directory structure
  IMPACT_Tree(&config1);
 

  //Set up variable objects which provide references to positions
  //in vectors
  IMPACT_Var f0(config1);
  IMPACT_Var B(config1,'b');
  IMPACT_Var E(config1,'e');
  IMPACT_Var f1(config1,'1');
  IMPACT_Var f2(config1,'2');
  IMPACT_Var f3(config1,'3');


  //Create parallel vectors (ParVec) which contain all 
  // the necessary data points
  IMPACT_ParVec *allvar; //i.e. "all variables"
  allvar= new IMPACT_ParVec(&MPIconf);
  IMPACT_ParVec *allvar_lagged;   //lagged data for iterative solver 
  allvar_lagged= new IMPACT_ParVec(&MPIconf);
  IMPACT_ParVec *allvar_lagged_old;   //lagged data from previous iteration 
  allvar_lagged_old= new IMPACT_ParVec(&MPIconf);


  // New - set electron density from Z and ni if ionization on
  if (IMPACTA_ions::ionization_on) {
    for (int i=MPIconf.istart();i<=MPIconf.iend();++i) {
      for (int j=1;j<=config1.Ny();++j) {
	Initial_Conditions::ne.set(i,j,Initial_Conditions::ni.get(&i,&j)*Initial_Conditions::Z.get(&i,&j));
      }
    }
    // Now initialize variables from data input from file:
    Initialize(allvar,&config1,&MPIconf,&f0,&f1,&f2,&f3,&E,&B);
    
  } else {
    // Now initialize variables from data input from file:
    Initialize(allvar,&config1,&MPIconf,&f0,&f1,&f2,&f3,&E,&B);
    
    if (!IMPACT_Cmp_CmdLine(argnum, argstr, "-no_calc_ni"))
      {
	//Initial_Conditions::ni.SwitchOffConst(config1.Nx(),config1.Ny());
	//Initial_Conditions::ni.ChangeSize(,config1.Ny());
	for (int i=MPIconf.istart();i<=MPIconf.iend();++i)
	  for (int j=1;j<=config1.Ny();++j)
	    {
	      Initial_Conditions::ni.set(i,j,Local_ne(allvar,&f0,&config1,&i,&j)/Initial_Conditions::Z.get(&i,&j));
	    }
	MPI_Barrier(MPI_COMM_WORLD);
	IMPACTA_Share_Moment(&config1,&MPIconf, &Initial_Conditions::ni);
	Initialize(allvar,&config1,&MPIconf,&f0,&f1,&f2,&f3,&E,&B);
      }
  }


  // Change temperature normalization to help matrix convergence 
  IMPACTA_Rescale_T(&config1,&MPIconf,allvar,&f0,&f1,&f2,&f3,&E,&B,
		    equation_switches::NEW_T_NORM);
 
//Now swap ghost cells - very important!
  SwapGhosts(&config1,&MPIconf,allvar);

  allvar_lagged->Duplicate(allvar);
  allvar_lagged_old->Duplicate(allvar);
  //Create stencil operators needed
  // Create stencil operators (d/dx etc.)
  IMPACT_StenOps stenops(&config1);
  
  // For the e-e collision operator, this space stores the flux
  Cee0_flux_store::flux = new IMPACT_VS_Array(&config1,&MPIconf);

  //for the IB contribution to f2
  Cee0_flux_store::f2_IB = new IMPACT_VS_Array(&config1,&MPIconf); 

  //Create Sparse Matrix
  IMPACT_ParSpa *S;
  S=new IMPACT_ParSpa(&MPIconf);

  //**********************************************************
  /*  
      This bit counts non-zero elements for a vector with all elements=1
      This is so that the total non-zero elements can be counted in
      the general case - then counting does not need to be repeated

      In future, can probably be replaced by a formula that just works 
      out the correct elements?
  */
  IMPACT_ParVec *unitvec; 
  unitvec= new IMPACT_ParVec(&MPIconf);
  unitvec->unitary();   
 // Dump t=0
//  if (!equation_switches::relaxtoeq)
 //   IMPACT_Dump(&config1,&MPIconf,allvar,&f0,&f1,&f2,&f3,&E,&B,0);
  MPI_Barrier(MPI_COMM_WORLD);

  int f0equationtemp=equation_switches::f0_equation_on;
  int f1equationtemp=equation_switches::f1_equation_on;
  int f2equationtemp=equation_switches::f2_equation_on;

  // equation_switches::f0_equation_on=1;
  if (config1.Nf1()>0) equation_switches::f1_equation_on=1;
  if (config1.Nf2()>0) equation_switches::f2_equation_on=1;  
  IMPACT_Cee0(&config1, &MPIconf,allvar_lagged,&f0,0);
  // change ion motion elements so count is correct
  IMPACTA_Set_Ci(&config1,1.0);

  IMPACT_Form_Sparse_count(&config1,&MPIconf,unitvec,S,&stenops,
			   &f0,&f1,&f2,&f3,&E,&B);
  
  if (IMPACTA_ions::ion_motion)
  {
  //Reset to zero
  IMPACTA_Set_Ci(&config1,0.0);

  // Now set Ci from boundaries
  IMPACTA_Initiate_Ci(&config1);
  //Update divergence of Ci
  IMPACTA_Update_DivCi(&config1,&MPIconf,&stenops);
  }
  //equation_switches::f0_equation_on=f0equationtemp;
  equation_switches::f1_equation_on=f1equationtemp;
  equation_switches::f2_equation_on=f2equationtemp;
  S->reset();

  delete unitvec;

  //**********************************************************
  const double dt=config1.dt();
 
 int n=0;
  IMPACT_Diagnostics::init_end_time=clock();
  init_time_message(&MPIconf);
  //**********************************************************
  //                 END OF INITIALIZATION
  //**********************************************************
   

  //**********************************************************
  //                 INITIAL CONDITION
  //**********************************************************
   
 /* 
     These are temporary stores for the switches
     for the initial step at time t="0"
   */
 
  double einertiatemp=equation_switches::e_inert_on;
  double df2bydttemp=equation_switches::e_viscosity_on;
  double dispjtemp=equation_switches::disp_j_on; 
  
  // For zero timestep, do not use fixed bounds - need equilibrium
  // E and f1

   // small timestep to start as a fraction
   //**************************************************************
  // Ready for main loop - first do 2 steps with very small dt
  // To bring in f1 then f2 on the second
  int Bequationtemp=equation_switches::B_equation_on;
  MPI_Barrier(MPI_COMM_WORLD);
  if (equation_switches::relaxtoeq>0)
   {
    //Initially a step with somefraction of dt is performed..
    // to get fields correctly
    // Also - d/dt terms are turned off for this step.
     
     config1.set_dt(dt*zerotolerance::equil_percent);
     
     equation_switches::f0_equation_on=0;
     equation_switches::B_equation_on=0;
     // equation_switches::f2_equation_on=0;
     equation_switches::e_inert_on=0.0;
     equation_switches::e_viscosity_on=0.0;
     equation_switches::disp_j_on=zerotolerance::equil_percent; //jan 2011 agrt - turns out it just can't be 0 - boo

     n=-1;
     
     // To get fields and f1
     if (equation_switches::relaxtoeq==2) {
	   equation_switches::f2_equation_on=0;
	}	
 
     IMPACT_Main_Loop(&config1,&MPIconf,S,&stenops,&argnum,
		      argstr,&f0,&f1,&f2,&f3,&E,&B,allvar,
		      allvar_lagged,allvar_lagged_old,&n,dt*zerotolerance::equil_percent);

     n=0;
     //To get f2 
     if (equation_switches::relaxtoeq==1&&f2equationtemp)
       {
	 equation_switches::f2_equation_on=f2equationtemp;
	 equation_switches::e_viscosity_on=zerotolerance::equil_percent; //jan 2011 agrt - turns out it just can't be 0 - boo	 
	 IMPACT_Main_Loop(&config1,&MPIconf,S,&stenops,&argnum,
			  argstr,&f0,&f1,&f2,&f3,&E,&B,allvar,
			  allvar_lagged,allvar_lagged_old,&n,2.0*dt*zerotolerance::equil_percent);
	
	  }
   }
  MPI_Barrier(MPI_COMM_WORLD);
 // Finally reset switches for main loop
  config1.set_dt(dt);
  equation_switches::f0_equation_on=f0equationtemp;
  equation_switches::f2_equation_on=f2equationtemp;
  equation_switches::B_equation_on=Bequationtemp;
  equation_switches::e_inert_on=einertiatemp;
  equation_switches::e_viscosity_on=df2bydttemp;
  equation_switches::disp_j_on=dispjtemp;
  
  IMPACT_Boundaries::fix_f0.switch_on();
  IMPACT_Boundaries::fix_f1_E.switch_on();
  IMPACT_Boundaries::fix_f2.switch_on();
  IMPACT_Boundaries::fix_B.switch_on();
  IMPACT_Boundaries::fix_Ci.switch_on();
  IMPACT_Boundaries::fix_ni.switch_on();
  IMPACT_Boundaries::fixed_any.switch_on();
  
  IMPACT_Dump(&config1,&MPIconf,allvar,&f0,&f1,&f2,&f3,&E,&B,n);
  
  //**************************************************************
  //                        *MAIN LOOP*
  //**************************************************************
 
  //Set clocks
  IMPACT_Diagnostics::start_time=clock();
  for (int i=0;i<=IMPACT_Diagnostics::graphbars;++i) 
    IMPACT_Diagnostics::timegraph[i]=0.0;
  


  //start of loop

  for (n=1;n<=config1.n_max();++n)
    {      
      MPI_Barrier(MPI_COMM_WORLD);
      IMPACT_Main_Loop(&config1,&MPIconf,S,&stenops,&argnum,
		       argstr,&f0,&f1,&f2,&f3,&E,&B,allvar,
		       allvar_lagged,allvar_lagged_old,&n,dt);


      for (int i=MPIconf.istart();i<=MPIconf.iend();++i) 
      {
          for (int j=1;j<=config1.Ny();++j) 
          {
            for (int k=1;k<=config1.Nv();++k) 
            {
              // if (f0.get(allvar,&i,&j,&k)<0) f0.set(allvar,0.0,&i,&j,&k);
            }
          } 
      }

      
      IMPACT_Dump(&config1,&MPIconf,allvar,&f0,&f1,&f2,&f3,&E,&B,n);
    

      n_time_message(&MPIconf,&n);
    } //end of main loop
  
  IMPACT_Exit_Main(&config1,&MPIconf);
    }
  catch (std::bad_exception&)
    {
      std::cout<<"\nIMPACT: Cleaning up...\n";
      IMPACT_Bad_Exit(0);
    }
  return 0; 
}
