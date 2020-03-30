/*
  Ver 3
  Functions for the main program
  updated 2011 Jan -AGRT - removed adaptive timestep, as it was not 
  reliabliy implemented. General cleanup also
 */
  //Ray Tracing
void tracerunner(IMPACT_Matrix &n_mat, IMPACT_Matrix &Bz_mat, IMPACT_StenOps *stenops, IMPACT_Config *c, IMPACT_MPI_Config *M);
void IMPACTA_Share_Moment(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_IMatrix *Moment);


//This function simply creates a larger command line character string
char** IMPACT_Extend_Args(char **argstr_in,int *argnum)
{
  char **argstr;
    argstr = new char*[*argnum+IMPACT_Input_Deck::Extra_cmdln];
  for (int i=0;i<*argnum;++i)
    argstr[i]=argstr_in[i];
  for (int i=*argnum;i<*argnum+IMPACT_Input_Deck::Extra_cmdln;++i)
    {    
      argstr[i]=new char[30];
      for (int j=0;j<30;++j)
	argstr[i][j]=0;
    }
  *argnum= *argnum+IMPACT_Input_Deck::Extra_cmdln;
  return argstr;
}



int IMPACT_Lagged_Loop(IMPACT_Config *config1, IMPACT_MPI_Config *MPIconf,
		       IMPACT_ParSpa *S, IMPACT_StenOps *stenops,
		        int *argnum,char** argstr,
		       IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,
		       IMPACT_Var *f3,IMPACT_Var *E,IMPACT_Var *B,
		       IMPACT_ParVec *allvar,IMPACT_ParVec* allvar_lagged,
		       IMPACT_ParVec *allvar_lagged_old, 
		       int *n,int *lagged)
{
  MPI::COMM_WORLD.Barrier();
  IMPACT_Soft_Exit(MPIconf,0);
  
  IMPACT_Cee0(config1, MPIconf,allvar_lagged,f0,*lagged);
  
  //____________________________________________________________
  S->reset(); //clean sparse matrix
  
  // This function modifies E in z to allow magnetic field generation in the plane
  IMPACT_update_Ez(config1,MPIconf,f0,E,allvar_lagged,n);
  
  //Switch on or off elements (if d/dt ->0)
  if (!(*lagged))
    {

      // For fixed boundary conditions, update known vector
      //IMPACT_Fixed_Bounds_Sub_Vector(config1,MPIconf,allvar,f0,B);
    
      if (zerotolerance::linear_solution_PC)
      IMPACT_Linear_Soln(allvar_lagged,config1,MPIconf,stenops,f0,f1,f2,f3,E,B);
      SwapGhosts(config1,MPIconf, allvar);
    }

  //____________________________________________________________
  //Form sparse matrix from lagged parallel vector
  IMPACT_Form_Sparse(config1,MPIconf,allvar_lagged,S,
  		     stenops, f0,f1,f2,f3,E,B);
  
  //____________________________________________________________
  //Diagnose Matrix
  if(IMPACT_Diagnostics::showmatrix&&MPIconf->size())
    {
      IMPACT_Matrix A=S2M(S);
      std::cout<<"Matrix:\n";
      A.Print();
    }
  if(IMPACT_Diagnostics::showvector&&MPIconf->size())
    {
      std::cout<<"RHS vector:\n";
      allvar->PrintArray();
    }
  int checkmatrix=0,check=0;
  //____________________________________________________________
  //Solve sparse matrix
  if (equation_switches::f0_equation_on||equation_switches::f1_equation_on
      ||equation_switches::f2_equation_on||equation_switches::E_equation_on
      ||equation_switches::B_equation_on) {
    
    checkmatrix=MatrixSolve(config1,MPIconf,S,allvar,allvar_lagged,
			    *argnum,argstr,lagged);
  }
  else { 
    allvar_lagged->Duplicate(allvar);
  }

  if(!checkmatrix)
    {
      //____________________________________________________________
      // Swap the local boundary cells with ghost cells in the new vectors
      SwapGhosts(config1,MPIconf, allvar_lagged);

      //____________________________________________________________
      // Test how close the old and new solutions are to each other
      // (returns 1 if they are close)
      
      check = IMPACT_Lag_Check(config1,MPIconf, allvar_lagged,
			       allvar_lagged_old,n,(*lagged)+1);
      
      //____________________________________________________________
      
      *lagged=*lagged+1;
   
    }
  else 
    {
      if (!MPIconf->rank()) std::cout<<"\nIMPACTA: Matrix solver error\nExiting\n";
      exit(0);
    }
  allvar_lagged_old->Duplicate(allvar_lagged);
  return check;
}

int IMPACT_Main_Loop(IMPACT_Config *config1, IMPACT_MPI_Config *MPIconf,
		     IMPACT_ParSpa *S, IMPACT_StenOps *stenops,
		     int *argnum,char** argstr,
		     IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,
		     IMPACT_Var *f3,IMPACT_Var *E,IMPACT_Var *B,
		     IMPACT_ParVec *allvar,IMPACT_ParVec* allvar_lagged,
		     IMPACT_ParVec *allvar_lagged_old, int *n,const double dt)
{
  
  
  // reset matrix iteration value
  zerotolerance::iterate_matrix=zerotolerance::iterate_matrix_orig_val;
  //____________________________________________________________
  
  // Do ion motion
  if (IMPACTA_ions::ion_motion&&*n>0)
    {

      if (!MPIconf->rank()) std::cout<<"\nConvecting ion density....";
      IMPACTA_Update_ni(config1,MPIconf,allvar,stenops,f0,n);
      if (!MPIconf->rank()) std::cout<<"\nConvecting ionization state....";
      IMPACTA_Update_Z(config1,MPIconf,allvar,stenops,f0,n);
      if (!MPIconf->rank()) std::cout<<" Finished\nUpdating ion velocity....";
		IMPACTA_Update_Ci(config1,MPIconf,allvar,stenops,f0,f1,f2,E,B,n);
      if (!MPIconf->rank()) std::cout<<" Finished\n";
    }

  // Do ionization
  // New - set electron density from Z and ni if ionization on
  if (IMPACTA_ions::ionization_on&&*n>0) {
    IMPACTA_Ionization(config1, MPIconf,allvar,f0,*n-1);
  }
  
  //now update heating term temporal information:
  if (!MPIconf->rank())
    std::cout<<"Updating heating grid....";
//  !!!!!!!!!!!!!!! //////////////////
  //  !!!!!!!!!!!!!!! //////////////////
  //  !!!!!!!!!!!!!!! //////////////////
  //  !!!!!!!!!!!!!!! //////////////////
  //  !!!!!!!!!!!!!!! //////////////////
  // RAYTRACE STEP!!!!!!!!!!!!
  int timel,timeg;
   if (equation_switches::tr_bool){
        struct IMPACT_Moment ne,BM;
        IMPACT_Dim x1=3;

        IMPACT_Matrix n_mat(config1->Nx(),config1->Ny());
        IMPACT_IMatrix nemat,Bmat,localint;
    
        MPI::COMM_WORLD.Barrier();
        new_ne(allvar,f0,config1,MPIconf,&ne);
        nemat.copy(&ne.values);
        IMPACTA_Share_Moment(config1,MPIconf,&nemat);

        MPI::COMM_WORLD.Barrier();
        new_B(allvar,B,config1,MPIconf,&BM,&x1);
        Bmat.copy(&BM.values);
        IMPACTA_Share_Moment(config1,MPIconf,&Bmat); 
        
        timel = std::time(0);

        tracerunner(nemat,Bmat,stenops,config1,MPIconf);
        timel = std::time(0)-timel;

        MPI_Reduce(&timel,&timeg,1, MPI_DOUBLE, MPI_SUM,0,MPI::COMM_WORLD);
        timeg=timeg/MPIconf->size();

        if (!MPIconf->rank()) std::cout << "\n\n Average time for raytracing: " << timeg << " seconds \n\n";

        MPI::COMM_WORLD.Barrier();

      }
//  !!!!!!!!!!!!!!! //////////////////
//  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////
      //  !!!!!!!!!!!!!!! //////////////////

  IMPACT_update_heating(allvar,config1,MPIconf,stenops,f0,f1,E,*n-1);
  if (!MPIconf->rank())	std::cout<<" Finished\n";
  

  allvar_lagged->Duplicate(allvar);
  allvar_lagged_old->Duplicate(allvar);
 
  IMPACT_Switch_Elements(config1,MPIconf,allvar,
  			 stenops,f0,f1,f2,f3,E,B,*n);

  int lagged=0,check=0; //lagged loop index and tolerance condition flag
 
  lagged=0;
  //____________________________________________________________
  // Print time step number
  float number=*n;
  if (*n==-1) number = zerotolerance::equil_percent;
  if (*n==0) number = 2.0*zerotolerance::equil_percent;
  std::ostringstream Imessage;
  if (!MPIconf->rank())
    {
      Imessage<<'\n'<<BCYAN<<T_ULINE<<BBLUE;
      Imessage<<"\n ******** n = "<<number<<"\n\n"<<BCYAN<<T_ULINE<<'\n'<<ENDFORMAT;
      std::cout<<Imessage.str();
    }
  MPI::COMM_WORLD.Barrier();
  //____________________________________________________________
  // ******* START OF LAGGED LOOP *******
  
  do {
    check=IMPACT_Lagged_Loop(config1,MPIconf,S,stenops,argnum,
				 argstr,f0,f1,f2,f3,E,B,allvar,
			     allvar_lagged,allvar_lagged_old,
			     n,&lagged);
    MPI::COMM_WORLD.Barrier();
    
    // Check maximum picard iterations is not exceeded
    if (lagged>zerotolerance::max_picard_its&&zerotolerance::max_picard_its) {
      if (!MPIconf->rank())
	{
	  std::cout<<"IMPACTA: Exceeding maximum Picard iterations - Continuing to next timestep";
	}
      check = 1;
    }
  } while(!check); 
  
  // ******* END OF LAGGED LOOP *******
  //____________________________________________________________

 //for diagnostics....
 IMPACT_Diagnostics::total_picard_its+=lagged;
 IMPACT_Diagnostics::picard_times+=1.0;
 
 //Update the main vector to the next timestep
 allvar->Duplicate(allvar_lagged);
 SwapGhosts(config1,MPIconf, allvar);
 return check;  

}
void n_time_message(IMPACT_MPI_Config *MPIconf,int *n)
{
  IMPACT_Diagnostics::n_end_time=clock();
  if(!MPIconf->rank())
    {
      std::ostringstream Imessage;
      Imessage<<"\n"<<ENDFORMAT<<ULINE;
      Imessage<<"  IMPACT: n = "<<*n<<" timestep took - " << 
	IMPACT_GetTime((IMPACT_Diagnostics::n_end_time-
			IMPACT_Diagnostics::n_start_time)/CLOCKS_PER_SEC)
	       <<"\n"<<ENDFORMAT;
      Imessage<<ULINE<<ENDFORMAT;
      std::cout<<Imessage.str();
    }
  IMPACT_Diagnostics::n_start_time=clock();
}
void init_time_message(IMPACT_MPI_Config *MPIconf)
{
   if (!MPIconf->rank())
    {
      std::ostringstream Imessage;
      Imessage<<"\n"<<ENDFORMAT<<"IMPACT: Initialization took - "<<
	IMPACT_GetTime((IMPACT_Diagnostics::init_end_time-
			IMPACT_Diagnostics::init_start_time)/CLOCKS_PER_SEC)
	       <<"\n"<<ULINE<<'\n';
      std::cout<<Imessage.str();}
}
void IMPACT_Exit_Main(IMPACT_Config *config1,IMPACT_MPI_Config *MPIconf)
{
  IMPACT_Diagnostics::end_time=clock();
  IMPACT_Diagnostics::n_end_time=((IMPACT_Diagnostics::end_time 
				   -IMPACT_Diagnostics::start_time)
				  /config1->n_max());
  //**********************************************************
  //Exit Program.
  IMPACT_Exit(MPIconf,1);
}
