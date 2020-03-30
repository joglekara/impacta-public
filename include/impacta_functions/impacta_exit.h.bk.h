/*
**********************************************************
Exit and cleanup.

Version 1.0.0
AGRT

5/3/07

**********************************************************
*/
int IMPACT_Exit(IMPACT_MPI_Config *M,int status)
{
/*
     Always call PetscFinalize() before exiting a program.  This routine
     - finalizes the PETSc libraries as well as MPI
     - provides summary and diagnostic information if certain runtime
     options are chosen (e.g., -log_summary). 
   */
   MPI::COMM_WORLD.Barrier();
  if (status&&!M->rank())
    {
      double totaltime=(IMPACT_Diagnostics::end_time
			-IMPACT_Diagnostics::start_time)/CLOCKS_PER_SEC;
      double averagetime=IMPACT_Diagnostics::n_end_time/CLOCKS_PER_SEC;
      double averageittime=(IMPACT_Diagnostics::n_end_time/CLOCKS_PER_SEC
			    *IMPACT_Diagnostics::picard_times
			    /IMPACT_Diagnostics::total_picard_its);
      std::cout<<BBLUE<<ULINE<<BWHITE<<"  IMPACT finished successfully:\n\n"
	       <<CYAN
	       <<"  Total simulation time: "<<IMPACT_GetTime(totaltime)<<'\n'
	       <<"  Timestep average time: "<<IMPACT_GetTime(averagetime)<<'\n'
	       <<"  Lagged iteration average time: "<<IMPACT_GetTime(averageittime)
      
	       <<"\n  Average no. non-linear matrix iterations: "<<IMPACT_Diagnostics::total_nl_its/IMPACT_Diagnostics::nl_times<<"\n  Average no. lagged (Picard) iterations: "<<IMPACT_Diagnostics::total_picard_its/IMPACT_Diagnostics::picard_times<<'\n'
	   <<BBLUE<<ULINE<<'\n'<<ENDFORMAT;
    }
  
  //PetscErrorCode ierr=0;
  //ierr = PetscFinalize();CHKERRQ(ierr);

  if (M->size()>1)
    MPI::Finalize();
  
  exit(0);
  return 0;
  }
int IMPACT_Bad_Exit(int petsc_on)
{
  
  MPI::COMM_WORLD.Barrier();
  int size=MPI::COMM_WORLD.Get_size();
  //PetscErrorCode ierr=0;
  //if (petsc_on) ierr = PetscFinalize();CHKERRQ(ierr);
  
  if (size>1) MPI::Finalize();
  exit(0);
  return 0;
}
void IMPACT_wait ( int seconds )
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}
  
int IMPACT_Soft_Exit(IMPACT_MPI_Config *M,int reason)
{
  
  MPI::COMM_WORLD.Barrier();
  int result=1;
	int *allresults;
	allresults = new int[M->size()];
  std::string exit_reason="";
  switch (reason)
    {
    case 0:
      exit_reason="User Requested";
      break;
    case 1:
      exit_reason="Matrix Solver Error";
    }
  std::string unix_com = "test -d "+IMPACT_Messages::Root_Dir+"kill_IMPACT";
  result = system(unix_com.c_str()); 
  //std::cout<<"\nResult="<<result<<'\n';
  for (int i=0;i<M->size();++i)
    allresults[i]=1;  


  MPI::COMM_WORLD.Barrier();
  MPI::COMM_WORLD.Allgather(&result,1,MPI::INT,&allresults[0],1,MPI::INT);

  result=1;
  
  for (int i=0;i<M->size();++i)
    {
      //std::cout<<"Rank: "<<M->rank()<<", source: "<<i<<", value="<<allresults[i]<<'\n';
      if (!allresults[i]) result=0;
    }  
  if (!result) //Stop the simulation
    { 
      //  PetscErrorCode ierr;
      // ierr = PetscFinalize();CHKERRQ(ierr); //close petsc

      if (M->size()>1)
	{
	  MPI::Finalize();// close mpi
	  std::cout<<ENDFORMAT<<"Rank: "<<M->rank()<<" Finalize OK\n";   
	}  
      if (!M->rank())
	{
	  unix_com = "rm -rf "+IMPACT_Messages::Root_Dir+"kill_IMPACT";
	  result = system(unix_com.c_str()); 
	  std::cout<<BGREEN<<ULINE<<"\tIMPACT: Exited because "<<CYAN
		   <<exit_reason<<"\n"<<BGREEN<<ULINE<<'\n'<<ENDFORMAT;
	} 
      
      exit(0);
    }
    if (reason==1) exit(0);
  return 0;
  }


