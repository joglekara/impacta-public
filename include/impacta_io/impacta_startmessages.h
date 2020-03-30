/*
**********************************************************
The header messages that print to screen

Version 1.1
AGRT

1/7/087 - pdated output
30/4/07
29/01/08 - now writes current step and time to file run_info
**********************************************************
*/
int IMPACTA_Sum_bounds(int i)
{
  int answer=0;
  answer+=IMPACT_Boundaries::fix_f0[i];
  answer+=IMPACT_Boundaries::fix_f1_E[i];
  answer+=IMPACT_Boundaries::fix_f2[i];
  answer+=IMPACT_Boundaries::fix_B[i];
  answer+=IMPACT_Boundaries::fix_ni[i];
  answer+=(IMPACT_Boundaries::fix_Ci[i]!=0.0);
  return answer;
}
std::string IMPACTA_str_bounds(int i)
{
  std::string answer="";
  if(IMPACT_Boundaries::fix_f0[i]) answer+="f0, ";
  if(IMPACT_Boundaries::fix_f1_E[i]) answer+="f1, E, ";
  if(IMPACT_Boundaries::fix_f2[i]) answer+="f2, ";
  if(IMPACT_Boundaries::fix_B[i]) answer+="B, ";
  if(IMPACT_Boundaries::fix_ni[i]) answer+="ni, ";
  if(IMPACT_Boundaries::fix_Ci[i]) answer+="Ci, ";
  return answer;
}
void IMPACT_Start_Messages()
{ 
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  // for wokring out bounds
  int ju=0,jl=0,iu=0,il=0;
  il=IMPACTA_Sum_bounds(0);
  iu=IMPACTA_Sum_bounds(1);
  jl=IMPACTA_Sum_bounds(2);
  ju=IMPACTA_Sum_bounds(3);

  // For outputting bounds info
  
  std::string ilstr="",iustr="",jlstr="",justr="";
  ilstr=IMPACTA_str_bounds(0);
  iustr=IMPACTA_str_bounds(1);
  jlstr=IMPACTA_str_bounds(2);
  justr=IMPACTA_str_bounds(3);

  //InDeck is initially (28/3/07) "imstdin" - defined in Environments
  if (!rank)
    {
      std::cout<<"\n\n"<<BGREEN<<ULINE<<"*\t\tIMPACTA "<<IMPACT_Messages::Version<<"\t\t\t *\n";
      std::cout<<"*\t\t"<<CYAN<<"AGRT University of Michigan\t\t"<<BGREEN<<" *\n*\t\t"<<PURPLE<<IMPACT_Messages::Date<<BGREEN<<"\t\t\t\t *\n"<<ULINE<<ENDFORMAT<<'\n';
    
      std::cout<< "Ionization is "<< onoroff(IMPACTA_ions::ionization_on);
if (IMPACTA_ions::ionization_on) {
	std::cout << ", model is "<<IMPACTA_ions::ionization_model[IMPACTA_ions::ionization_on]<<'\n';
} else {
      std::cout <<'\n';
}
std::cout<< "Quasineutral ion model is "<< onoroff(IMPACTA_ions::quasin_on)<<'\n';
	std::cout<< "Number of particles in a Debye Sphere is "<<  globalconsts::N_par_Debye<<"\n";
      std::cout<< "Normalizing electron density is "<<  ne24*1e24<<" cm^-3\n";
      std::cout<< "Normalizing electron temperature is "<<  real_T_J/1.602e-19/equation_switches::NEW_T_NORM<<" eV\n";
      std::cout<< "Normalizing mean-free-path is "<<  2.998e8/c_L/nuei/1.0e-6<<" microns\n";
//      std::cout<< "w_p is "<<  globalconsts::real_wp << "\n1";
//      std::cout<< "nu_ei is "<<  globalconsts::nuei << "\n";
      
      double logL_QM = logLambda*sqrt(36.0/(oneover_atomic_Z*oneover_atomic_Z*real_T_J/equation_switches::NEW_T_NORM/1.602e-19));
      if (logLambda<logL_QM) {
	logL_QM = logLambda;
      }
      std::cout<< "Normalizing Coulomb logarithm (inc. QM correction where appropriate) is "<< logL_QM<<"\n\n";
    
      std::cout<<BPURPLE<<"Boundary Conditions:\n";
      if ((!iu && !il) || !strcmp(boundaries[0],"periodic"))
	std::cout<<CYAN<<"x boundaries are "<<GREEN<<boundaries[0]<<'\n';
      else
	{
	  std::cout<<CYAN<<"x upper boundary is "<<GREEN;
	  std::cout<<boundaries[0];
	  if (iu) std::cout<< ", fixed for "<<iustr;
	  std::cout<<'\n';
	  std::cout<<CYAN<<"x lower boundary is "<<GREEN;
	  std::cout<<boundaries[0];
	  if (il) std::cout<< ", fixed for "<<ilstr;
	  std::cout<<'\n';
	}
      if ((!ju && !jl)|| !strcmp(boundaries[1],"periodic"))
	std::cout<<CYAN<<"y boundaries are "<<GREEN<<boundaries[1]<<'\n';
 else
	{
	  std::cout<<CYAN<<"y upper boundary is "<<GREEN;
	  std::cout<<boundaries[1];
	  if (ju) std::cout<< ", fixed for "<<justr;
	  std::cout<<'\n';
	  std::cout<<CYAN<<"y lower boundary is "<<GREEN;
	  std::cout<<boundaries[1];
	  if (jl) std::cout<< ", fixed for "<<jlstr;
	  std::cout<<'\n';
	}
      std::cout<<BPURPLE<<"\nDifferencing:";
      std::cout<<CYAN<<"\nx is "<<GREEN<< xdiff<<CYAN <<" differenced";
      std::cout<<CYAN<<"\ny is "<<GREEN<<ydiff<<CYAN <<" differenced\n";
      std::cout<<BPURPLE<<"\nDiagnostics:\n"<<ENDFORMAT;
      std::cout<<CYAN<<"Sparse row-checking is "<<GREEN<<onoroff(sparserowcheck)<<'\n';
      std::cout<<CYAN<<"System vector zero checking is "<<GREEN<<onoroff(withzerovectorchecking)<<"\n";
      std::cout<<CYAN<<"Vector element checking is "<<GREEN<<onoroff(vecerrchkon)<<"\n\n"<<ENDFORMAT;

      /* 
	 New bit: LAX stencil operator
      */
      std::cout<<CYAN<<"f1 and E, stencil operator is:\n";
      std::cout<<GREEN<<"\t\t"<<LAX.get(3)<<"\n\t"<<LAX.get(2)<<"\t"<<LAX.get(0)<<"\t"<<LAX.get(1)<<"\n\t\t"<<LAX.get(4)<<"\n\n";
    }
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  initial_time_str=asctime (timeinfo);
}
std::string IMPACTA_Time_Graph(int *n)
{
  std::string answer ="";
  int nmax=*n,nmin=*n-IMPACT_Diagnostics::graphbars;
  if (nmin<0) nmin=0;
  if (nmax<IMPACT_Diagnostics::graphbars) nmax=IMPACT_Diagnostics::graphbars;
  std::ostringstream nminstr,nmaxstr,nhalfstr;
  nminstr<<nmin;nhalfstr<<((nmax+nmin)/2);nmaxstr<<nmax;
  double max_time = 0.0;
  double time_graph[IMPACT_Diagnostics::graphbars];
  for (int i=0;i<IMPACT_Diagnostics::graphbars;++i) 
    time_graph[i]=(IMPACT_Diagnostics::timegraph[i+1]
		   -IMPACT_Diagnostics::timegraph[i]);
  for (int i=0;i<IMPACT_Diagnostics::graphbars;++i) 
    if (max_time<time_graph[i]) max_time=time_graph[i];
  double yaxis[6]={max_time,max_time*0.8,max_time*0.6,
		   max_time*0.4,max_time*0.2,0.0};
  
  for (int j=0;j<IMPACT_Diagnostics::graphbars;++j) 
    answer=answer+"_";
  for (int i=0;i<5;++i) 
    {
      answer=answer+" _ ";
      answer=answer+IMPACT_GetTime(yaxis[i]);
      if (i==3) answer=answer+" time/step";
      answer=answer+'\n';
       for (int j=0;j<IMPACT_Diagnostics::graphbars;++j) 
	 {
	   if (time_graph[j]>yaxis[i+1]) 
	     answer=answer+"|";
	   else answer=answer+" ";
	 }   
    }
  answer=answer+" _ 0.0 s\n";
  for (int j=0;j<IMPACT_Diagnostics::graphbars;++j) 
    answer=answer+"_";
  answer=answer+'\n';
  answer=answer+nminstr.str();
  for (int j=0;j<IMPACT_Diagnostics::graphbars/2-2;++j) 
    answer=answer+" ";
  answer=answer+nhalfstr.str();
  for (int j=0;j<IMPACT_Diagnostics::graphbars/2-2;++j) 
    answer=answer+" ";
  answer=answer+nmaxstr.str()+'\n';
  for (int j=0;j<IMPACT_Diagnostics::graphbars/2-2;++j) 
    answer=answer+" ";
  answer=answer+"n\n";
  return answer;
  
}
void IMPACT_File_Message(IMPACT_Config *c,int *n,int lagged, double maxvalue)
{
  time_t rawtime;
      struct tm * timeinfo;
 
      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      std::string dateStr=asctime (timeinfo);
      IMPACT_Diagnostics::end_time=clock();
      IMPACT_Diagnostics::n_end_time=((IMPACT_Diagnostics::end_time 
				       -IMPACT_Diagnostics::start_time)
				      /(*n));
      double averagetime=IMPACT_Diagnostics::n_end_time/CLOCKS_PER_SEC;
      
      if (*n<=IMPACT_Diagnostics::graphbars&&*n>=0) 
	{
	  IMPACT_Diagnostics::timegraph[0]= IMPACT_Diagnostics::start_time;
	  IMPACT_Diagnostics::timegraph[*n]=clock()/CLOCKS_PER_SEC;
	}
      if (*n>IMPACT_Diagnostics::graphbars&&lagged==1) 
	{
	  for (int j=0;j<IMPACT_Diagnostics::graphbars;++j) 
	    IMPACT_Diagnostics::timegraph[j]=IMPACT_Diagnostics::timegraph[j+1];
	  IMPACT_Diagnostics::timegraph[IMPACT_Diagnostics::graphbars]=clock()/CLOCKS_PER_SEC;
	}
     
	

      if (averagetime<0.0) averagetime=0.0;
      std::ofstream outfile;
      std::string dir;
      std::string str;
      dir = IMPACT_Messages::Data_Directory;
      str = "runinfo";
      std::string name = dir+str;
      outfile.open(name.c_str());
      if (!outfile) IMPACT_ferr(name.c_str());
      outfile<<"IMPACTA current run information\n";
      outfile<<"-------------------------------\n\n";
      outfile<<"Current timestep = "<<*n;
      outfile<<"\nCurrent iteration step = "<<lagged<<'\n';
      outfile<<"Current iteration difference/tolerance = "<<maxvalue<<"/"
	     <<zerotolerance::laggedtolerance<<'\n';
      outfile<<"\nAverage time per timestep = "<<IMPACT_GetTime(averagetime);
      outfile<<"\nPredicted total length of simulation = "
	     <<IMPACT_GetTime(averagetime*c->n_max())<<'\n';
      outfile<<"\nRun started: "<<initial_time_str<<'\n';
      outfile<<"Last update: "<<dateStr<<"\n\nPerformance History (time/timestep)\n";
      outfile<<IMPACTA_Time_Graph(n);
     
      outfile.close();
}
