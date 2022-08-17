/*
**********************************************************
For checking the closeness of the lagged vector to the previous
iteration

Version 2.9.4

AGRT

2/4/07 
march 08 - various updates including now outputting max value properly

jan 2011 - got rid of adaptive timestep
**********************************************************
*/
void IMPACT_File_Message(IMPACT_Config *c,int *n,int lagged, double maxvalue);

int IMPACT_Lag_Check(IMPACT_Config *c,IMPACT_MPI_Config* MPIc,IMPACT_ParVec *vec1, IMPACT_ParVec* vec2, int *n,int lag)
{
  double lagcheck = zerotolerance::laggedtolerance;
  double maxvalue = 0.0; //the largest value in the vector
  double value=0.0;
  int result = 0; //whether the lagcheck is successful or not
  // vec1->PrintArray();
  //vec2->PrintArray();
  for (int i=MPIc->start();i<=MPIc->end();++i)
    {
      value = (vec1->Get(i)-vec2->Get(i));
      if (fabs(value)>maxvalue) maxvalue = value;
    }
  if (fabs(maxvalue)<lagcheck) result = 1;
   //now make sure this is true on all processes
  int whichproc=0;
  if (MPIc->size()>1)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      int *all_results;
	all_results=new int [MPIc->size()];
      MPI_Allgather ( &result, 1, MPI_INT,all_results, 1, MPI_INT,
		      MPI_COMM_WORLD);
      result = 0;
      for (int i=0;i<MPIc->size();++i)
	result += all_results[i];
      result/=MPIc->size(); //integer divide so it will be zero unless all are 1
      
      MPI_Barrier(MPI_COMM_WORLD);
      //now get all maxvalues
      double *all_maxvalues;
      all_maxvalues = new double [MPIc->size()];
      MPI_Allgather(&maxvalue, 1, MPI_DOUBLE,all_maxvalues, 1, MPI_DOUBLE,
		    MPI_COMM_WORLD);
      for (int i=0;i<MPIc->size();++i)
	if (all_maxvalues[i]>maxvalue){
	  maxvalue = all_maxvalues[i];
	  whichproc=i;
	}
      delete[] all_maxvalues;
      delete[] all_results;
    }
  std::ostringstream Imessage;
  Imessage<<BCYAN<<"\n\n ******** \t l = "<<lag<<", max diff = "<<maxvalue<<" on rank "<<whichproc<<" (tol = "<<lagcheck<<")\n\n"<<ENDFORMAT;
  if (!MPIc->rank()) 
    {
      std::cout<<Imessage.str();
      IMPACT_File_Message(c,n,lag,maxvalue);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  // for fast convergence of the matrix
  double multiplier = floor(log10(fabs(maxvalue)*1e5));
  multiplier = pow(10,multiplier);
  if (multiplier>1.0) multiplier=1.0;
  zerotolerance::iterate_matrix*=multiplier;

  //check for divergence
  if (maxvalue>zerotolerance::picard_divergence && *n>0 && lag>1 &&IMPACT_Diagnostics::divcheck)
    {
      if (!MPIc->rank()) std::cout<<"\nIMPACTA: Divergence detected in Picard iteration\nExiting\n";
      exit(0);
    }
  return result;
}
