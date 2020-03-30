/*
**********************************************************
Function for exchanging inner cells from one processor 
with ghosts from another.

Version 3.0
AGRT

21/3/07

Note, when it is established that it works (!) we can 
get rid of one of the data arrays

28/8/07 - problems with LAM on brighid so debugging has to commence.
5/10/07 - MPI_BARRIER positions critical! new ones put in
because Ghost cells were not swapping

15/10/07 - MPI problems - I think it can't pass such a large buffer,
so I'll have to pass it cell by cell.

8/08 - New put in MPI::Request which ensures all messages are sent/recieved

5/2010 - Changed send and receive to sendrev call - why didn't I do that before!?
**********************************************************
*/
inline int SwapGhosts(IMPACT_Config *c,IMPACT_MPI_Config *M, IMPACT_ParVec *v)
{
  MPI::COMM_WORLD.Barrier();
  if (M->size()>1)
    {
      
      double timestart,timeend; //for timing SwapGhosts.
      
      timestart=clock();
      int istart=M->start(),iend=M->end(); //start and end of vector domain with NG
      //(NG - no ghosts)
      
      int rank = M->rank();
      int last_proc=M->size()-1; //last processor rank
      int Npoints_per_i = c->cellpoints()*c->Ny(); //number of points in an i block
      int mylowerfriend=rank-1, myupperfriend=rank+1; //neighbouring processors
      if (mylowerfriend<0) mylowerfriend=last_proc;
      if (myupperfriend>last_proc) myupperfriend=0; // wrap 
      // So rank 0 wants to communicate with ranks size-1 and 1, rank 1 - 0&2 etc.
      
	//needs delete!
      double *datalower;
	datalower=new double[Npoints_per_i]; // double array for sending lower bound data.
      double *dataupper;
	dataupper=new double[Npoints_per_i]; // double array for sending upper bound data.
      
      

      //Lower bboundary first - I send out my boundary i block (not ghost) with
      // unique code....
      //_______________________________________________________________________
      
      for (int i=0;i<Npoints_per_i;++i)
	{
	datalower[i]=v->Get(i+istart);
	dataupper[i]=0.0;
	}
      MPI::COMM_WORLD.Barrier();
      
      MPI::COMM_WORLD.Sendrecv(&datalower[0],Npoints_per_i,MPI::DOUBLE,
			       mylowerfriend,300+rank,
			       &dataupper[0],Npoints_per_i,MPI::DOUBLE,
			       myupperfriend,300+myupperfriend);
      
      MPI::COMM_WORLD.Barrier();
      
      //now we put this data into our upper ghost cells....
      //  these start at iend+1
      for (int i=0;i<Npoints_per_i;++i)
	v->Set(i+iend+1,dataupper[i]);
      
      //_______________________________________________________________________
      // COMPLETED! ... now upper boundary.

      for (int i=0;i<Npoints_per_i;++i)
	{
	  dataupper[i]=v->Get(i+iend-Npoints_per_i+1); 
	  datalower[i]=0.0;
	}
      MPI::COMM_WORLD.Barrier();
      
      MPI::COMM_WORLD.Sendrecv(&dataupper[0],Npoints_per_i,MPI::DOUBLE,
			       myupperfriend,400+rank,
			       &datalower[0],Npoints_per_i,MPI::DOUBLE,
			       mylowerfriend,400+mylowerfriend);
      
      MPI::COMM_WORLD.Barrier();

      //now we put this data into our lower ghost cells....
      for (int i=0;i<Npoints_per_i;++i)
	v->Set(i+M->start_WG(),datalower[i]);
      //_______________________________________________________________________
      timeend=clock();
      MPI::COMM_WORLD.Barrier();
      // MESSAGES
      std::ostringstream Imess;
      Imess<<"\nIMPACT: Processor "<<rank<<" of "<<last_proc+1<<" - Transfer of ghosts took "<<IMPACT_GetTime(timeend-timestart)<<" s";
      std::cout<<Imess.str();
	delete[] datalower;
	delete[] dataupper; 
   }
  else
    {
      int iend=M->end(); //start and end of vector domain with NG
      int Npoints_per_i = c->cellpoints()*c->Ny(); //number of points in an i block
      for (int i=1;i<=Npoints_per_i;++i)
	v->Set(i+iend,v->Get(i));

      for (int i=1;i<=Npoints_per_i;++i)
	v->Set(i+M->start_WG()-1,v->Get(i+iend-Npoints_per_i));
   
    }
  MPI::COMM_WORLD.Barrier();
  return 0;
}
