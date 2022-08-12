/*
**********************************************************
Some Parallel stuff for Impacta. Includes MPI_Config class,
 Parallel Vector class and Gather function

Version 1.3
AGRT

5/3/07
12/3 - Many changes! Most important - Gather for sending
vectors between processors.

13/307 - Important to get the right domain in i when calling start/end or 
start_WG etc. as to whether to include ghost cells or not.

Also ghost cells now include 0 and Nx+1 - i.e. beyond the bounds
of the simulation. This is to eventually deal with periodic boundary conditions.

10/05/07 - update - because the Gather function results in
unnessary memory usage, GatherMoments was written
this gathers the data from multiple processors and sends just
the moment data to rank 0
**********************************************************
*/

class IMPACT_MPI_Config
{
 private:
  int my_rank; //MPI rank of this machine
  int num_procs; //number of processors
  int my_i_lower, my_i_upper; //lower and upper i bounds WITHIN ghost boundaries
  //i.e. not including ghost cells.
  int my_i_lower_noghosts,my_i_upper_noghosts;
  int *vecupper;
  int *veclower; //upper and lower bounds of all processors WITH gc's
  int *iupper;
  int *ilower; //upper and lower bounds in i of all processors
  int myvecupper_wghosts;
  int myveclower_wghosts;//upper&lower with no ghost cells included
 public:
  ~IMPACT_MPI_Config()
    {
      delete[] vecupper; 
      delete[] veclower; 
      delete[] iupper; 
      delete[] ilower; 
      }
  IMPACT_MPI_Config()
    {
      my_rank=0;
      num_procs=1;
      my_i_lower= my_i_upper=1;
      my_i_upper_noghosts=1;
      my_i_lower_noghosts=1;
      vecupper=new int[1];
      veclower=new int[1];
      vecupper[0]=1;
      veclower[0]=1;
      myvecupper_wghosts=1;
      myveclower_wghosts=1;
      iupper=new int[1];
      ilower=new int[1];
    }
  IMPACT_MPI_Config(IMPACT_Config *c,int *xi,int *xf,int xi_WG, int xf_WG,int myrank,int numprocs)
    {
      my_rank=myrank;
      num_procs=numprocs;
      my_i_upper_noghosts=xf[myrank];
      my_i_lower_noghosts= xi[myrank];
      my_i_lower=xi_WG;
      my_i_upper=xf_WG;
      vecupper=new int[numprocs];
      veclower=new int[numprocs];
      for (int i=0;i<num_procs;++i)
	{
	  vecupper[i]=xf[i]*c->cellpoints()*c->Ny();//xf lots of Ny cells
	  veclower[i]=(xi[i]-1)*c->cellpoints()*c->Ny()+1;//xi lots of Ny cells
	}
      iupper=new int[numprocs];
      ilower=new int[numprocs];
      for (int i=0;i<num_procs;++i)
	{
	  iupper[i]=xf[i];//xf 
	  ilower[i]=xi[i];//xi 
	}
      myvecupper_wghosts=xf_WG*c->cellpoints()*c->Ny();
     myveclower_wghosts=(xi_WG-1)*c->cellpoints()*c->Ny()+1;
    }
  //access methods.
  int istart_WG() //WG- with ghosts
  {
    return my_i_lower;
  }
  int iend_WG()
  {
    return my_i_upper;
  }
 int istart()
  {
    return my_i_lower_noghosts;
  }
  int iend()
  {
    return my_i_upper_noghosts;
  }
  int start() 
  {
    return veclower[my_rank];
  }
  int end()
  {
    return vecupper[my_rank];
  }
  int start_WG() //WG - with ghosts
  {
    return myveclower_wghosts;
  }
 int end_WG()
  {
    return myvecupper_wghosts;
  }
  int start(int i)
  {
    return veclower[i];
  }
  int end(int i)
  {
    return vecupper[i];
  }
  int N()
  {
    return vecupper[my_rank]-veclower[my_rank]+1;
  }
  int N_WG() //size with ghosts 
  {
    return myvecupper_wghosts-myveclower_wghosts+1;
  }

  int rank()
  {						
    return my_rank;
  }
  int rank(int *i)
  {		
    int answer=0;
    for (int j=0;j<num_procs;++j)
      {
	if (*i>=ilower[j]&&*i<=iupper[j]) answer=j;
      }
    return answer;
  }
  int size()
  {
    return num_procs;
  }
};

//Configure MPI
IMPACT_MPI_Config IMPACT_MPI_Setup(int * argnum, char** argstr,IMPACT_Config *c)
{  //MPI Stuff
  int rank, size;
  int Numx=c->Nx();
  
  //size=MPI::COMM_WORLD.Get_size();
  //rank=MPI::COMM_WORLD.Get_rank();
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if (size>Numx)
    {
      std::cout<<"IMPACT: ERROR more processors than x gridpoints!"<<std::endl;
      exit(0);
    }
// changed to be ISO C++ compliant - needs delete
  int * xfirst;
  int * xlast; //the bounds of all domains.
  xfirst = new int[size];	
  xlast = new int[size];	
  
  int xfirst_WG=0,xlast_WG=0; //the bounds of this domain with no ghosts
  //First we work out all processor values and check for consistency
  //as well as assigning the above bounds.
  if (rank==0)
    {
  std::cout <<std::endl<<BCYAN<<"IMPACTA: Distributing i grid over "<<size<<" processor(s)"<<CYAN<<std::endl<<ULINE;
  std::cout<<"\tGrid Division\t\t|\tWith Ghosts"<<std::endl;
  std::cout <<"Rank\timin\timax\tNx\t|\timin\timax\tNx"<<std::endl<<ULINE<<ENDFORMAT;
    }
  int ranksumx=0; //to check consistency
  for(int ranktemp=0;ranktemp<size;ranktemp++)
    {

  //Now we must split the domain in i and j.
  int Nxtemp=Numx/size; //integer divide...

  int diff_x = Numx-size*Nxtemp; // find difference...
  int my_x_start=1; //start of local domain
  int myghoststart=1; //start with ghost cells
  int myghostend=1; //end ditto
  
  if (ranktemp<diff_x) 
    {
      ++Nxtemp;
      my_x_start=ranktemp*Nxtemp+1;
    }
  else
    {
      my_x_start=diff_x*(Nxtemp+1)+1;
      my_x_start+=Nxtemp*(ranktemp-diff_x);
    }  
/*Now we make the start and end points one larger for all processors not 
    at the endpoints in i and j - this is so that they can perform
    differentials without having to communicate with other processors
    Obviously unfortunately this is not the case for the edge processors
    
    13/3 - now it is the case, there are ghost cells for all processors.
  */
 myghoststart=my_x_start-1;
 myghostend=my_x_start+Nxtemp;

  if (rank==0)
  {  std::cout<<CYAN<<ranktemp<<"\t"<<my_x_start<<"\t "<<my_x_start+Nxtemp-1<<"\t"<<Nxtemp<<"\t|\t"<<myghoststart<<"\t"<<myghostend<<"\t"<<myghostend-myghoststart+1<<std::endl;}

  xfirst[ranktemp]=my_x_start;
  xlast[ranktemp]=my_x_start+Nxtemp-1;
  if (ranktemp==rank)
    {
      xfirst_WG=myghoststart; //points with ghost cells.
      xlast_WG=myghostend;
    }
      ranksumx+=Nxtemp;
    }
  if (rank==0) std::cout<<ULINE<<ENDFORMAT;
  
  //Check whether domain divide is self-consistent  
  if ((ranksumx-Numx)!=0&&rank==0)
    {
      std::cout<<"IMPACTA: ERROR - in IMPACT_MPI_Setup domain not correctly divided"<<std::endl<<std::endl;
      exit(0);
    }
  //
  
  IMPACT_MPI_Config MPIc(c,xfirst,xlast,xfirst_WG,xlast_WG,rank,size);

  delete[] xfirst;
  delete[] xlast;
  return MPIc;
}

/*
Parallel vector class - is basically a Vector but instead of 
taking the element given to it it calls the element in IMPACT_Vector
shifted by its start point which is to do with its ownership
 */
class IMPACT_ParVec : public IMPACT_Vector
{
 private:
  int mystart,myend;  
 public:
 IMPACT_ParVec() : IMPACT_Vector()
    {
      mystart=1;myend=1;
      reset();
    }
  //in the next constructor it is +1 because end-start=N-1
 IMPACT_ParVec(int startin,int endin) : IMPACT_Vector(endin-startin+1)
    {
      mystart=startin;
      myend=endin;
      reset();
    }
 IMPACT_ParVec(IMPACT_MPI_Config *MPIconf) : IMPACT_Vector(MPIconf->N_WG())
    {
      mystart=MPIconf->start_WG();
      myend=MPIconf->end_WG();
      reset();
    }

 double Get(int element)
 {
   //error check;
   int elemshift=element-mystart+1;
   // chkerror(elemshift);
   return IMPACT_Vector::Get(elemshift);
 }
double* Get_add(int element)
 {
   //error check;
   int elemshift=element-mystart+1;
   // chkerror(elemshift);
   return IMPACT_Vector::Get_add(elemshift);
 }
 void Set(int element, double value)
 {
   //error check;
   int elemshift=element-mystart+1;
   //  chkerror(elemshift);
   return IMPACT_Vector::Set(elemshift,value);
 }
 int length()
 {
   return IMPACT_Vector::length();
 }
int start()
 {
   return mystart;
 }
int end()
 {
   return myend;
 }
void ChangeSize(int array_size)
{
  IMPACT_Vector::ChangeSize(array_size);
}
void PrintArray()
{
  IMPACT_Vector::PrintArray();
}
void PrintNonZero()
{
  IMPACT_Vector::PrintNonZero();
}
 void unitary()
 {
   IMPACT_Vector::unitary();
 }
void reset()
 {
   IMPACT_Vector::reset();
 }
void operator=(IMPACT_Vector vec)
  {
    IMPACT_Vector::operator=(vec);
  }
  int chkzero(int element)
{
  int elemshift=element-mystart+1;
  //chkerror(elemshift);
 return IMPACT_Vector::chkzero(elemshift);
}
 void chkerror(int element)
 {
   if (element<1||element>length())
     {
       std::cout << "In ParVec: element not in my domain!"<<std::endl;
       std::cout <<"element = "<<element<<std::endl;
       std::cout<< "Mystart = "<<mystart<<", Myend = "<<myend<<std::endl;
       exit(0);
     }
 }
 double VecSum()
 {
   return IMPACT_Vector::VecSum();
 }
};

IMPACT_ParVec Gather(IMPACT_ParVec *P,IMPACT_Config *c,IMPACT_MPI_Config *MPIc)
{
  /*
    This function gathers a series of ParVecs and combines
    them into one long ParVec, which it returns.
   */
  IMPACT_ParVec answer(1,c->totalpoints());
  //NOw gather from all nodes:
  MPI_Barrier(MPI_COMM_WORLD);
  if (MPIc->rank()>0)
    {
      //needs delete!
      double * data;
      data = new double[MPIc->N()];
      for (int i=0;i<MPIc->N();++i)
	data[i]=P->Get(i+MPIc->start());
      //MPI::COMM_WORLD.Send(&data,MPIc->N(),MPI::DOUBLE,0,200+MPIc->rank());
      MPI_Send(&data,MPIc->N(),MPI_DOUBLE,0,200+MPIc->rank(),MPI_COMM_WORLD);
      delete[] data;
     }
  MPI_Barrier(MPI_COMM_WORLD);
  if (!MPIc->rank())
    {
      double *datarec;
      datarec = new double[c->totalpoints()];
      for (int i=0;i<MPIc->N();++i)
	datarec[i]=P->Get(i+1); //The start of the 0 rank processor at 1.
      
      for (int ranktemp=1;ranktemp<MPIc->size();ranktemp++)
	{
	  //MPI::COMM_WORLD.Recv(&datarec[MPIc->start(ranktemp)-1],MPIc->end(ranktemp)-MPIc->start(ranktemp)+1,MPI::DOUBLE,MPI::ANY_SOURCE,200+ranktemp);
    MPI_Recv(&datarec[MPIc->start(ranktemp)-1],MPIc->end(ranktemp)-MPIc->start(ranktemp)+1,MPI_DOUBLE,MPI_ANY_SOURCE,200+ranktemp,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
      
      for (int i=0;i<c->totalpoints();++i)
	answer.Set(i+1,datarec[i]);
      delete[] datarec;
    }
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  return answer;
}

void Gather_kstring(IMPACT_ParVec *P,IMPACT_Config *c,IMPACT_MPI_Config *MPIc,int *i, int *j,double *data)
{
  /*
    This function gathers a series of ParVecs and combines
    them into a single k-string, which it returns.
  */ 
  
 
  int klength=c->cellpoints();
  for (int k=0;k<klength;++k)
    data[k]=0.0;

  int istart=MPIc->istart();
  int iend=MPIc->iend();
  
  //NOw gather from all nodes:
  //MPI_Barrier(MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (*i>=istart&&*i<=iend)
    {
      int startofkstring=(*j-1+c->Ny()*(*i-1))*klength+1;
      if (MPIc->rank()>0)
	{

	  for (int ii=0;ii<klength;++ii)
	    data[ii]=P->Get(ii+startofkstring);
	   //MPI::COMM_WORLD.Send(&data[0],klength,MPI::DOUBLE,0,700+*i*100000+*j);
    MPI_Send(&data[0],klength,MPI_DOUBLE,0,700+*i*100000+*j,MPI_COMM_WORLD);

	}
      if (!MPIc->rank())
	{
	  for (int ii=0;ii<klength;++ii)
	    data[ii]=P->Get(startofkstring+ii);
	}
    }
  else if (!MPIc->rank())
    {
      //MPI::COMM_WORLD.Recv(&data[0],klength,MPI::DOUBLE,MPI::ANY_SOURCE,700+*i*100000+*j);
      MPI_Recv(&data[0],klength,MPI_DOUBLE,MPI_ANY_SOURCE,700+*i*100000+*j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }  
  MPI_Barrier(MPI_COMM_WORLD);
}
