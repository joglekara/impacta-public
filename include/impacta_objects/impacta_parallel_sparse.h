/*
**********************************************************
IMPACT Version 1.0
Parallel Sparse matrix class

Version 1.2
AGRT

13/2/07
This class creates a normal sparse matrix, but has an 
indexed start and end point for the sparse matrix
interms of rownumber. The rownumber vector is smaller
via a new constructor in the sparse class.

9/4/07 
New functions nullequation and count_nullequation added
which stick a 1 on the diagonal if the equation is not
wanted in the matrix
**********************************************************
*/
class IMPACT_ParSpa : public IMPACT_Sparse
{
 private:
  int rowstart,rowend; //rowvector goes from rowstart to rowend +1!!!!  
  int rowlength;
 public:
 IMPACT_ParSpa() : IMPACT_Sparse()
    {
      rowstart=1;
      rowend=1;
      rowlength=1;
    }

 IMPACT_ParSpa(IMPACT_MPI_Config *MPIc) : IMPACT_Sparse(MPIc->N())
   //i.e. the necessary rows - no ghosts
    {
      rowstart=MPIc->start();
      rowend=MPIc->end();
      rowlength=rowend-rowstart+1;
    }
 inline int length()
  {
    return rowlength;
  }
  //now need a series of access methods that translate the
  //inputted rows to the rows actaully stored.
 inline  int getrow(int index)//returns value in row index
  {
    int indshift=index-rowstart+1;
    //chkerror(indshift);
    return IMPACT_Sparse::getrow(indshift);
  }
inline  int getrow_nocheck(int index)//returns value in row index
  {
    int indshift=index-rowstart+1;
    return IMPACT_Sparse::getrow_nocheck(indshift);
  }
inline void setrow(int index,int value)//sends value to row index
{
  int indshift=index-rowstart+1;
  //chkerror(indshift);
  IMPACT_Sparse::setrow(indshift,value);
}
inline void incrow(int index,int value)//sends value to row index
{
  int indshift=index-rowstart+1;
  //chkerror(indshift);
  IMPACT_Sparse::incrow(indshift,value);
}
inline void Resize() // after count, resizes column and row vectors
  {
    IMPACT_Sparse::Resize();
  }
 inline void chkerror(int element)
 {
   /*if (element<1||element>length()+1)
     {
       std::cout << "In ParSpa: element not in my domain!"<<std::endl;
       std::cout <<"element = "<<element<<std::endl;
       std::cout<< "Rowstart = "<<rowstart<<", Rowend = "<<rowend<<std::endl;
       exit(0);
       }*/
   }
inline void nullequation(IMPACT_Config *c, int type,int rownumber)
 {
   int index = getrow(rownumber);
   setcol(index,rownumber);
   /*
I don't think there is any need for this now.
     int check_equation_switch=0;
   switch (type)
     {
     case 0:
       if (equation_switches::f0_equation_on) check_equation_switch=1; break;
     case 1:
       if (equation_switches::f1_equation_on) check_equation_switch=1; break;
     case 2:
       if (equation_switches::f2_equation_on) check_equation_switch=1; break;
     case 4:
       if (equation_switches::E_equation_on) check_equation_switch=1; break;
     case 5:
       if (equation_switches::B_equation_on) check_equation_switch=1; break;
     }
     if (check_equation_switch) setval(index,c->idt());
   if (!check_equation_switch) */
   setval(index,1.0);

 }
 inline void count_nullequation(IMPACT_Config *c, int rownumber)
 {
   setrow(rownumber+1,1);
 }
 inline void flipsign()
 {
   for (int i=1;i<=getSl();++i)
     setval(i,-1.0*getval(i));
 }
};
//For Sparse Mtrix checking
 inline void chksparse_chk(IMPACT_ParSpa *S, int *rownumber,int *runningtotal)
  {
    int a=*runningtotal;
    int b=S->getrow_nocheck(*rownumber+1);
    if (a>b) //numinrows[*rownumber])
    {
	std::cout<<"IMPACT: ERROR - insufficient memory reserved in Sparse row "<<*rownumber<<std::endl;
	std::cout<< "Tried to access Sparse vector element "
		 <<a<<std::endl<<"Start of next row is at - "
		 <<b<<'\n';
	exit(0);
      }
  }
namespace no_sparserow_checking
{
  int sparserowcheck=0;
  inline  void chksparse(IMPACT_ParSpa *S, int *rownumber,int *runningtotal)
  {}
}
namespace with_sparserow_checking
{
  int sparserowcheck=1;
 inline  void chksparse(IMPACT_ParSpa *S, int *rownumber,int *runningtotal)
 {chksparse_chk(S,rownumber,runningtotal);}
}