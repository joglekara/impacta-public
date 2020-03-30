/*
**********************************************************
Variable class in IMPACTA code

Version 1.6
1
AGRT

12/2/07

5/3/07 - IMPORTANT CHANGE -IMPACT_Vector replaced with IMPACT_ParVec
-> parallel vector class.
and skip index improved slightly for vec objects,

NOW DIMENSION CLASS for clarity for the overloading
note - f2 and 3 need updating too.

18/4/07 - velocity differential stencil insert added

19/4/07 - get modified for f0 as it is always present - make code 
slightly quicker for Cee0 term

2/5/07 - extract_k added to speed up Cee0 term.

4/6/07 - Code for f2 put in - there are two forms, one
   has the format x1,x2 to identify the component -e.g xy, zz
the other is a single number from one to Nf2

27/7/07 - f2 code changed to be slightly more efficient (set method)

11/7/08 - Added method to rescale normalized variable depending 
          on new temperature normalization
**********************************************************

*/

/* This class is the cell of the simulation
each cell contains e fields, bfields and distrib
ution functions. These I guess can be added to

In this the construction and referencing for the ParVec containing
all the components of the distribution functions and field is
dealt with. One makes a variable class by 
e.g. IMPACT_Var E( nx,ny,nv,type)
nx - number of x points in grid
ny - ditto y, nv ditto v
type is a character defining whether the reference is 
E - efield
B - bfield
0,1,2,3 - f0,f1,f2 and f3
Then to access a particular cell at i(x),j(y) and k(v) in the 
ParVec one uses:
f0.get(i,j,k)
E.get(i,j,k,x1) where x1 is either x or y (similar for f1)
f2.get(i,j,k,x1,x2) etc.

x1,x2,x3 are integers where 1 is x, 2 is y and 3 is z.

19/1/07

This has now been made better - a config class with various switches to 
turn off componnents of the ParVec has been added that also takes NX,Ny 
and Nv arguments IMPACT_Config config(Nx,Ny,Nv,f0on,f1on,f2on,f3on,Eon,Bon)

This can also be defined IMPACT_Config config(Nx,Ny,Nv) which means that 
all variables are "on".
This means a variable can be defined simply like IMPACT_Var E(config,type)

22/1/07

Differential methods need a namespace declaration which gives the type of 
differential, whether upwind, cell centred etc...

1/2/07 the above (differential methods) have been removed as they were nonsense
That means dx,dy,dv probably can be removed at some point as a relic.
*/


class IMPACT_Var
{
 private:
  int Nv,Nx,Ny;  //Number of v,x,y points
  int objecttype; //defines whether f or field
  int totalincell; //total points in a single cell at i,j
  int N; //total length of ParVec
  /*This labels whether the object is  
  a field component or distribution function etc.
  The object type codes will be:
  0 - f0
  1 - f1
  2 - f2
  3 - f3
  E,4 - E
  B,5 - B
  pqr
  000 - f0
  100 - f1x
  010 - f1y
  001 - f1z
  200 - f2xx
  110 - f2xy
  101 - f2xz
  011 - f2yz
  020 - f2yy
  002 - f2zz (=-f2xx-f2yy)
  300 - f3xxx
  210 - f3xxy
  201 - f3xxz
  120 - f3xyy 
  102 - f3xzz (f3xzz = ?f3xxx ? f3xyy )
  111 - f3xyz
  030 - f3yyy (f3yyy = ?f3xxy ? f3yzz) 
  021 - f3yyz 
  012 - f3yzz
  003 - f3zzz (f3zzz = -f3xxz-f3yyz)
  400 - Ex
  401 - Ey
  500 - Bx
  501 - By
  502 - Bz*/
  int N_f0,N_f1,N_f2,N_f3,NE,NB; // These are numbers of 
        //components in the f's and B

  int ** vecskipindex; //in an attempt to speed up calculation.
  //vecskipindex [1][4][5]have 3 components with the appropriate skip from
  //earlier code
public:                                
  IMPACT_Var()
    {
      totalincell = (N_f0+N_f1+N_f2+N_f3)*Nv+NE+NB;
      Nv=1;Nx=1;Ny=1;objecttype=0; //really quite meaningless// so::
      std::cout<<"No cells specified!"<<std::endl;
      exit(0);
    }
~IMPACT_Var()
    {
      delete[] vecskipindex;
    }
  IMPACT_Var(int Num_x,int Num_y, int Num_v) // constructor for f0
    {
      N_f0=1;N_f1=3;
      N_f2=5;N_f3=7;
      NE=2;NB=3;
      objecttype=0; // create f0
      Nv=Num_v;
      Nx=Num_x;
      Ny=Num_y;
      totalincell = (N_f0+N_f1+N_f2+N_f3)*Nv+NE+NB;
      N=totalincell*Nx*Ny;
      
    }
  IMPACT_Var(IMPACT_Config config1) // constructor for f0
    { // from configuration object
      objecttype=0; // create f0
      N_f0=config1.conf_N_f0;N_f1=config1.conf_N_f1;
      N_f2=config1.conf_N_f2;N_f3=config1.conf_N_f3;
      NE=config1.conf_NE;NB=config1.conf_NB;
      Nv=config1.conf_Nv;
      Nx=config1.conf_Nx;
      Ny=config1.conf_Ny;
      totalincell = (N_f0+N_f1+N_f2+N_f3)*Nv+NE+NB;
      N=totalincell*Nx*Ny;
     
    }
  IMPACT_Var(int Num_x,int Num_y,char object_char) // constructor for E and B
    {
      switch (toupper(object_char))
	{
	case 'E':
	  objecttype=4;
	  break;
	case 'B':
	  objecttype=5;
	  break;
	default:
	  std::cout<<"No v cells specified!"<<std::endl;
	  exit(0);
	}
      N_f0=1;N_f1=3;
      N_f2=5;N_f3=7;
      NE=2;NB=3;
      Nv=1;
      Nx=Num_x;
      Ny=Num_y;
      
      totalincell = (N_f0+N_f1+N_f2+N_f3)*Nv+NE+NB;
      N=totalincell*Nx*Ny;
      vecskipindex=new int*[6];
      vecskipindex[objecttype]=new int[4];
      for (int x1=1;x1<4;++x1)
	{
	  switch(objecttype)
	    {
	    case 4:
	      vecskipindex[4][x1]=Nv*(N_f0+N_f1+N_f2+N_f3)+x1; //skip f0,1,2,3
	      // x1=either 1(x) or 2(y) for E 
	      break;
	    case 5:
	      vecskipindex[5][x1]=Nv*(N_f0+N_f1+N_f2+N_f3)+NE+(NB-3+x1);
	      break;
	    default:
	      break;
	    }
	}
    }
  IMPACT_Var(int Num_x,int Num_y, int Num_v,char object_char) // constructor for f0,1 2 nd 3
    {
       switch (toupper(object_char))
	{
	case 'E':
	  objecttype=4;
	  break;
	case 'B':
	  objecttype=5;
	  break;
	  case '0':
	  objecttype=0;
	  break;
	case '1':
	  objecttype=1;
	  break;
	case '2':
	  objecttype=2;
	  break;
	case '3':
	  objecttype=3;
	  break;
	default:
	   std::cout<<"IMPACT_Var (char,int) must be 'E','B' '1','2' or '3'"<<std::endl;
	}
       N_f0=1;N_f1=3;
       N_f2=5;N_f3=7;
       NE=2;NB=3;
       Nv=Num_v;
       Nx=Num_x;
       Ny=Num_y;
    
      totalincell = (N_f0+N_f1+N_f2+N_f3)*Nv+NE+NB;
      N=totalincell*Nx*Ny;
      setskip();
    }
IMPACT_Var(IMPACT_Config config1, char object_char)
 // constructor for f0,1 2 and 3 and E and B from config object
    {
       switch (toupper(object_char))
	{
	case 'E':
	  objecttype=4;
	  break;
	case 'B':
	  objecttype=5;
	  break;
	  case '0':
	  objecttype=0;
	  break;
	case '1':
	  objecttype=1;
	  break;
	case '2':
	  objecttype=2;
	  break;
	case '3':
	  objecttype=3;
	  break;
	default:
	  std::cout<<"IMPACT_Var (char,int) must be 'E','B' '1','2' or '3'"<<std::endl;
	}
       N_f0=config1.conf_N_f0;N_f1=config1.conf_N_f1;
       N_f2=config1.conf_N_f2;N_f3=config1.conf_N_f3;
       NE=config1.conf_NE;NB=config1.conf_NB;
       Nv=config1.conf_Nv;
       Nx=config1.conf_Nx;
       Ny=config1.conf_Ny;
       totalincell = (N_f0+N_f1+N_f2+N_f3)*Nv+NE+NB;
       N=totalincell*Nx*Ny;
       setskip();
    }
 void setskip()
 {
if (objecttype==1||objecttype==4||objecttype==5)
	{vecskipindex=new int*[6];
	  vecskipindex[objecttype]=new int[4];
      for (int x1=1;x1<4;++x1)
	{
	  switch(objecttype)
	    {
	    case 1:
	      vecskipindex[1][x1]=Nv*(N_f0+x1-1); //skips f0
	      break;
	    case 4:
	      vecskipindex[4][x1]=Nv*(N_f0+N_f1+N_f2+N_f3)+x1; //skip f0,1,2,3
	      // x1=either 1(x) or 2(y) for E 
	      break;
	    case 5:
	      vecskipindex[5][x1]=Nv*(N_f0+N_f1+N_f2+N_f3)+NE+(NB-3+x1);
	      break;
	    default:
	      break;
	    }
	}}
 }
   //________________________________________________________
 // New line of code for rescaling normalized variable to new temp
 void scale_T(IMPACT_Config *c, IMPACT_MPI_Config *MPIc, IMPACT_ParVec *v, double newT)
 {
  int istart=MPIc->istart();
  int iend=MPIc->iend();
  double lambda_s=1.0/(newT*newT); // scaling for distance
  double v_s=1.0/sqrt(newT); // scaling for velocity
  double t_s=1.0/sqrt(newT*newT*newT); // scaling for distance
  double scalar=1.0;
  // First work out scaling factor
  switch(objecttype)
    {
    case 0:
      //Normalization for f's
      scalar = 1.0/(v_s*v_s*v_s); break;
    case 1:
      scalar = 1.0/(v_s*v_s*v_s); break;
    case 2:
      scalar = 1.0/(v_s*v_s*v_s); break;
    case 3:
      scalar = 1.0/(v_s*v_s*v_s); break;
    case 4:
      scalar = lambda_s/(t_s*t_s); break;
    case 5:
      scalar = 1.0/t_s; break;
    default:
      break;
    }
  MPI::COMM_WORLD.Barrier();
  // Now multiply values
  for (int i=istart;i<=iend;++i) 
    for (int j=1;j<=c->Ny();++j)
      switch(objecttype)
	{
	case 0:
	  for (int k=1;k<=c->Nv();++k) set(v,get(v,&i,&j,&k)*scalar,&i,&j,&k);
	  break;
	case 1:
	  if (c->Nf1()>0)
	  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
	    for (int k=1;k<=c->Nv();++k) 
	      set(v,get(v,&i,&j,&k,&x1)*scalar,&i,&j,&k,&x1);
	  break;
	case 2:
	  if (c->Nf2()>0)
	  for (IMPACT_Dim x1=1;x1<=c->Nf2();++x1)
	    for (IMPACT_Dim x2=x1.get();x2<=c->Nf2();++x2)
	      for (int k=1;k<=c->Nv();++k) 
	      set(v,get(v,&i,&j,&k,&x1,&x2)*scalar,&i,&j,&k,&x1,&x2);
	  break;
	case 4:
	  if (c->NE()>0)
	  for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
	    set(v,get(v,&i,&j,&x1)*scalar,&i,&j,&x1);
	  break;
	case 5:
	  if (c->NB()>0)
	  for (IMPACT_Dim x1=3;x1>3-c->NB();--x1)
	    set(v,get(v,&i,&j,&x1)*scalar,&i,&j,&x1);
	  break;
	default:
	std::cout<<"IMPACTA: ERROR - in scale_T\n";
	exit(0);
	}
 }
 inline double get(IMPACT_ParVec *theParVec,int *i,int *j,int *k) // access for f0
 {
   // For wrong syntax:
   //________________________________________________________
   /*if (objecttype!=0)
     {std::cout<<"exiting - Not f0 being referenced"<<std::endl;exit(0);}*/
   //________________________________________________________
   
   int index;
   /* double answer;
      switch (N_f0)
      {
      case 1:
      index = *k+totalincell*(*j-1)+totalincell*Ny*(*i-1);
      // in other words after i-1 total blocks 
      // then j-1 sub Ny blocks then k in the cell
      answer= theParVec->Get(index);
      break;
      case 0:
	answer= 0.0;
	break;
	}
	return answer; */
    index = *k+totalincell*(*j-1+Ny*(*i-1));
    return theParVec->Get(index);
  }

   //________________________________________________________

inline  double get(IMPACT_ParVec *theParVec,int *i,int *j,int *k, IMPACT_Dim *x1) // access for f1 ONLY
 {
   /* if (objecttype!=1) {std::cout << "ERROR: not f1"<<std::endl;exit(0);}
   
   if (x1->get()>N_f1){std::cout<<"ERROR - f1.get - trying to access x1>N_f1"<<std::endl;exit(0);}*/
   //int skipindex=0; //ok - this skips over the f components to the appropriate
                    //component
   /*switch(objecttype)
     {
     case 1:
       skipindex=k+Nv*(N_f0+x1-1); //skips f0
       break;
     case 4:
       skipindex=Nv*(N_f0+N_f1+N_f2+N_f3)+x1; //skip f0,1,2,3
       // x1=either 1(x) or 2(y) for E 
       break;
     case 5:
       skipindex=Nv*(N_f0+N_f1+N_f2+N_f3)+NE+(NB-3+x1); //skip f0,1,2,3 and E
       break;
     default:
  // For wrong syntax:
     std::cout<<" Not f1, E or B being referenced"<<std::endl;
     exit(0);
     break;
     }*/
   int index;
   double answer=0.0;
   index = totalincell*((*j-1)+Ny*(*i-1))+vecskipindex[1][x1->get()]+*k;
   answer=theParVec->Get(index);
   /* in other words after i-1 total blocks 
   then j-1 sub Ny blocks then k in the cell
    */
   /*switch (NB)
     {
     case 3:
       if (objecttype==5&&x1+NB<=3) {answer=0.0; break;}
       if (objecttype==4&&x1>NE) {answer=0.0; break;}
       if (objecttype==1&&x1>N_f1) {answer=0.0; break;}
       answer=theParVec->Get(index);
       break;
     case 1:
       if (objecttype==5&&x1+NB<=3) {answer=0.0; break;}
       if (objecttype==4&&x1>NE) {answer=0.0; break;}
       if (objecttype==1&&x1>N_f1) {answer=0.0; break;}
       // if (x1==3) {answer=0.0; break;}
       answer=theParVec->Get(index);
       break;
     case 0:
       if (objecttype==5&&x1+NB<=3) {answer=0.0; break;}
       if (objecttype==4&&x1>NE) {answer=0.0; break;}
       if (objecttype==1&&x1>N_f1) {answer=0.0; break;}
       // if (x1==3) {answer=0.0; break;}
       answer=theParVec->Get(index);
       break;
       }*/
   return answer;
  }
inline double get(IMPACT_ParVec *theParVec,int *i,int *j, IMPACT_Dim *x1) // access for E&B ONLY
 {
   //Error checking - can be removed eventually for efficiency
   /*switch(objecttype) 
     {
     case 4:
       if (x1->get()>NE){std::cout<<"IMPACT ERROR - E.get() - trying to access x1>N_E"<<std::endl;exit(0);}
       break;
     case 5:
       if (x1->get()<=3-NB){std::cout<<"IMPACT ERROR - B.get() - trying to access x1>N_B"<<std::endl;exit(0);}
       break;
     default:
       std::cout << "ERROR in Var.get: not E or B->OT="<<objecttype<<std::endl;exit(0);
       }
   */
   double answer=0.0;
   int index;
   index = totalincell*((*j-1)+Ny*(*i-1))+vecskipindex[objecttype][x1->get()];
   answer=theParVec->Get(index);
   return answer;
 }
   //________________________________________________________

 double get(IMPACT_ParVec *theParVec,int *i,int *j,int *k,IMPACT_Dim *x1, IMPACT_Dim *x2) // access for f2
  {
    // For wrong syntax:
    //________________________________________________________
    /*if (objecttype!=2)
      {std::cout<<"Exiting - Not f2 being referenced"<<std::endl;exit(0);}*/
    //________________________________________________________
    //--(*x1); //so 1,2,3 beome 0,1,2
    //--(*x2);
    int x1i=x1->get()-1;//so 1,2,3 beome 0,1,2
    int x2i=x2->get()-1;
    if (x1i==2) x1i=3;//the reason for this is simply to allow the use of the
    if (x2i==2) x2i=3;//code x=0,y=1,z=3 so that x1+x2=0,1,2,3,4 -i.e. the
    //                components of 2l+1
    double answer;
    int index, skipindex;
    /* Slightly complicated, have to translate x1 and x2 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xx=0,xy=1,yy=2,xz=3,yz=4.... zz=-xx-yy
     */
    switch(x1i+x2i)
      {
      case 6:
	if (N_f2==0){answer=0.0;break;} //may be able to remove
	skipindex=Nv*(N_f0+N_f1); 
	index=*k+totalincell*((*j-1)+Ny*(*i-1))+skipindex;
	answer= -(theParVec->Get(index))-(theParVec->Get(index+2*Nv));
	//2*Nv is yy which is the 2nd f2
	break;
      default:
	if (x1i+x2i>(N_f2-1)){answer=0.0;break;}
	skipindex=Nv*(N_f0+N_f1+x1i+x2i);
	index = *k+totalincell*((*j-1)+Ny*(*i-1))+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	answer=theParVec->Get(index);
	break;      
      }
    // if (objecttype==4&&x1>1)answer=0.0;
    return answer;
    
  }
 /*
  //Ok, shit - now I have to do f3 - argh!
  double get(IMPACT_ParVec *theParVec,int i,int j,int k,IMPACT_Dim x1, IMPACT_Dim x2,IMPACT_Dim x3) 
// access for f3
  {
    // For wrong syntax:
    //________________________________________________________
    if (objecttype!=3)
      {std::cout<<"Exiting - Not f3 being referenced"<<std::endl;
	exit(0);}
    //________________________________________________________
   
    double answer;
    int index, skipindex;
 */
 /* Slightly complicated, have to translate x1 and x2 and x3 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xxx=0,xxy=1,xyy=2,yyy=3,xxz=4,xyz=5,zyy=6.... 
       xzz=-xxx-xyy
       yzz=-yxx-yyy
       zzz=-zyy-zxx
 */
 /*
    --x1; //so 1,2,3 beome 0,1,2
    --x2;
    --x3;
    if (x1==2) x1=4;//the reason for this is simply to allow the use of the
    if (x2==2) x2=4;//code x=0,y=1,z=4 so that x1+x2+x3=0,1,2,3,4 -i.e. the
    if (x3==2) x3=4;// components of 2l+1
    switch((x1+x2+x3).get())
      {
      case 8:
	if (N_f3==0){answer=0.0;break;}
	skipindex=Nv*(N_f0+N_f1+N_f2); 
	index=k+totalincell*(j-1)+totalincell*Ny*(i-1)+skipindex;
	answer= -(theParVec->Get(index))-(theParVec->Get(index+2*Nv));
	//2*Nv is yy which is the 2nd f2
	break;
      case 9:
	if (N_f3==0){answer=0.0;break;}
	skipindex=Nv*(N_f0+N_f1+N_f2); 
	index=k+totalincell*(j-1)+totalincell*Ny*(i-1)+skipindex;
	answer= -(theParVec->Get(index+Nv))-(theParVec->Get(index+3*Nv));
	//2*Nv is yy which is the 2nd f2
	break;
      case 12:
	if (N_f3==0){answer=0.0;break;}
	if (NB<3){answer=0.0;break;}
	skipindex=Nv*(N_f0+N_f1+N_f2); 
	index=k+totalincell*(j-1)+totalincell*Ny*(i-1)+skipindex;
	answer= -(theParVec->Get(index+4*Nv))-(theParVec->Get(index+6*Nv));
	//2*Nv is yy which is the 2nd f2
	break;
      default:
	if (x1+x2+x3>(N_f3-1)){answer=0.0;break;}
	skipindex=Nv*(N_f0+N_f1+N_f2)+(x1+x2+x3)*Nv;
	index = k+totalincell*(j-1)+totalincell*Ny*(i-1)+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	answer=theParVec->Get(index);
	break;      
      }
    if (objecttype==4&&x1>1)answer=0.0;
    return answer;
    
    }*/
 inline  void set(IMPACT_ParVec *theParVec,double value, int *i,int *j,int *k) // set method for f0
  {
    // For wrong syntax:
    //________________________________________________________
    /* if (objecttype!=0)
       {std::cout<<"exiting - Not f0 being referenced"<<std::endl;exit(0);}*/
    //________________________________________________________
    
    int index;
    /* switch (N_f0)
      {
      case 1:*/
    index = *k+totalincell*((*j-1)+Ny*(*i-1));
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	theParVec->Set(index,value);
	/*	break;
      case 0:
	break;
	}*/
  }
 inline void set(IMPACT_ParVec *theParVec, double value, int *i,int *j,int *k, IMPACT_Dim *x1)
 // set method for f1
 {
   /* if (objecttype!=1) {std::cout << "ERROR: not f1"<<std::endl;exit(0);}
      if (x1->get()>N_f1){std::cout<<"ERROR - f1.set - trying to access x1>N_f1"<<std::endl;exit(0);}*/

   int index;
   index = totalincell*((*j-1)+Ny*(*i-1))+vecskipindex[1][x1->get()]+*k;
   theParVec->Set(index,value);
  }
inline void set(IMPACT_ParVec *theParVec,double value,int *i,int *j, IMPACT_Dim *x1) // access for E&B ONLY
 {
   //Error checking - can be removed eventually for efficiency
   /*  switch(objecttype) 
     {
     case 4:
       if (x1->get()>NE){std::cout<<"ERROR - E.set() - trying to access x1>N_E"<<std::endl;exit(0);}
       break;
     case 5:
       if (x1->get()<=3-NB){std::cout<<"ERROR - B.set() - trying to access x1>N_B"<<std::endl;exit(0);}
       break;
     default:
       std::cout << "ERROR in Var.set: not E or B->OT="<<objecttype<<std::endl;exit(0);
       }*/
   int index;
   index = totalincell*((*j-1)+Ny*(*i-1))+vecskipindex[objecttype][x1->get()];
   theParVec->Set(index,value);
 }
   //________________________________________________________



 inline void set(IMPACT_ParVec *theParVec,double value,int *i,int *j,int *k,IMPACT_Dim *x1, IMPACT_Dim *x2)
 // set method for f2
  {
    // For wrong syntax:
    //________________________________________________________
    /* if (objecttype!=2)
       {std::cout<<"Exiting - Not f2 being referenced"<<std::endl;exit(0);}*/
    //________________________________________________________
    int x1i=x1->get()-1; //so 1,2,3 beome 0,1,2
    int x2i=x2->get()-1;
    if (x1i==2) x1i=3;//the reason for this is simply to allow the use of the
    if (x2i==2) x2i=3;//code x=0,y=1,z=3 so that x1+x2=0,1,2,3,4 -i.e. the
    //                components of 2l+1
   
    int index, skipindex;
    /* Slightly complicated, have to translate x1 and x2 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xx=0,xy=1,yy=2,xz=3,yz=4.... zz=-xx-yy
*/
    switch(x1i+x2i)
      {
      case 6:
	break;
      default:
	if (x1i+x2i>(N_f2-1)){;break;} //could remove
	skipindex=Nv*(N_f0+N_f1)+(x1i+x2i)*Nv;
	index = *k+totalincell*((*j-1)+Ny*(*i-1))+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	theParVec->Set(index,value);
	break;      
      } 
  }
inline void set(IMPACT_ParVec *theParVec,double value,int *i,int *j,int *k,int *position)
 // set method for f2 -alternative form
  {
    // For wrong syntax:
    //________________________________________________________
    /* if (objecttype!=2)
       {std::cout<<"Exiting - Not f2 being referenced"<<std::endl;exit(0);}*/
    //________________________________________________________
    
   
    int index, skipindex;
    /* Slightly complicated, have to translate x1 and x2 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xx=0,xy=1,yy=2,xz=3,yz=4.... zz=-xx-yy
*/
    /*switch(*position)
      {
      case 6:
	break;
	default:*/
    //if (*position>(N_f2-1)){;break;} //could remove
    if (*position<N_f2){
	skipindex=Nv*(N_f0+N_f1)+(*position)*Nv;
	index = *k+totalincell*((*j-1)+Ny*(*i-1))+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	theParVec->Set(index,value);}
	//	break;      
	//} 
  }
/*
 //Ok, shit - now I have to do f3 - argh!
  void set(IMPACT_ParVec *theParVec,double value,int i,int j,int k,IMPACT_Dim x1, IMPACT_Dim x2,IMPACT_Dim x3) 
// access for f3
  {
    // For wrong syntax:
    //________________________________________________________
    if (objecttype!=3)
      {std::cout<<"Exiting - Not f3 being referenced"<<std::endl;
	exit(0);}
    //________________________________________________________
   
    int index, skipindex;*/
    /* Slightly complicated, have to translate x1 and x2 and x3 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xxx=0,xxy=1,xyy=2,yyy=3,xxz=4,xyz=5,zyy=6.... 
       xzz=-xxx-xyy
       yzz=-yxx-yyy
       zzz=-zyy-zxx
  *//*
    --x1; //so 1,2,3 beome 0,1,2
    --x2;
    --x3;
    if (x1==2) x1=4;//the reason for this is simply to allow the use of the
    if (x2==2) x2=4;//code x=0,y=1,z=4 so that x1+x2+x3=0,1,2,3,4 -i.e. the
    if (x3==2) x3=4;// components of 2l+1
    switch((x1+x2+x3).get())
      {
      case 8:
	break;
      case 9:
	break;
      case 12:
	break;
      default:
	if (x1+x2+x3>(N_f3-1)){break;}
	skipindex=Nv*(N_f0+N_f1+N_f2)+(x1+x2+x3)*Nv;
	index = k+totalincell*(j-1)+totalincell*Ny*(i-1)+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	theParVec->Set(index,value);
	break;      
      }
    
      }*/
 //inc - incremental set methods
inline void inc(IMPACT_ParVec *theParVec,double value,int *i,int *j,int *k) 
 {
   set(theParVec,get(theParVec,i,j,k)+value,i,j,k);
 }
inline void inc(IMPACT_ParVec *theParVec,double value,int *i,int *j,int *k,IMPACT_Dim *x1) 
 {
   set(theParVec,get(theParVec,i,j,k,x1)+value,i,j,k,x1);
 }
inline void inc(IMPACT_ParVec *theParVec,double value,int *i,int *j,IMPACT_Dim *x1) 
 {
   set(theParVec,get(theParVec,i,j,x1)+value,i,j,x1);
 }
inline void inc(IMPACT_ParVec *theParVec,double value,int *i,int *j,int *k,IMPACT_Dim *x1, IMPACT_Dim *x2) 
 {
   set(theParVec,get(theParVec,i,j,k,x1,x2)+value,i,j,k,x1,x2);
 }
/*
void inc(IMPACT_ParVec *theParVec,double value,int i,int j,int k,IMPACT_Dim x1, IMPACT_Dim x2,IMPACT_Dim x3) 
 {
   set(theParVec,get(theParVec,i,j,k,x1,x2,x3)+value,i,j,k,x1,x2,x3);
   }*/
 //_____________________________________________________________________________

 //Now follows methods for getting the row number in the ParVec

 //_____________________________________________________________________________

inline int getrow(int *i,int *j,int *k) // access for f0
  {
     // For wrong syntax:
    //________________________________________________________
    /*  if (objecttype!=0)
	{std::cout<<"exiting - Not f0 being referenced"<<std::endl;exit(0);}*/
      //________________________________________________________
    
    int index;
    /* switch (N_f0)
      {
      case 1:*/
    index = *k+totalincell*((*j-1)+Ny*(*i-1));
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	/*	break;
      case 0:
	index= 0;//returns 0 if does not exist
	break;
	}*/
	return index; 
    //return  *k+totalincell*(*j-1)+totalincell*Ny*(*i-1);
  }

   //________________________________________________________

inline int getrow(int *i,int *j,int *k, IMPACT_Dim *x1) // access for f1 ONLY
 {
   /* if (objecttype!=1) {std::cout << "ERROR: not f1"<<std::endl;exit(0);}
   if (x1->get()>N_f1){std::cout<<"ERROR - f1.getrow - trying to access x1>N_f1"<<std::endl;exit(0);}
   */
   int index;
   index = totalincell*((*j-1)+Ny*(*i-1))+vecskipindex[1][x1->get()]+*k;
   return index;
   //return totalincell*(*j-1)+totalincell*Ny*(*i-1)+vecskipindex[1][x1->get()]+*k;
  }

 inline int getrow(int* i,int* j, IMPACT_Dim* x1) // access for E&B ONLY
 {
   //Error checking - can be removed eventually for efficiency
   /* switch(objecttype) 
     {
     case 4:
       if (x1->get()>NE){std::cout<<"ERROR - E.getrow() - trying to access x1>N_E"<<std::endl;exit(0);}
       break;
     case 5:
       if (x1->get()<=3-NB){std::cout<<"ERROR - B.getrow() - trying to access x1>N_B"<<std::endl;exit(0);}
       break;
     default:
       std::cout << "ERROR in Var.getrow: not E or B ->OT="<<objecttype<<std::endl;exit(0);
       }*/

   int index;
   index = totalincell*((*j-1)+Ny*(*i-1))+vecskipindex[objecttype][x1->get()];
   return index;
   // return totalincell*(*j-1)+totalincell*Ny*(*i-1)+vecskipindex[objecttype][x1->get()];
 }
   //________________________________________________________

 int getrow(int *i,int *j,int *k,IMPACT_Dim *x1i, IMPACT_Dim *x2i) // access for f2
  {
    // For wrong syntax:
    //________________________________________________________
    /*if (objecttype!=2)
      {std::cout<<"Exiting - Not f2 being referenced"<<std::endl;exit(0);}*/
    //________________________________________________________
    int x1=x1i->get()-1; //so 1,2,3 beome 0,1,2
    int x2=x2i->get()-1;
    if (x1==2) x1=3;//the reason for this is simply to allow the use of the
    if (x2==2) x2=3;//code x=0,y=1,z=3 so that x1+x2=0,1,2,3,4 -i.e. the
    //                components of 2l+1
    int index, skipindex;
    /* Slightly complicated, have to translate x1 and x2 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xx=0,xy=1,yy=2,xz=3,yz=4.... zz=-xx-yy
 */
    switch((x1+x2))
      {
      case 6:
	index=0; 
	break;
      default:
	if (x1+x2>=(N_f2)){index=0;break;}
	skipindex=Nv*(N_f0+N_f1)+(x1+x2)*Nv;
	index = *k+totalincell*((*j-1)+Ny*(*i-1))+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	break;      
      }
    // if (objecttype==4&&x1>1) index=0;
    if (index==0) {
      std::cout<<"IMPACT: ERROR - getrow trying to access non-existent f2 element\n"; exit(0);
    }
    return index;
    
  }
 int getrow(int *i,int *j,int *k,int *position) // access for f2
  {
    // For wrong syntax:
    //________________________________________________________
    /*if (objecttype!=2)
      {std::cout<<"Exiting - Not f2 being referenced"<<std::endl;exit(0);}*/
    //________________________________________________________

    int index, skipindex;
    /* Slightly complicated, have to translate x1 and x2 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xx=0,xy=1,yy=2,xz=3,yz=4.... zz=-xx-yy
 */
    switch(*position)
      {
      case 6:
	index=0; 
	break;
      default:
	if (*position>(N_f2-1)){index=0;break;}
	skipindex=Nv*(N_f0+N_f1)+(*position)*Nv;
	index = *k+totalincell*((*j-1)+Ny*(*i-1))+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	break;      
      }
  if (index==0) {
      std::cout<<"IMPACT: ERROR - getrow trying to access non-existent f2 element\n"; exit(0);
    }
    // if (objecttype==4&&x1>1) index=0;
    return index;
    
  }
 /*
  //Ok, shit - now I have to do f3 - argh!
  int getrow(int i,int j,int k,IMPACT_Dim x1, IMPACT_Dim x2,IMPACT_Dim x3) 
// access for f3
  {
    // For wrong syntax:
    //________________________________________________________
    if (objecttype!=3)
      {std::cout<<"Exiting - Not f3 being referenced"<<std::endl;
	exit(0);}
    //________________________________________________________
   

    int index, skipindex;*/
    /* Slightly complicated, have to translate x1 and x2 and x3 to pqr
       going to use system that x1x2 is equivalent to x1+x2=reference. 
       Remember that matrix is symmetric and traceless hence:
       xxx=0,xxy=1,xyy=2,yyy=3,xxz=4,xyz=5,zyy=6.... 
       xzz=-xxx-xyy
       yzz=-yxx-yyy
       zzz=-zyy-zxx
   *//*
    --x1; //so 1,2,3 beome 0,1,2
    --x2;
    --x3;
    if (x1==2) x1=4;//the reason for this is simply to allow the use of the
    if (x2==2) x2=4;//code x=0,y=1,z=4 so that x1+x2+x3=0,1,2,3,4 -i.e. the
    if (x3==2) x3=4;// components of 2l+1
    switch((x1+x2+x3.get())
      {
      case 8:
	index=0; 
	break;
      case 9:
	index=0; 
	break;
      case 12:
	index=0; 
	break;
      default:
	if (x1+x2+x3>(N_f3-1)){index=0; break;}
	skipindex=Nv*(N_f0+N_f1+N_f2)+(x1+x2+x3)*Nv;
	index = k+totalincell*(j-1)+totalincell*Ny*(i-1)+skipindex;
	// in other words after i-1 total blocks 
	// then j-1 sub Ny blocks then k in the cell
	break;      
      }
    if (objecttype==4&&x1>1)index=0; 
    return index;
    }*/
 //To extract the variable as a function of k at a point i,j
 inline void Extract_k(double *dataout, IMPACT_ParVec *v,int *i,int *j)
 {
   for (int k=1;k<=Nv;++k)
     dataout[k-1] = get(v,i,j,&k);
 }
  //_______________________________________________________________________
  //These functions are defined outside the class to allow the definition of
  //namespaces which contain the appropriate boundary conditions
  // for example for periodic f0 in x and y there is a unique piece of code
  // in a namespace, which is defined.
  inline void insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k);
  inline void insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1);
  inline void insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1);
  inline void insert_yonly(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1);
  inline void insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1,IMPACT_Dim *x2);
  /*
    void insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int i, int j, int k,IMPACT_Dim x1,IMPACT_Dim x2,IMPACT_Dim x3);*/
  //so these next six are the adjustable bits of code specifically designed for the bounds.
  inline void insert_BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k);
  inline void insert_f1BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1);
inline void insert_EBC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1);
inline void insert_yonly_EBC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1);
inline void insert_BBC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1);
inline void insert_BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1,IMPACT_Dim *x2);

// void insert_BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int i, int j, int k,IMPACT_Dim x1,IMPACT_Dim x2,IMPACT_Dim x3);

  //Overloaded versions of INSERT for a IMPACT_vint object - i.e. integral over v space.
 inline void insert(IMPACT_vint *vint, IMPACT_ParVec *v, int *i, int *j);
  inline void insert(IMPACT_vint *vint, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1);
  //and for velocity differential
  inline void insert(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k);
  inline void insert_BC(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k);
  //differential operators wrt v
  inline double ddv(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k);
  inline double ddv_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k);
  inline double ddv_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k, IMPACT_Dim *x1);
  inline double ddv1overv(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k,IMPACT_Dim *x1);
  inline double ddv1overv_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k,IMPACT_Dim *x1);
  inline double ddv_v2(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k, IMPACT_Dim* x1);
  inline double ddv_v2_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int* k, IMPACT_Dim *x1);

  inline double ddv_v3(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k, IMPACT_Dim* x1, IMPACT_Dim* x2);
  inline double ddv_v3_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int* k, IMPACT_Dim *x1, IMPACT_Dim* x2);
  //****************************************************************************
  // EFFICIENT ZEROING OF MATRIX ELEMENTS
  inline void clean(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, int *k); // for f0
  inline void clean(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, int *k, IMPACT_Dim * x1); //for f1
  inline void clean(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, IMPACT_Dim* x1); // for E/B
  inline void clean(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j,int *k,  int *position); // for f2
  inline void clean(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j,int *k,  IMPACT_Dim *x1,IMPACT_Dim *x2); // for f2
  //
 inline void cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, int *k); // for f0
  inline void cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, int *k, IMPACT_Dim * x1); //for f1
 inline void cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, IMPACT_Dim* x1); // for E/B
 inline void cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j, int *k,int *position); //for f2
  inline void cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j,int *k,  IMPACT_Dim *x1,IMPACT_Dim *x2); // for f2
 //Curl operator for E or B
 inline double Curl(IMPACT_ParVec * v,IMPACT_Config *config1,IMPACT_StenOps * O, int *i, int *j,  IMPACT_Dim * x1);

//print method for diagnostic purposes.
  void print(IMPACT_ParVec * v);
};

//for diagnostic purposes:
void IMPACT_Var::print(IMPACT_ParVec * v)
{
  IMPACT_Dim x1(1),x2(2),x3(3);
  switch (objecttype)
    {
    case 0:
      
      std::cout<<std::endl<<"f0"<<std::endl;
      for(int k=1;k<=Nv;++k){
	//	std::cout<<"k="<<k<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&k);//<<", ";
	//std::cout<<std::endl;
      }
      std::cout<<std::endl;
      }
      break;
    case 1:
      
      std::cout<<std::endl<<"f1x"<<std::endl;
for(int k=1;k<=Nv;++k){
  std::cout<<"k="<<k<<std::endl;
      for (int j=1;j<=Ny;++j){
	for (int i=1;i<=Nx;++i)
	  std::cout<< get(v,&i,&j,&k,&x1)<<", ";
	std::cout<<std::endl;
      }
std::cout<<std::endl;
      }
/*for(int k=1;k<=Nv;++k){
	std::cout<<"k="<<k<<std::endl;
      std::cout<<std::endl<<"f1y"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&k,&x2)<<", ";
	std::cout<<std::endl;
      }
std::cout<<std::endl;
      }
for(int k=1;k<=Nv;++k){
	std::cout<<"k="<<k<<std::endl;
      std::cout<<std::endl<<"f1z"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&k,&x3)<<", ";
	std::cout<<std::endl;
      }
std::cout<<std::endl;
}*/
      break;
      
    case 4:
      
      std::cout<<std::endl<<"Ex"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&x1)<<", ";
	std::cout<<std::endl;
      }
      std::cout<<std::endl<<"Ey"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&x2)<<", ";
	std::cout<<std::endl;
      }
      std::cout<<std::endl<<"Ez"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&x3)<<", ";
	std::cout<<std::endl;
      }
      break;
    case 5:
      
      std::cout<<std::endl<<"Bx"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&x1)<<", ";
	std::cout<<std::endl;
      }
      std::cout<<std::endl<<"By"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&x2)<<", ";
	std::cout<<std::endl;
      }
      std::cout<<std::endl<<"Bz"<<std::endl;
      for (int i=1;i<=Nx;++i){
	for (int j=1;j<=Ny;++j)
	  std::cout<< get(v,&i,&j,&x3)<<", ";
	std::cout<<std::endl;
      }
      break;
    }
}

