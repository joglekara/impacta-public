/*
**********************************************************
IMPACTA Configuration info class

1
AGRT

2/4/07
**********************************************************

10/4/08 - x,y positions updated to start at dx/2
*/
/*
In this the construction and referencing for the vector containing
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
vector one uses:
f0.get(i,j,k)
E.get(i,j,k,x1) where x1 is either x or y (similar for f1)
f2.get(i,j,k,x1,x2) etc.

x1,x2,x3 are integers where 1 is x, 2 is y and 3 is z.

19/1/07

This has now been made better - a config class with various switches to 
turn off componnents of the vector has been added that also takes NX,Ny 
and Nv arguments IMPACT_Config config(Nx,Ny,Nv,f0on,f1on,f2on,f3on,Eon,Bon)

This can also be defined IMPACT_Config config(Nx,Ny,Nv) which means that 
all variables are "on".
This means a variable can be defined simply like IMPACT_Var E(config,type)

22/1/07

Now added is dx,dy and dv to the config class

23/1/07

Now a Matrix class IMPACT_Matrix M(N)
M.get(i,j)
and M.set(i,j,value) - equivalent to vector member functions

24/2 - Matrix class moved to own header file

5/3/07 - New addition - local lower and upper bounds for i and j
         -> this allows the node to know what its local limits are.

2/4/07 - new additions include xmin and max and ndump and nmax
2/5/07 - v2 function altered so that it returns from a new array
- to speed up Cee0 term.

25/7/07 - new access method to get number of components of a particular f

11/7/08 - Added temperature scaling. This scales v by T^0.5, dx by T^2
and dv by T^1.5

*/

class IMPACT_Config //Configuration data variable
{
 private:
  friend class IMPACT_Var;
  friend class IMPACT_Sparse;
  friend class IMPACT_Matrix;
  int conf_Nx,conf_Ny,conf_Nv;
  int conf_N_f0,conf_N_f1,conf_N_f2,conf_N_f3,conf_NE,conf_NB;
  double *conf_v; // values for v
  double *conf_v2; // values for v squared
  double *conf_dx;
  double *conf_dy;
  double *conf_dv;
  double conf_dt;
  double *conf_dx_inv;//these are the inverted forms of the above
  double *conf_dy_inv;//to speed up calculations (* better than /)
  double *conf_dv_inv;
  double conf_dt_inv;  
  
  double xinit,yinit; //the initial values of x and y so the grid can be c
  //alculated
  int Ndump, NMAX; // iteration information
// remember these int N_f0=1,N_f1=3,N_f2=5,N_f3=7,E=2,NB=3;
 public:
IMPACT_Config()
    {
      conf_N_f0=1;conf_N_f1=3;
      conf_N_f2=5;conf_N_f3=7;
      conf_NE=3;conf_NB=3;
      conf_Nx=1;
      conf_Ny=1;
      conf_Nv=1;
      conf_dx=new double[conf_Nx+2];
      conf_dy=new double[conf_Ny+2];
      conf_dv=new double[conf_Nv+2];
      conf_dx_inv=new double[conf_Nx+2];
      conf_dy_inv=new double[conf_Ny+2];
      conf_dv_inv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx_inv[i]=conf_dx[i]=1.0;
      for (int i =0;i<conf_Ny+2;++i) conf_dy_inv[i]=conf_dy[i]=1.0;
      for (int i =0;i<conf_Nv+2;++i) conf_dv_inv[i]=conf_dv[i]=1.0;
      conf_dt=1.0;
      conf_dt_inv=1.0;
      conf_v=new double[conf_Nv+2];
      conf_v2=new double[conf_Nv+2];
      conf_v[0]=0;
      conf_v2[0]=0;
      for (int i =1;i<conf_Nv+2;++i)
	{
	conf_v[i]=conf_v[i-1]+0.5*(conf_dv[i-1]+conf_dv[i]);
	conf_v2[i]=conf_v[i]*conf_v[i];
	}
      xinit=yinit=0.0;
      Ndump =NMAX=1;
    }

  IMPACT_Config(int Numx,int Numy,int Numv)
    {
      conf_N_f0=1;conf_N_f1=3;
      conf_N_f2=5;conf_N_f3=7;
      conf_NE=3;conf_NB=3;
      conf_Nx=Numx;
      conf_Ny=Numy;
      conf_Nv=Numv;
      conf_dx=new double[conf_Nx+2];
      conf_dy=new double[conf_Ny+2];
      conf_dv=new double[conf_Nv+2];
      conf_dx_inv=new double[conf_Nx+2];
      conf_dy_inv=new double[conf_Ny+2];
      conf_dv_inv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx_inv[i]=conf_dx[i]=1.0;
      for (int i =0;i<conf_Ny+2;++i) conf_dy_inv[i]=conf_dy[i]=1.0;
      for (int i =0;i<conf_Nv+2;++i) conf_dv_inv[i]=conf_dv[i]=1.0;
      conf_dt=1.0;
      conf_dt_inv=1.0;
      conf_v=new double[conf_Nv+2];
      conf_v2=new double[conf_Nv+2];
      conf_v[0]=0;
      conf_v2[0]=0;
      for (int i =1;i<conf_Nv+2;++i)
	{
	conf_v[i]=conf_v[i-1]+0.5*(conf_dv[i-1]+conf_dv[i]);
	conf_v2[i]=conf_v[i]*conf_v[i];
	}
      xinit=yinit=0.0;
      Ndump =NMAX=1;
    }
  IMPACT_Config(int Numx,int Numy,int Numv, int f0,int f1,int f2, int f3, int E,int B) 
    /*Switches for the various components
    f0 -> E can be 1 or 0 (on or off)
    B can be 0,1 or 3 (i.e. 0, z, or xyz)
    The Numx etc. are the size of the grid in x,y, and v
    This will allow easy switching and information transfer
    */
    {
      try{
      if (f0>1||f0<0) throw 11;
      if (f1>1||f1<0) throw 11;
      if (f2>1||f2<0) throw 11;
      if (f3>1||f3<0) throw 11;
      if (E>1||E<0) throw 11;
      if (B>3||B<0) throw 12;
      if(B==2) throw 13;
      }
      catch (int err)
	{
	  if (err==11){std::cout<<"f or E Switch must be 1 or 0"<<std::endl;exit(0);}
	  if (err==12){std::cout<<"B components limited to 3!"<<std::endl;exit(0);}
	  if (err==13){std::cout<<"B components can only be 0,1 or 3"<<std::endl;exit(0);}
	}
      conf_N_f0=1*f0;
      switch (B)
	{
	case 3:
	  conf_N_f1=3*f1;
	  conf_N_f2=5*f2;
	  conf_N_f3=7*f3;
	  conf_NE=3*E;
	  break;
	case 1:
	  conf_N_f1=2*f1;
	  conf_N_f2=3*f2;
	  conf_N_f3=4*f3;
	  conf_NE=2*E;
	  break;
	case 0:
	  conf_N_f1=2*f1;
	  conf_N_f2=3*f2;
	  conf_N_f3=4*f3;
	  conf_NE=2*E;
	  break;
	}     
      conf_NB=B;
      conf_Nx=Numx;
      conf_Ny=Numy;
      conf_Nv=Numv;
      conf_dx=new double[conf_Nx+2];
      conf_dy=new double[conf_Ny+2];
      conf_dv=new double[conf_Nv+2];
      conf_dx_inv=new double[conf_Nx+2];
      conf_dy_inv=new double[conf_Ny+2];
      conf_dv_inv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx_inv[i]=conf_dx[i]=1.0;
      for (int i =0;i<conf_Ny+2;++i) conf_dy_inv[i]=conf_dy[i]=1.0;
      for (int i =0;i<conf_Nv+2;++i) conf_dv_inv[i]=conf_dv[i]=1.0;
      conf_dt=1.0;
      conf_dt_inv=1.0;
      conf_v=new double[conf_Nv+2];
      conf_v2=new double[conf_Nv+2];
      conf_v[0]=0;
      conf_v2[0]=0;
      for (int i =1;i<conf_Nv+2;++i)
	{
	conf_v[i]=conf_v[i-1]+0.5*(conf_dv[i-1]+conf_dv[i]);
	conf_v2[i]=conf_v[i]*conf_v[i];
	}
      xinit=yinit=0.0;
      Ndump =NMAX=1;
    }
  IMPACT_Config(int Numx,int Numy,int Numv,double *dxin,double *dyin, double *dvin,double dtin, int f0,int f1,int f2, int f3, int E,int B) 
    /*Switches for the various components
    f0 -> E can be 1 or 0 (on or off)
    B can be 0,1 or 3 (i.e. 0, z, or xyz)
    The Numx etc. are the size of the grid in x,y, and v
    This will allow easy switching and information transfer

    Then also dxin,dyin and dvin are the dx's
    */
    {
      try{
      if (f0>1||f0<0) throw 11;
      if (f1>1||f1<0) throw 11;
      if (f2>1||f2<0) throw 11;
      if (f3>1||f3<0) throw 11;
      if (E>1||E<0) throw 11;
      if (B>3||B<0) throw 12;
      if(B==2) throw 13;
      }
      catch (int err)
	{
	  if (err==11){std::cout<<"f or E Switch must be 1 or 0"<<std::endl;exit(0);}
	  if (err==12){std::cout<<"B components limited to 3!"<<std::endl;exit(0);}
	  if (err==13){std::cout<<"B components can only be 0,1 or 3"<<std::endl;exit(0);}
	}
      conf_N_f0=1*f0;
      switch (B)
	{
	case 3:
	  conf_N_f1=3*f1;
	  conf_N_f2=5*f2;
	  conf_N_f3=7*f3;
	  conf_NE=3*E;
	  break;
	case 1:
	  conf_N_f1=2*f1;
	  conf_N_f2=3*f2;
	  conf_N_f3=4*f3;
	  conf_NE=2*E;
	  break;
	case 0:
	  conf_N_f1=2*f1;
	  conf_N_f2=3*f2;
	  conf_N_f3=4*f3;
	  conf_NE=2*E;
	  break;
	}
      
      conf_NB=B;
      conf_Nx=Numx;
      conf_Ny=Numy;
      conf_Nv=Numv;
      conf_dx=new double[conf_Nx+2];
      conf_dy=new double[conf_Ny+2];
      conf_dv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx[i]=dxin[i];
      for (int i =0;i<conf_Ny+2;++i) conf_dy[i]=dyin[i];
      for (int i =0;i<conf_Nv+2;++i) conf_dv[i]=dvin[i];
      conf_dx_inv=new double[conf_Nx+2];
      conf_dy_inv=new double[conf_Ny+2];
      conf_dv_inv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx_inv[i]=1.0/conf_dx[i];
      for (int i =0;i<conf_Ny+2;++i) conf_dy_inv[i]=1.0/conf_dy[i];
      for (int i =0;i<conf_Nv+2;++i) conf_dv_inv[i]=1.0/conf_dv[i];
      conf_dt=dtin;
      conf_dt_inv=1.0/conf_dt;
      conf_v=new double[conf_Nv+2];
      conf_v2=new double[conf_Nv+2];
      conf_v[0]=0;
      conf_v2[0]=0;

      conf_v[0]=-0.5*conf_dv[0];
      conf_v2[0]=conf_v[0]*conf_v[0];

      for (int i =1;i<conf_Nv+2;++i)
	{
	conf_v[i]=conf_v[i-1]+0.5*(conf_dv[i-1]+conf_dv[i]);
	conf_v2[i]=conf_v[i]*conf_v[i];
	}
      xinit=yinit=0.0;
      Ndump =NMAX=1;
    }
  IMPACT_Config(int Numx,int Numy,int Numv,double *dxin,double *dyin, double *dvin,double dtin)
    {
      conf_N_f0=1;conf_N_f1=3;
      conf_N_f2=5;conf_N_f3=7;
      conf_NE=3;conf_NB=3;
      conf_Nx=Numx;
      conf_Ny=Numy;
      conf_Nv=Numv;
      conf_dx=new double[conf_Nx+2];
      conf_dy=new double[conf_Ny+2];
      conf_dv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx[i]=dxin[i];
      for (int i =0;i<conf_Ny+2;++i) conf_dy[i]=dyin[i];
      for (int i =0;i<conf_Nv+2;++i) conf_dv[i]=dvin[i];
      conf_dx_inv=new double[conf_Nx+2];
      conf_dy_inv=new double[conf_Ny+2];
      conf_dv_inv=new double[conf_Nv+2];
      for (int i =0;i<conf_Nx+2;++i) conf_dx_inv[i]=1.0/conf_dx[i];
      for (int i =0;i<conf_Ny+2;++i) conf_dy_inv[i]=1.0/conf_dy[i];
      for (int i =0;i<conf_Nv+2;++i) conf_dv_inv[i]=1.0/conf_dv[i];
      conf_dt=dtin;
      conf_dt_inv=1.0/conf_dt;
      conf_v=new double[conf_Nv+2];
      conf_v2=new double[conf_Nv+2];
      conf_v[0]=0;
      conf_v2[0]=0;
      for (int i =1;i<conf_Nv+2;++i)
	{
	conf_v[i]=conf_v[i-1]+0.5*(conf_dv[i-1]+conf_dv[i]);
	conf_v2[i]=conf_v[i]*conf_v[i];
	}
      xinit=yinit=0.0;
      Ndump =NMAX=1;
    }

  //*************************************
  //Access methods
  //*************************************

  int cellpoints() // the total number of points in a particular cell x1,y1
  {
    return (conf_Nv*(conf_N_f1+conf_N_f0+conf_N_f2+conf_N_f3)+conf_NE+conf_NB);
  }
  int totalpoints() // total number of points in all variables in x,y and v
  {
    return cellpoints()*conf_Nx*conf_Ny;
  }
  int Nx()
  {
    return conf_Nx;
  }
  int Ny()
  {
    return conf_Ny;
  }
  int Nv()
  {
    return conf_Nv;
  }
  double dx(int *index)
  {
    checkerr(index,conf_Nx);
    return conf_dx[*index];
  }
  double dy(int *index)
  {
    checkerr(index,conf_Ny);
    return conf_dy[*index];
  }
  double dv(int *index)
  {
    checkerr(index,conf_Nv);
    return conf_dv[*index];
  }
 double idx(int *index)
  {
     checkerr(index,conf_Nx);
    return conf_dx_inv[*index];
  }
  double idy(int *index)
  {
     checkerr(index,conf_Ny);
    return conf_dy_inv[*index];
  }
  double idv(int *index)
  {
     checkerr(index,conf_Nv);
    return conf_dv_inv[*index];
  }
double v(int *index)
  {
     checkerr(index,conf_Nv);
    return conf_v[*index];
  }
 double v2(int *index) // v squared
  {
    //checkerr(index,conf_Nv);
     
     // return conf_v[*index]*conf_v[*index];
     return conf_v2[*index];
  }
double v3(int *index) // v cubed
  {
     checkerr(index,conf_Nv);
    return conf_v2[*index]*conf_v[*index];
  }
 double dt()
 {
   return conf_dt;
 }
double idt()
 {
   return conf_dt_inv;
 }
  void checkerr(int *index, int Ncheck)
  {
    if(*index<0||*index>Ncheck+1){std::cout<<"Outside range - config.dx,y,v(int)\nelement = "<<*index<<'\n';exit(0);}
  }
  int Nf(int n)
  {
    int N_of_f=0;
    switch (n)
      {
      case 0:
	N_of_f=1; break;
      case 1:
	N_of_f=conf_N_f1; break;
      case 2:
	N_of_f=conf_N_f2; break;
      case 3:
	N_of_f=conf_N_f3; break;
      default:
	std::cout<<"Error in Nf in IMPACT_Config\n";
	exit(0);
      }
    return N_of_f;
  }
  int Nf1()
  {
    return conf_N_f1;
  }
  int Nf2()
  {
    return conf_N_f2;
  }
int N3f2()
  {
    return (conf_N_f2+1)/2;
  }
  int Nf3()
  {
    return conf_N_f3;
  }
  int NE()
  {
    return conf_NE;
  }
  int NB()
  {
    return conf_NB;
  }
  int n_max()
  {
    return NMAX;
  }
  int N_Dump()
  {
    return Ndump;
  }
  //Set functions rather than  including in the constructor
  void set_xmin(double xmin)
  {
    xinit=xmin;
  }
  void set_ymin(double ymin)
  {
    yinit=ymin;
  }
  double get_xmin()
  {
    return xinit-0.5*conf_dx[0];
  }
  double get_ymin()
  {
    return yinit-0.5*conf_dy[0];
  }
  double get_xmax()
  {
    double pos=xinit-0.5*conf_dx[0];
    
    for (int i=0;i<conf_Nx;++i)
       pos+=0.5*(conf_dx[i]+conf_dx[i+1]); //remeber that conf_dx[0] is out of the domain

    return pos;
  }
  double get_ymax()
  {
    double pos=yinit-0.5*conf_dy[0];
    
     for (int i=0;i<conf_Ny;++i)
       pos+=0.5*(conf_dy[i]+conf_dy[i+1]); //remeber that conf_dx[0] is out of the domain
    
    return pos;
  }
  double xb(int index)
  {
   double pos=xinit-0.5*conf_dx[0];
    
     for (int i=1;i<=index;++i)
       pos+=(conf_dx[i]);  //remeber that conf_dx[0] is out of the domain

     // [0] is the left boundary [Nx] is the right boundary,

    return pos;
  }
  double yb(int index)
  {
    double pos=yinit-0.5*conf_dy[0];
    
     for (int i=1;i<=index;++i)
       pos+=(conf_dy[i]); //remeber that conf_dx[0] is out of the domain
    // [0] is the left boundary [Ny] is the right boundary,
    return pos;
  }
  void set_nmax(int N)
  {
    NMAX=N;
  }
  void set_ndump(int N)
  {
    Ndump = N;
  }
  //to get x and y grids:
  double xpos(int index)
  {
    double pos=xinit-conf_dx[0];
    
     for (int i=0;i<index;++i)
       pos+=0.5*(conf_dx[i]+conf_dx[i+1]); //remeber that conf_dx[0] is out of the domain
    
    return pos;
  }
  double ypos(int index)
  {
    double pos=yinit-conf_dy[0];
    
     for (int i=0;i<index;++i)
       pos+=0.5*(conf_dy[i]+conf_dy[i+1]); //remeber that conf_dx[0] is out of the domain
    
    return pos;
  }

  //Need to be able to change dt....
  void set_dt(double newdt)
  {
    conf_dt=newdt;
    if (newdt) {
      conf_dt_inv=1.0/conf_dt;
    } 
    else {
      conf_dt_inv=0.0;
    }
  }
  //This rescales the problem with new normalizing temperature
  void scale_T(double newT)
  {
    double lambda_s=1.0/(newT*newT); // scaling for distance
    double v_s=1.0/sqrt(newT); // scaling for velocity
    double t_s=1.0/sqrt(newT*newT*newT); // scaling for distance
    
    for (int i =0;i<conf_Nx+2;++i) conf_dx[i]*=lambda_s;
    for (int i =0;i<conf_Ny+2;++i) conf_dy[i]*=lambda_s;
    for (int i =0;i<conf_Nv+2;++i) conf_dv[i]*=v_s;
 
    for (int i =0;i<conf_Nx+2;++i) conf_dx_inv[i]=1.0/conf_dx[i];
    for (int i =0;i<conf_Ny+2;++i) conf_dy_inv[i]=1.0/conf_dy[i];
    for (int i =0;i<conf_Nv+2;++i) conf_dv_inv[i]=1.0/conf_dv[i];
    
    conf_dt*=t_s;
    
    conf_dt_inv=1.0/conf_dt;
 
    for (int i =1;i<conf_Nv+2;++i)
      {
	conf_v[i]=conf_v[i-1]+0.5*(conf_dv[i-1]+conf_dv[i]);
	conf_v2[i]=conf_v[i]*conf_v[i];
      }
  }
};
   //_____________________________________________
