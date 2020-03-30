/*
**********************************************************
INITILIZATION

This is at present (22/2/07) very rudimentary code for testing.
Version 1.3
AGRT

22/2/07

24/3/07 - Now more involved - takes data from .imd file
4/4/07 - New class - IMatrix class which is an initial 
condition (i.e. ni, Z etc) when given a reference it either 
has a constant value or a matrix the size of the 2d grid
10/4/07 - multiplyall added to IMatrix which allows all
values to be multiplied by a constant
26/4/07 - added fixed boundary code
8/5 - stored RB coefficents here
9/8/07 - IMatrix modified to include new set method.

25/9/07 - An Ez initiator added to simulate grad n x grad T
 b field generation in the plane.
 basically acts like the heating operator 

29/1/08 - Now added capability to initialize with DLM distribution

9/7/08 - Added function for initializing Ci at walls
**********************************************************
*/

double Gamma_Func(double a);

namespace Initial_Conditions
{
  IMPACT_IMatrix Z(1.0);
  IMPACT_IMatrix ni(1.0);
  IMPACT_IMatrix ne(1.0);
  IMPACT_IMatrix Te(1.0); 
  IMPACT_IMatrix B[3];
  IMPACT_IMatrix C_i[3]; // New - ion flow velocity
  IMPACT_IMatrix C_istar[3]; // New - ion flow velocity for intermediate step
  IMPACT_IMatrix DivC_i;
  double T_mult = 1.0;
  double f0delta = 0.0;
}
namespace IMPACT_Heating
{
  // These fields are all concerned with IB laser absorption
  IMPACT_IMatrix Total_heating_contribution(0.0);
  IMPACT_IMatrix Heating_xy(0.0);
  IMPACT_Matrix i_mat;
  IMPACT_Vector Heating_t;

  double IB_type_on = 0.0;
  double MX_type_on = 0.0;
  double vosc_squared = 0.01;


  //bool tr_bool = 0;
  double beam_width = 0.0;
  int shape = 0;
  int beam_res = 0;
  double n_c = 1.0;
  int direction = 0;
  double theta = 0.0;
  double dtheta = 0.0;
  double ds = 0.0;
  int diffsteps = 0;


  double ray_x0 = 0.0;
  double ray_y0 = 0.0;

  char polarization = 'y'; //polarization of laser
  IMPACT_Matrix vosc_hat(3,3); // This matrix corresponds to the 
  // stress unit tensor of the laser field
  //----------------------------
  // these matrices are to do with the effect of static magnetic
  //fields on the IB absorption
  IMPACT_Matrix vosc_hat_Bx(3,3);
  IMPACT_Matrix vosc_hat_By(3,3);
  IMPACT_Matrix vosc_hat_Bz(3,3);
  //----------------------------
  IMPACT_IMatrix Dnz_xy(0.0);
  IMPACT_Vector Dnz_t;
  IMPACT_IMatrix DTz_xy(0.0);
  IMPACT_Vector DTz_t;
  IMPACT_IMatrix DTbydzgrid(0.0);
  IMPACT_IMatrix Dnbydzgrid(0.0);
}
namespace Cee0_flux_store
{
  IMPACT_VS_Array *flux;
  IMPACT_VS_Array *f2_IB;
}
void Initialize(IMPACT_ParVec * v, IMPACT_Config * c, IMPACT_MPI_Config * MPIc,
		IMPACT_Var * f0,IMPACT_Var * f1,IMPACT_Var * f2,
		IMPACT_Var * f3,IMPACT_Var * E,IMPACT_Var * B)
{
  /* //Input electron density profile
  struct IMPACT_Moment ne; 
  IMPACT_Read(&ne,"ne_init.imd");
  
  //input e Temperatur profile
  struct IMPACT_Moment Ue;
  IMPACT_Read(&Ue,"Ue_init.imd");
  */
  //Now fill in values
  double Local_Maxwellian,Local_kbTe; //To initialize f0...
  //for test 4.4
  /*double d=1.0,x0=1.0,x,y,xl,yl;
  Initial_Conditions::Te.SwitchOffConst();
  Initial_Conditions::Te.ChangeSize(c->Nx(),c->Ny());
  
  for (int i=1;i<=c->Nx();++i)
    for (int j=1;j<=c->Ny();++j)
      {
	xl=c->xpos(c->Nx());
	yl=c->ypos(c->Ny());
	x=c->xpos(i);
	y=c->ypos(j);
	d=xl/6.0;
	x0 = (xl/2) + (yl/8)*cos(pi*y/yl);
	Local_kbTe =1.0+ 0.005*(1.0-tanh((x-x0)/d));
       Initial_Conditions::Te.set(i,j,Local_kbTe);
       }  */

  // for DLM distribution 
  double alpha = sqrt(3.0/2.0*Gamma_Func(3.0/equation_switches::Gamma_A)
		      /Gamma_Func(5.0/equation_switches::Gamma_A));
  double Const = equation_switches::Gamma_A/fourpi/alpha/alpha/alpha
    /Gamma_Func(3.0/equation_switches::Gamma_A);
  
  IMPACT_Dim x2(2),x3(3);
  double localEz=0.0;

  //double localf1,curlB;
  for (int i=MPIc->start_WG();i<=MPIc->end_WG();++i) v->Set(i,0.0);
  
  for (int i=MPIc->istart();i<=MPIc->iend();++i)
    {
    for (int j=1;j<=c->Ny();++j)
      {
	// modified to include the temperature rescaling
	Local_kbTe = Initial_Conditions::Te.get(&i,&j);//equation_switches::NEW_T_NORM;

	if (c->NB()>0)
	  for (IMPACT_Dim x1=3;x1>3-c->NB();--x1)
	    B->set(v,Initial_Conditions::B[x1.get()-1].get(&i,&j),&i,&j,&x1);
	//Initialize Ez according to given gradients
	if (c->NE()==3)
	  {
	    localEz=-0.5*Local_kbTe*Initial_Conditions::ni.get(&i,&j)
	      *(IMPACT_Heating::Dnz_xy.get(&i,&j)*IMPACT_Heating::Dnz_t.Get(1)
		+2.5*IMPACT_Heating::DTz_xy.get(&i,&j)
		*IMPACT_Heating::DTz_t.Get(1));
	    E->set(v,localEz,&i,&j,&x3);
	   

	    }
	int one=1;
  double T1m = 2.0/(1.0+Initial_Conditions::T_mult);
  double T2m = Initial_Conditions::T_mult*T1m;

	if (equation_switches::Gamma_A==2.0)
	  if (sqrt(Local_kbTe)>2.0*c->dv(&one)) { // otherwise it is too narrow for grid
	    for (int k=1;k<=c->Nv();++k)
	      {
	      	Local_Maxwellian = Initial_Conditions::ne.get(&i,&j)*(
            (1.0-Initial_Conditions::f0delta)/pow(twopi*T1m*Local_kbTe,1.5)*exp(-(c->v2(&k)/(2.0*T1m*Local_kbTe)))+
            Initial_Conditions::f0delta/pow(twopi*T2m*Local_kbTe,1.5)*exp(-(c->v2(&k)/(2.0*T2m*Local_kbTe)))); // /2.0
	      	f0->set(v,Local_Maxwellian,&i,&j,&k);
      //     Local_Maxwellian = Initial_Conditions::ne.get(&i,&j)
      // /pow(twopi*Local_kbTe,1.5)
      // *exp(-(c->v2(&k)/(2.0*Local_kbTe))); // /2.0
          // f0->set(v,Local_Maxwellian,&i,&j,&k);
	      }
	  } else {
	    f0->set(v,Initial_Conditions::ne.get(&i,&j)/(fourpi*c->v2(&one)*c->dv(&one)),&i,&j,&one);
	    for (int k=2;k<=c->Nv();++k) {
	      f0->set(v,0.0,&i,&j,&k);
	    }
	  } 
	else
	  
	  // Sets to a DLM distribution
	  if (sqrt(Local_kbTe)>2.0*c->dv(&one)) { // otherwise it is too narrow for grid
	    for (int k=1;k<=c->Nv();++k)
	      {
		Local_Maxwellian = Initial_Conditions::ne.get(&i,&j)
		  *Const/pow(Local_kbTe,1.5)
		  *exp(-pow(c->v(&k)/alpha/sqrt(Local_kbTe),
			    equation_switches::Gamma_A)); // /2.0
		f0->set(v,Local_Maxwellian,&i,&j,&k);
	      }
	  } else {
	    f0->set(v,Initial_Conditions::ne.get(&i,&j)/(fourpi*c->v2(&one)*c->dv(&one)),&i,&j,&one);
	    for (int k=2;k<=c->Nv();++k) {
	      f0->set(v,0.0,&i,&j,&k);
	    }
	  }
      }//end of j loop
    }//end of i loop
  //********************************************************************
 
 SwapGhosts(c,MPIc,v);
}
/*double ne_func(int i, int j, IMPACT_Config *c)
  {
  return 1.0;
  }
  double Ue_func(int i, int j, IMPACT_Config *c)
  {
  return 0.5+0.1*sin(1.0*pi*i/c->Nx())*sin(1.0*pi*i/c->Nx());
}
void IMPACT_Math_init(IMPACT_Config *c)
{
  struct IMPACT_Moment ne;
  int Nx=c->Nx();
  int Ny=c->Ny();
  strcpy(ne.name,"n_e");
  strcpy(ne.coords,IMPACT_Coords);
  ne.Nx=Nx;
  ne.Ny=Ny;
  ne.values.ChangeSize(Nx,Ny);
  
  for (int i=1;i<=Nx;++i)
    for (int j=1;j<=Ny;++j)
      ne.values.set(i,j,ne_func(i,j,c));
  GetGrid(c,&ne); // gets x and y axis
  struct IMPACT_Moment Ue=Duplicate(&ne);
  strcpy(Ue.name,"U_e");
  for (int i=1;i<=Nx;++i)
    for (int j=1;j<=Ny;++j)
      Ue.values.set(i,j,Ue_func(i,j,c));
  IMPACT_Write(&ne,"ne_init.imd",IMPACT_Messages::Input_Directory);
  IMPACT_Write(&Ue,"Ue_init.imd",IMPACT_Messages::Input_Directory);
  }*/
void IMPACTA_Initiate_Ci(IMPACT_Config *c)
{
  if (IMPACTA_ions::ion_motion)
    {
      int Nx=c->Nx(),Ny=c->Ny();
      if(IMPACT_Boundaries::fix_Ci[0])
	for (int j=1;j<=Ny;++j)
	  Initial_Conditions::C_i[0].set(1,j,IMPACT_Boundaries::fix_Ci[0]
					 /sqrt(equation_switches::NEW_T_NORM));
      
      if(IMPACT_Boundaries::fix_Ci[1])
	for (int j=1;j<=Ny;++j)
	  Initial_Conditions::C_i[0].set(Nx,j,-IMPACT_Boundaries::fix_Ci[1]
					 /sqrt(equation_switches::NEW_T_NORM));
     if(IMPACT_Boundaries::fix_Ci[2])
	for (int i=1;i<=Nx;++i)
	  Initial_Conditions::C_i[1].set(i,1,IMPACT_Boundaries::fix_Ci[2]
					 /sqrt(equation_switches::NEW_T_NORM));
  
     if(IMPACT_Boundaries::fix_Ci[3])
       for (int i=1;i<=Nx;++i)
	 Initial_Conditions::C_i[1].set(i,Ny,-IMPACT_Boundaries::fix_Ci[3]
					/sqrt(equation_switches::NEW_T_NORM));
     
    }
}
void IMPACTA_Multiply_Ci(IMPACT_Config *c,double value);
void IMPACTA_Rescale_T(IMPACT_Config *c,IMPACT_MPI_Config *M, IMPACT_ParVec *v,
		       IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,
		       IMPACT_Var *f3,IMPACT_Var *E,IMPACT_Var *B,
		       double newT_norm)
{
  if (newT_norm!=1.0)
    {
      IMPACTA_Multiply_Ci(c,1.0/sqrt(newT_norm));
      c->scale_T(newT_norm);
      f0->scale_T(c,M,v,newT_norm);
      f1->scale_T(c,M,v,newT_norm);
      f2->scale_T(c,M,v,newT_norm);
      E->scale_T(c,M,v,newT_norm);
      B->scale_T(c,M,v,newT_norm);
    }
}
double Gamma_Func(double a)
{
  int Nsteps=2000000; // number of steps in the integral
  double dt=1e-5;  // step size in the integral
  double t=dt, ans=0.0;
  if (a==1.0) ans=1.0;
  else if (equation_switches::Gamma_A==2.0) ans=1.0;
  else
    for (int i=0;i<Nsteps;++i)
      {
	//	ans+=exp(-t)*pow(t,a-1.0)*dt;
	t+=dt;
      } 
  return ans;
}
