/*
**********************************************************
Ionization package for IMPACTA

AGRT

March 2011 - very simple to start with - Saha model
according to Drake HEDP
**********************************************************
*/

double Saha_eos(double TkeV, double ne_local);
double TF_eos(double TkeV, double ne_local);

int IMPACTA_Ionization(IMPACT_Config *c, IMPACT_MPI_Config *M,IMPACT_ParVec *v,IMPACT_Var *f0,int n)
{
  double Ztemp=0.0,TkeV = 0.0, ne_local=0.0, Delta_Z=0.0, EoverEl =0.0;
  int one=1;
   int n_adjusted=n;
  if (n_adjusted<1) n_adjusted=1;
  for (int i=M->istart();i<=M->iend();++i) {
    for (int j=1;j<=c->Ny();++j) {
      TkeV = keVnorm*Local_Te(v,f0,c,&i,&j)/c_L/c_L;
	ne_local = Local_ne(v,f0,c,&i,&j);

switch (IMPACTA_ions::ionization_on)
{
	case 1:
	Ztemp = Saha_eos(TkeV,ne_local);
	break;	
	case 2:
	Ztemp = TF_eos(TkeV,ne_local);
	break;
	default:
	break;
} 

    // now add contribution from tunneling ionizarion - laser only at present
      EoverEl =sqrt(IMPACT_Heating::Heating_xy.get(&i,&j)*IMPACT_Heating::Heating_t.Get(n_adjusted)*IMPACT_Heating::vosc_squared);
      if (EoverEl){
	EoverEl=IMPACTA_ions::ionization_v_osc/ EoverEl;
      }
      
      Ztemp += (4.0*IMPACTA_ions::ionization_time_frequency*c->dt()) *EoverEl*exp(-(2.0/3.0)*EoverEl);
      
      if (Ztemp<Z_minimum) Ztemp=Z_minimum;
      if (Ztemp>1.0) Ztemp=1.0;
      /*	
		Initial_Conditions::C_istar[0].set(i,j,Ztemp);  
		}
		}
		
  MPI_Barrier(MPI_COMM_WORLD);
  IMPACTA_Share_Moment(c,M, &Initial_Conditions::C_istar[0]);
  
  // Now apply 5 point smoothing filter to Z: every N_Z_smooth steps
  
  if (!(n%IMPACTA_ions::N_Z_smooth)&&IMPACTA_ions::N_Z_smooth&&n>0) {
  int iplus,iminus,jplus,jminus;
  IMPACT_stencil temp_sten;
  for (int i=M->istart();i<=M->iend();++i) {
  for (int j=1;j<=c->Ny();++j){
	iplus=i+1,iminus=i-1,jplus=j+1,jminus=j-1;
	IMPACT_ni_bound(&temp_sten,c->Nx(),c->Ny(),&iplus,&iminus,&jplus,&jminus);
	if (iplus>c->Nx()) iplus=1; //for periodic bounds
	if (iminus<1) iminus=c->Nx();
	// Using 0.2, 0.2,0.2,0.2,0.2, smoothing stencil
	Ztemp = (Initial_Conditions::C_istar[0].get(&i,&j))*0.2;
	Ztemp += (Initial_Conditions::C_istar[0].get(&iplus,&j))*0.2;
	Ztemp += (Initial_Conditions::C_istar[0].get(&iminus,&j))*0.2;
	Ztemp += (Initial_Conditions::C_istar[0].get(&i,&jplus))*0.2;
	Ztemp += (Initial_Conditions::C_istar[0].get(&i,&jminus))*0.2;
	
	Initial_Conditions::C_istar[1].set(i,j,Ztemp); 
	}
	}
	for (int i=M->istart();i<=M->iend();++i) {
	for (int j=1;j<=c->Ny();++j){
	Initial_Conditions::C_istar[0].set(i,j,Initial_Conditions::C_istar[1].get(&i,&j)); 
	}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	IMPACTA_Share_Moment(c,M, &Initial_Conditions::C_istar[0]);
	
	}
	
	
  // Finally set Z and f0
  for (int i=M->istart();i<=M->iend();++i) {
  for (int j=1;j<=c->Ny();++j){
     
  Ztemp = Initial_Conditions::C_istar[0].get(&i,&j);
  
  if (Ztemp<Z_minimum) Ztemp=Z_minimum;
  if (Ztemp>1.0) Ztemp=1.0;
      */
      Delta_Z = Ztemp - Initial_Conditions::Z.get(&i,&j);
      
      if (Delta_Z>0.0){ // only allow ionization, not recombination
	Initial_Conditions::Z.set(i,j,Ztemp);
	
	// Now increment f0 to add DeltaZ*ni electrons at v=0
	f0->inc(v,Delta_Z*Initial_Conditions::ni.get(&i,&j)/(fourpi*c->v2(&one)*c->dv(&one)),&i,&j,&one);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  IMPACTA_Share_Moment(c,M, &Initial_Conditions::Z);
  return 0;
}

double Saha_eos(double TkeV, double ne_local)
{
	// Saha ionization model
      // see Drake 'High Energy Density Physics"
     double Ztemp = 1.0+0.19*log(pow(TkeV,1.5)/(ne_local*ne24));
      
      if (Ztemp>0.0) {
	Ztemp = (19.7*sqrt(TkeV*(Ztemp))-0.5)*oneover_atomic_Z;
      } else {
	Ztemp = 0.0;
      } 
	return Ztemp;
}

/* EOS 
   Calculation of Zfree using the Thomas-Fermi model (Salzmann, 1998) 
	(Gregori's implemntation - returns DeltaZ)
	I have modified to be the same as in Salzmann as I don't
	understand his 
*/

double TF_eos(double TkeV, double ne_local)
{
  //	double Te,ne;
	double rho,R,T0,alpha,beta;
	double TF,a1,a2,a3,a4,A,b0,b1,b2,B,c1,c2,CC,Q1,Q,x,Ztemp;

	
// constants
const double amu=1.66053873E-27;

// in gm/cc
	rho   = ne_local*ne24*1.0e24*oneover_atomic_Z*IMPACTA_ions::atomic_mass_number*amu*1.0E3;

// Parameters in Salzmann
	R     = rho*oneover_atomic_Z/IMPACTA_ions::atomic_mass_number;
	T0    = TkeV/1000.0*pow(oneover_atomic_Z,4.0/3.0);

	alpha = 14.3139;
	beta  = 0.6624;
	TF    = T0/(1.0+T0);
	a1    = 3.323E-3;
	a2    = 0.971832;
	a3    = 9.26148E-5;
	a4    = 3.10165;
	A     = a1*pow(T0,a2)+a3*pow(T0,a4);
	b0    = -1.7630;
	b1    = 1.43175;
	b2    = 0.315463;
	B     = -exp(b0+b1*TF+b2*pow(TF,7.0));
	c1    = -0.366667;
	c2    = 0.983333;
	CC    = c1*TF+c2;
	Q1    = A*pow(R,B);
	Q     = pow(pow(R,CC)+pow(Q1,CC),1.0/CC);
	x     = alpha*pow(Q,beta);
	
	// old->	ans = Z-ZA*x/(1.0+x+sqrt(1.0+2.0*x));
	
	 Ztemp = x/(oneover_atomic_Z*(1.0+x+sqrt(1.0+2.0*x)));
	
	return Ztemp;
}
