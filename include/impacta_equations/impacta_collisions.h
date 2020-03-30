/*
**********************************************************
Rosenbluth coefficients for rhs collision operators

also heating terms
IMPACT
Version 2.9

AGRT
11/2/08 
Updated to include IB contribution to f2
18/4/07
Updated to be more efficient 6/5/07
8/5/07 Today I will try and make it so that the RB coefficents
are stored (or rather the flux differential stencil) so
that they don't necessarily have to be iterated.

25/9/07 - Includes code for updating Ez from Dn and DT
as form is similar to Heating 

16/10/07 - Update Heating changed so that f0 is updated by
Maxwellian Heating operator if f0 equation switched off

17/10/07 - Updated Maxwellian heating operator as it was wrong before!

july 2011 AGRT - added code for calculating  lambda_ei

July 2015 AGRT - fixed subtle error in delta calculation - 
for small W in Chang-Cooper iteration, the 1/W-1/(exp(W)-1)
calulation can be VERY wrong (not converge to 1/2), presumably because
of round off error. Fixed by adding condition if(W<1e-5) W=0.0;

Also removed exp calculation - pointless (and may be source of error)
**********************************************************

*/
double Calculate_Lambda()
{
  // Shkarofsky definition
  double Constant=sqrt(M_PI/2.0)/omega_p_nuei_n*pow(equation_switches::NEW_T_NORM,1.5);
  double soln=1.0/Constant, old_soln=2.0,delta=1.0;
  const double epsilon=1.0e-12;
    int its=0;
  while ((delta>epsilon)&&its<1000)
    {
      old_soln = soln;
          if (soln>0.0) {
  soln = soln - (Constant*soln*soln-soln*log(soln))/(log(soln)-1.0);
      } else {
  soln = soln - (exp(Constant*soln)-soln)/(Constant*exp(Constant*soln)-1.0);
      }

      delta=fabs(soln-old_soln);
        ++its;  
    }
  return soln;  
}

void update_real_quantities()
{
  globalconsts::Lambda_ei = Calculate_Lambda();
  globalconsts::logLambda = log( Lambda_ei); 
  globalconsts::N_par_Debye =  Lambda_ei/(9.0*oneover_atomic_Z);
  globalconsts::real_T_J = 9.1094e-31/c_L/c_L*2.9979e8*2.9979e8; // real temperature in J
  globalconsts::ne24 = 16.0/9.0*pi*pi*pow(8.8542e-12/1.6022e-19/1.6022e-19,3)/N_par_Debye/N_par_Debye*pow(real_T_J/equation_switches::NEW_T_NORM,3)/1.0e6/1.0e24;
  globalconsts::real_wp = sqrt(1.6022e-19*1.6022e-19*ne24*1.0e6*1.0e24/8.8542e-12/9.1094e-31);
  globalconsts::nuei = sqrt(2.0/pi)*logLambda/ Lambda_ei*real_wp;
 
  IMPACTA_ions::ionization_time_frequency = 4.11e16/nuei; // omega_a normalized to e-i collision time
  IMPACTA_ions::ionization_v_osc = 0.168 *(1.78e15/laser_frequency)*c_L; // 0.168 is a0 for 1 micron laser.
  
  // I don't take into account Z or T dependance as is logarithmic anyway, 
  // and calculation would really slow things down
  globalconsts::ionize_coll_log = log(1e5/c_L/c_L*oneover_atomic_Z)/logLambda; // normalized collision logarithm from Bethe formula
   
}
// Simpler exponential function - accurate enough in the combination in
// the delta function (worst is 0.1%, but <10^-5% for x>30 or x<1)
// this is becasue of the exp-1 combination
//2015 - removed - pointless
/*inline double IMPACT_exp(double x)
{
  return (1+x+0.5*x*x+0.1666666666666667*x*x*x
    +0.04166666666666667*x*x*x*x
    +0.00833333333333333*x*x*x*x*x
    +0.00138888888888889*x*x*x*x*x*x
    +0.0001984126984127 *x*x*x*x*x*x*x
    +2.48015873e-5      *x*x*x*x*x*x*x*x
    +2.75573192e-6      *x*x*x*x*x*x*x*x*x);
    }*/


//warn what is coming
double RB_C(IMPACT_Config *c, IMPACT_ParVec *v, IMPACT_Var *f0,
      int *i,int *j,int *k);
double RB_D(IMPACT_Config *c, IMPACT_ParVec *v, IMPACT_Var *f0,
      int *i,int *j,int *k,double delta);
//Delta function here:
inline double IMPACT_Calc_Cee0_delta(IMPACT_Config *c,double *R_C,double *R_D,
             int *i,int *j,int *k,int *accuracy)
{
  double answer=0.5; int kplus=*k+1;
  //first calculate W
  double W = 0.5*(c->dv(k)+c->dv(&kplus))*(*R_C)/(*R_D);
  //now calculate delta and return
  /*  if (!(*accuracy)) // removed 2015
    answer = 1.0/W - 1.0/(IMPACT_exp(W)-1);
  if (*accuracy) 
  */
      /* 
   I chose 1e-5 because series expansion gives 1/2-W/12+hots.. which is 
   only a 1e-6 error for 1e-5
      */
  if (W>1e-5) {
    answer = 1.0/W - 1.0/(exp(W)-1);
  }
  //check delta is a sensible value - can be removed later
  /* if (answer<0.0||answer>1.0) 
    {
      std::cout <<"IMPACT: Warning - In Calc delta - delta < 0 or > 1\n";
      }*/
  //  if(update) Cee0_delta::delta[*k-1]->set(*i,*j,answer); //update delta
  return answer;
}
//Delta function here:
inline double IMPACT_Calc_Cee0_delta_BC(IMPACT_Config *c,double *R_C,
          double *R_D,int *i,int *j,int *k
          ,int *accuracy)
{
  double answer=0.5; int kplus=*k+1;
  if (*k>0&&*k<c->Nv())
    {
      //first calculate W
      double W = 0.5*(c->dv(k)+c->dv(&kplus))*(*R_C)/(*R_D);
      //now calculate delta and return
      /* if (!(*accuracy)) // removed 2015
  answer = 1.0/W - 1.0/(IMPACT_exp(W)-1);
      if (*accuracy) 
      */

      /* 
   I chose 1e-5 because series expansion gives 1/2-W/12+hots.. which is 
   only a 1e-6 error for 1e-5
      */
      if (W>1.0e-5) { 
  answer = 1.0/W - 1.0/(exp(W)-1);
      }
       //check delta is a sensible value - can be removed later
      /*if (answer<0.0||answer>1.0) 
  {
    std::cout <<"IMPACT: Warning - In Calc delta - delta < 0 or > 1\n";
    }*/
      // if(update) Cee0_delta::delta[*k-1]->set(*i,*j,answer); //update delta
    }
  return answer;
}
  //Rosenbluth C coefficient (k+1/2)
inline double RB_C(IMPACT_Config *c, IMPACT_ParVec *v, IMPACT_Var *f0,
      int *i,int *j,int *k)
{
  // C_(k+1/2) = 4pi *sum_l=0^k (f0_l) v_l^2 deltav_l
  double answer=0.0;
  for (int l=1;l<=*k;++l)
    answer+=f0->get(v,i,j,&l)*c->v2(&l)*c->dv(&l);
  answer*=fourpi;
  return answer;
}
//Boundary checking variant
inline double RB_C_BC(IMPACT_Config *c, IMPACT_ParVec *v, IMPACT_Var *f0,
      int *i,int *j,int *k)
{
  // C_(k+1/2) = 4pi *sum_l=0^k (f0_l) v_l^2 deltav_l
  double answer=0.0;
  if (*k>0&&*k<=c->Nv())
    {
      for (int l=1;l<=*k;++l)
  answer+=f0->get(v,i,j,&l)*c->v2(&l)*c->dv(&l);
      answer*=fourpi;
    }
  return answer;
}

//Rosenbluth D coefficient (k+1/2)
inline double RB_D(IMPACT_Config *c, double *f0,int *i,int *j,int *k,
       double *delta)
{
  double answer=0.0;//,innerloopsoln=0.0;
  int kkplus,kplus=*k+1; //m + 1,k+1
  double oneminusdelta=1-*delta;
  double *theinnerequation;
  theinnerequation=new double[c->Nv()];
  
  int kk=c->Nv()-1; kkplus=kk+1;
  theinnerequation[kk]=((oneminusdelta*f0[kk]+(*delta)*f0[kk-1])
      *(c->v2(&kkplus)-c->v2(&kk)));
  for (kk=c->Nv()-2;kk>0;--kk)
    {
      kkplus=kk+1;
      theinnerequation[kk]=theinnerequation[kkplus]
  +((oneminusdelta*f0[kk]+(*delta)*f0[kk-1])
    *(c->v2(&kkplus)-c->v2(&kk)));
    }
  //Outer loop (sum)
  for (int l=1;l<=*k;++l)
    {
      answer+=c->v2(&l)*c->dv(&l)*theinnerequation[l];
    }
  answer*=fourpi;
  answer*=1.0/(c->v(k)+c->v(&kplus));
  delete[] theinnerequation; 
  return answer;
}
//Boundary checking
inline double RB_D_BC(IMPACT_Config *c, double *f0,
          int *i,int *j,int *k,double *delta)
{
  //Dk+1/2=4pi/vk+1/2 sum_k vk^2 *(sum_k f0 vm+1/2 deltavm+1/2)*deltav_k
  double answer=0.0;//,innerloopsoln=0.0;
  int kkplus,kplus=*k+1; //m + 1,k+1
  double oneminusdelta=1-*delta;
  double *theinnerequation;
  theinnerequation=new double[c->Nv()];

  if (*k>0&&*k<c->Nv())
    {
      int kk=c->Nv()-1; kkplus=kk+1;
      theinnerequation[kk]=((oneminusdelta*f0[kk]+(*delta)*f0[kk-1])
              *(c->v2(&kkplus)-c->v2(&kk)));
      for (kk=c->Nv()-2;kk>0;--kk)
  {
    kkplus=kk+1;
        theinnerequation[kk]=theinnerequation[kkplus]
        +((oneminusdelta*f0[kk]+(*delta)*f0[kk-1])
        *(c->v2(&kkplus)-c->v2(&kk)));
  }
      //Outer loop (sum)
      for (int l=1;l<=*k;++l)
  {
     answer+=c->v2(&l)*c->dv(&l)*theinnerequation[l];
  }
      answer*=fourpi;
      answer*=1.0/(c->v(k)+c->v(&kplus));
    }
  delete[] theinnerequation;
  return answer;
}
inline void IMPACT_Iterate_Delta(IMPACT_Config *c, double *f0,int *i, int *j, 
         int *k, double *C, double *D, double *Delta) 
{
  int iterations=0;
  int iteration_check=0; //to check condition
  int accuracy=0; 
  *Delta = 0.5;
  double Dold=*D; //D from last iteration
  do
    {
      //first find D from delta
      *D = RB_D(c,f0,i,j,k,Delta);
      if (fabs(*D-Dold)<zerotolerance::RB_D_tolerance*(1+fabs(*D+Dold))) 
  iteration_check=1;
      //now get new delta
      *Delta = IMPACT_Calc_Cee0_delta(c,C,D,i,j,k,&accuracy);
      ++iterations;
      if (iterations>zerotolerance::RB_D_itmax/2) accuracy=1;
      if (iterations>zerotolerance::RB_D_itmax) iteration_check=1;
  /*  {
  std::cout<<"IMPACT: ERROR - Exceeded maximum delta iterations\n";
  std::cout<<"its = "<<iterations<<",k = "<<*k<<", deltaD = "<<*D-Dold<<'\n';
  exit(0);
  }*/
      Dold = *D;
    } while(!iteration_check);
  //for diagnose average iterations:
  IMPACT_Diagnostics::total_delta_its+=iterations;
  IMPACT_Diagnostics::delta_times+=1.0;
}
inline void IMPACT_Iterate_Delta_BC(IMPACT_Config *c, double *f0,int *i,
            int *j, int *k, double *C, double *D, 
            double *Delta) 
{
  int iterations=0;
  int iteration_check=0; //to check condition
  int accuracy =0;
  *Delta = 0.5;
  double Dold=*D; //D from last iteration
   do
    {
      //first find D from delta
      *D = RB_D_BC(c,f0,i,j,k,Delta);
 
      if (fabs(*D-Dold)<zerotolerance::RB_D_tolerance*(1+fabs(*D+Dold))) 
  iteration_check=1;
      //now get new delta
      *Delta = IMPACT_Calc_Cee0_delta_BC(c,C,D,i,j,k,&accuracy);
      ++iterations;
      if (iterations>zerotolerance::RB_D_itmax/2) accuracy=1;
      if (iterations>zerotolerance::RB_D_itmax) iteration_check=1;
  /*  {
    std::cout<<"IMPACT: ERROR - Exceeded maximum delta iterations\n";
    std::cout<<"its = "<<iterations<<",k = "<<*k<<", deltaD = "<<*D-Dold<<'\n';
    exit(0);
    }*/
      Dold = *D;
      } while(!iteration_check);
  //for diagnose average iterations:
  IMPACT_Diagnostics::total_delta_its+=iterations;
  IMPACT_Diagnostics::delta_times+=1.0;
}
// This function gets the flux F and iterates delta coefficient
inline IMPACT_Vel_Sten Cee0_Flux(IMPACT_Config *c, IMPACT_ParVec *v, 
            IMPACT_Var *f0,double *f0_k_temp,
            int *i, int *j, int *k)
{
  int kplus=*k+1;
  double C=0.0,D=0.0,delta=0.5;
  
  //first get C, D and delta
  if (equation_switches::Cee0_on)
    {
      C=RB_C(c,v,f0,i,j,k);
      IMPACT_Iterate_Delta(c,f0_k_temp,i,j,k, &C, &D, &delta);
    }

 
  C*=oneover_atomic_Z;
  D*=oneover_atomic_Z;
 
   //Now Heating term
    D=D*equation_switches::Cee0_on+(IMPACT_Heating::IB_type_on
  *IMPACT_Heating::Total_heating_contribution.get(i,j)/c->v(k)
      +IMPACT_Heating::MX_type_on
  *IMPACT_Heating::Total_heating_contribution.get(i,j)*c->v2(k)*pibysix);
  //now work out coefficents
  double Doverdv_k_plus_half = 2.0*D/(c->dv(k)+c->dv(&kplus));
  double f0kplus1 = C*(1-delta)+Doverdv_k_plus_half;
  double f0k = C*delta-Doverdv_k_plus_half;
 
  //make stencil
  IMPACT_Vel_Sten answer(0.0,f0k,f0kplus1);
  return answer;
}
// This function gets the flux F and iterates delta coefficient at boundary
inline IMPACT_Vel_Sten Cee0_Flux_BC(IMPACT_Config *c, IMPACT_ParVec *v, 
            IMPACT_Var *f0,double *f0_k_temp,
            int *i, int *j, int *k)
{
  int kplus=*k+1;
 
  double C=0.0,D=0.0,delta=0.5;
  
  //first get C, D and delta
  if (equation_switches::Cee0_on)
    {
      C=RB_C_BC(c,v,f0,i,j,k); 
      IMPACT_Iterate_Delta_BC(c,f0_k_temp,i,j,k, &C, &D, &delta);
    }

  C*=oneover_atomic_Z;
  D*=oneover_atomic_Z;
 
  //Now Heating term
  if (*k>0)
    D=D*equation_switches::Cee0_on+(IMPACT_Heating::IB_type_on
  *IMPACT_Heating::Total_heating_contribution.get(i,j)/c->v(k)
      +IMPACT_Heating::MX_type_on
      *IMPACT_Heating::Total_heating_contribution.get(i,j)*c->v2(k)*pibysix);

   if (kplus>c->Nv()) --kplus;
  //now work out coefficents
   double Doverdv_k_plus_half = 2.0*D/(c->dv(k)+c->dv(&kplus));
  
  double f0kplus1 = C*(1-delta)+Doverdv_k_plus_half;
  double f0k = C*delta-Doverdv_k_plus_half;
 
    //make stencil
  IMPACT_Vel_Sten answer(0.0,f0k,f0kplus1);
  
  return answer;
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
// f2 bits
inline IMPACT_Vel_Sten IB_Flux_f2(IMPACT_Config *c, IMPACT_ParVec *v, 
           IMPACT_Var *f0, int *i, int *j, int *k)
{
  int kplus=*k+1;
  double D=0.0;
   //Now Heating term
  if (*k>0)
    D=3.0*IMPACT_Heating::Total_heating_contribution.get(i,j)/c->v2(k)/c->v2(k);
 
  //now work out coefficents
  double Doverdv_k_plus_half = 2.0*D/(c->dv(k)+c->dv(&kplus));
  
  double f0kplus1 = Doverdv_k_plus_half;
  double f0k = -Doverdv_k_plus_half;
  
  //make stencil
  IMPACT_Vel_Sten answer(0.0,f0k,f0kplus1);
  
  return answer;
}
// f2 bits
inline IMPACT_Vel_Sten IB_Flux_f2_BC(IMPACT_Config *c, IMPACT_ParVec *v, 
           IMPACT_Var *f0, int *i, int *j, int *k)
{
  int kplus=*k+1;
  double D=0.0;
   //Now Heating term
  if (*k>0)
    D=3.0*IMPACT_Heating::Total_heating_contribution.get(i,j)/c->v2(k)/c->v2(k);
  if (kplus>c->Nv()) --kplus;
  //now work out coefficents
  double Doverdv_k_plus_half = 2.0*D/(c->dv(k)+c->dv(&kplus));
  
  double f0kplus1 = Doverdv_k_plus_half;
  double f0k = -Doverdv_k_plus_half;
  
  //make stencil
  IMPACT_Vel_Sten answer(0.0,f0k,f0kplus1);
  
  return answer;
}
//_-___-___------___---_-_-------____---_-----______---_-_-_-----_

inline IMPACT_Vel_Sten Cee0_dF(IMPACT_Config *c, IMPACT_ParVec *v, 
          IMPACT_Var *f0,double *f0_k_temp,
          int *i, int *j, int *k)
{
  int kminus=*k-1;
  double Cee0_const = -c->dt()/c->v2(k)*c->idv(k);
  IMPACT_Vel_Sten Fkplushalf = Cee0_Flux(c,v,f0,f0_k_temp,i,j,k);
  IMPACT_Vel_Sten Fkminushalf = Cee0_Flux(c,v,f0,f0_k_temp,i,j,&kminus);
  Fkplushalf*(Cee0_const);
  Fkminushalf*(-Cee0_const);
  
  IMPACT_Vel_Sten Fk(0.0,0.0,0.0);
  Fk.inc(0,Fkminushalf.get(1));
  Fk.inc(1,Fkminushalf.get(2)+Fkplushalf.get(1));
  Fk.inc(2,Fkplushalf.get(2));
  Fk*Cee0_const;
  return Fk;
}
inline IMPACT_Vel_Sten Cee0_dF_BC(IMPACT_Config *c, IMPACT_ParVec *v, 
          IMPACT_Var *f0,double *f0_k_temp,
          int *i, int *j, int *k)
{
  int kminus=*k-1;
  double Cee0_const = -c->dt()/c->v2(k)*c->idv(k);
  IMPACT_Vel_Sten Fkplushalf = Cee0_Flux_BC(c,v,f0,f0_k_temp,i,j,k);
  IMPACT_Vel_Sten Fkminushalf = Cee0_Flux_BC(c,v,f0,f0_k_temp,i,j,&kminus);
  Fkplushalf*(Cee0_const);
  Fkminushalf*(-Cee0_const);
  
  IMPACT_Vel_Sten Fk(0.0,0.0,0.0);
  Fk.inc(0,Fkminushalf.get(1));
  Fk.inc(1,Fkminushalf.get(2)+Fkplushalf.get(1));
  Fk.inc(2,Fkplushalf.get(2));
  Fk*Cee0_const;
  return Fk;
}
inline void IMPACT_Update_Cee0(IMPACT_Config *c, IMPACT_ParVec *v, 
             IMPACT_Var *f0,int *i, int *j)
{
  double *f0_k_temp;
  f0_k_temp=new double[c->Nv()+1];
  f0->Extract_k(f0_k_temp,v,i,j);
  IMPACT_Vel_Sten tempvs;
  int k=1;
  tempvs=Cee0_Flux_BC(c,v,f0,f0_k_temp,i,j,&k);
  Cee0_flux_store::flux->set(&tempvs,i,j,&k);
  // now for f2
  tempvs=IB_Flux_f2_BC(c,v,f0,i,j,&k);
  Cee0_flux_store::f2_IB->set(&tempvs,i,j,&k);

  for (k=2;k<c->Nv();++k)
    {
      tempvs=Cee0_Flux(c,v,f0,f0_k_temp,i,j,&k);
      Cee0_flux_store::flux->set(&tempvs,i,j,&k);
      // now for f2
      tempvs=IB_Flux_f2(c,v,f0,i,j,&k);
      Cee0_flux_store::f2_IB->set(&tempvs,i,j,&k);
    }

  k=c->Nv();
  tempvs=Cee0_Flux_BC(c,v,f0,f0_k_temp,i,j,&k);
  Cee0_flux_store::flux->set(&tempvs,i,j,&k);
  // now for f2
  tempvs=IB_Flux_f2_BC(c,v,f0,i,j,&k);
  Cee0_flux_store::f2_IB->set(&tempvs,i,j,&k);
  delete[] f0_k_temp;
}

inline void IMPACT_Cee0(IMPACT_Config *c, IMPACT_MPI_Config *M,
      IMPACT_ParVec *v, IMPACT_Var *f0,int lagged_step)
{

      int check =1;
      if (!zerotolerance::RB_iterate&&lagged_step>0) check=0;
      if (check&&equation_switches::f0_equation_on)
  {
    if (!M->rank())
      std::cout<<"\nIMPACT: Calculating Cee0 terms";
    double time1=clock();
    for (int i=M->istart();i<=M->iend();++i)
      for (int j=1;j<=c->Ny();++j)
        IMPACT_Update_Cee0(c,v,f0,&i,&j);
    double time2=clock();
    time1=(time2-time1)/CLOCKS_PER_SEC;
    if (!M->rank())
      std::cout<<"                            - took "<<IMPACT_GetTime(time1)<<'\n';
  }

}

inline void IMPACT_Maxwellianf0evolve(IMPACT_ParVec *v, IMPACT_Config *c,
              IMPACT_MPI_Config *M, IMPACT_StenOps *O, 
              IMPACT_Var *f0, IMPACT_Var *f1, 
              IMPACT_Var *E)
{
     //first work out deltaT - change in temperature
      // deltaT=4/3*deltaUe/ne = 4/3*dt*(div q + j.E + dU/dt|heating)
    
      //At the moment only works with reflective bounds
      IMPACT_Dim x1(1),x2(2),x3(3);
      double deltaT=0.0,divq=0.0,jdotE=0.0,T=0.0,ne=0.0,heatingdT=0.0;
      int iplus,iminus,jplus,jminus;
      IMPACT_stencil tempsten;
      for (int i=M->istart();i<=M->iend();++i)
  for (int j=1;j<=c->Ny();++j)
    { 
      
      iplus=i+1; iminus=i-1; jplus=j+1; jminus=j-1;
      tempsten = *(O->div(&i,&j));
      IMPACT_f1_E_bound(&tempsten,c->Nx(),c->Ny(),&iplus,&iminus,
            &jplus,&jminus,&x1);
      /*if (iplus<M->istart()||iminus>M->iend())
        {
    std::cout<<"IMPACT: ERROR - EXMX does not yet work in parallel with periodic bounds\n";
    exit(0);
    }*/
      divq=Local_qT(v,f1,c,&iplus,&j,&x1)*tempsten(1)
        +Local_qT(v,f1,c,&iminus,&j,&x1)*tempsten(2);
      iplus=i+1; iminus=i-1; jplus=j+1; jminus=j-1;
      tempsten = *(O->div(&i,&j));
      IMPACT_f1_E_bound(&tempsten,c->Nx(),c->Ny(),&iplus,&iminus,
            &jplus,&jminus,&x2);
      /*if (iplus<M->istart()||iminus>M->iend())
        {
    std::cout<<"IMPACT: ERROR - EXMX does not yet work in parallel with periodic bounds\n";
    exit(0);
    }*/
      divq +=Local_qT(v,f1,c,&jplus,&j,&x2)*tempsten(3)
        +Local_qT(v,f1,c,&jminus,&j,&x2)*tempsten(4);
      
      if (c->NE()>0&&c->Nf1()>0)
        jdotE =E->get(v,&i,&j,&x1)*Local_je(v,f1,c,&i,&j,&x1);
      if (c->NE()>1&&c->Nf1()>1)
        jdotE+=E->get(v,&i,&j,&x2)*Local_je(v,f1,c,&i,&j,&x2);
      if (c->NE()>2&&c->Nf1()>2)      
        jdotE+=E->get(v,&i,&j,&x3)*Local_je(v,f1,c,&i,&j,&x3);
      T=Local_Te(v,f0,c,&i,&j);
      heatingdT=(2*IMPACT_Heating::IB_type_on
           *IMPACT_Heating::Total_heating_contribution.get(&i,&j)
           /T/sqrt(T)*0.85
           +2*IMPACT_Heating::MX_type_on
           *IMPACT_Heating::Total_heating_contribution.get(&i,&j));
      deltaT=4/3*c->dt()*(heatingdT-divq-jdotE);
      T=T+deltaT;
      ne=Local_ne(v,f0,c,&i,&j);
      // now update f0 with new maxwellian

      for (int k=1;k<=c->Nv();++k)
        f0->set(v,ne/pow(twopi*T*0.5,1.5)*exp(-(c->v2(&k)/T)),&i,&j,&k);
    } //end of i j loop
}
//update heating matrix temporally
inline void IMPACT_update_heating(IMPACT_ParVec *v, IMPACT_Config *c,
          IMPACT_MPI_Config *M, IMPACT_StenOps *O, 
          IMPACT_Var *f0, IMPACT_Var *f1, 
          IMPACT_Var *E, int n)
{
  int n_adjusted=n;
  if (n_adjusted<1) n_adjusted=1;
 
  if (!equation_switches::tr_bool){
  IMPACT_Heating::Total_heating_contribution.copy(&IMPACT_Heating::Heating_xy);}
  else{

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  IMPACT_Heating::Total_heating_contribution.copy(&IMPACT_Heating::i_mat); // Ray Tracing LINE

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  }
   IMPACT_Heating::Total_heating_contribution.SwitchOffConst(c->Nx(),c->Ny());
   IMPACT_Heating::Total_heating_contribution.multiplyall
    (IMPACT_Heating::Heating_t.Get(n_adjusted)*IMPACT_Heating::vosc_squared);
  double heating_val=0.0;
    
  for (int i=1;i<=c->Nx();++i)
    for (int j=1;j<=c->Ny();++j)
      {
  heating_val=Initial_Conditions::Z.get(&i,&j)*
    Initial_Conditions::Z.get(&i,&j)*
    Initial_Conditions::ni.get(&i,&j)*oneover6;
  IMPACT_Heating::Total_heating_contribution.Iset(IMPACT_Heating::Total_heating_contribution.get(&i,&j)*heating_val,&i,&j);
  }
  
   // Now if f0 equation is off, and heating maxwellian, update distribution
  // due to heat flow and heating.
  if (!equation_switches::f0_equation_on&&equation_switches::evolvef0 )
    {
      IMPACT_Maxwellianf0evolve(v,c,M,O,f0,f1,E);
    }
}

//update dn/dz matrix temporally
inline void IMPACT_update_Ez(IMPACT_Config *c, IMPACT_MPI_Config *MPIc,
           IMPACT_Var *f0,IMPACT_Var *E,
           IMPACT_ParVec *v,int *n)
{
  /*double localEz=0;
  IMPACT_Dim x3(3);
  double deltan=IMPACT_Heating::Dnz_t.Get(n+1)-IMPACT_Heating::Dnz_t.Get(n);
  double deltaT=IMPACT_Heating::DTz_t.Get(n+1)-IMPACT_Heating::DTz_t.Get(n);
  if (deltan!=0.0 && deltaT!=0.0)
    {
      for (int i=MPIc->istart();i<=MPIc->iend();++i)
  for (int j=1;j<=c->Ny();++j)
    {
      localEz=-0.5*(Local_Te(v,f0,c,&i,&j)*
        IMPACT_Heating::Dnz_xy.get(&i,&j)*deltan
        +1.5*Initial_Conditions::ni.get(&i,&j)*
        IMPACT_Heating::DTz_xy.get(&i,&j)*deltaT);
      E->inc(v,localEz,&i,&j,&x3);
    }
    }
    }*/
  int nadd=*n;
  if (nadd<1) nadd=1;
  double gridtemp=0.0; 
  for (int i=MPIc->istart();i<=MPIc->iend();++i)
    for (int j=1;j<=c->Ny();++j)
      {
  Initial_Conditions::ne.Iset(Local_ne(v,f0,c,&i,&j),&i,&j);
  Initial_Conditions::Te.Iset(Local_Te(v,f0,c,&i,&j),&i,&j);
  gridtemp=IMPACT_Heating::Dnz_xy.get(&i,&j)
    *IMPACT_Heating::Dnz_t.Get(nadd);
  // /Initial_Conditions::ne.get(&i,&j);
  
  IMPACT_Heating::Dnbydzgrid.Iset(gridtemp,&i,&j);
  
  gridtemp=IMPACT_Heating::DTz_xy.get(&i,&j)
    *IMPACT_Heating::DTz_t.Get(nadd);
  // /Initial_Conditions::Te.get(&i,&j);
  
  IMPACT_Heating::DTbydzgrid.Iset(gridtemp,&i,&j);

     }
}
