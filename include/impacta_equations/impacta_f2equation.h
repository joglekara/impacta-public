/*
**********************************************************
 f2 equation IMPACTA code for inner and outer cells


AGRT
4/6/07
11/2/08 - IB Heating term added in 
24/5/08 - IB Heating term modified to include stat B field effect
 NB NB NB NB NB So far this only includes the effects of Bz

2/7/08 - NB NB Weibel damping included - i.e. diffusion of f2

**********************************************************

*/

inline void f2equation_inner(IMPACT_Config *c, IMPACT_ParVec *vlagged,IMPACT_ParVec *vtemp, IMPACT_StenOps * O,IMPACT_Var* f0,IMPACT_Var* f1,IMPACT_Var* f2,IMPACT_Var* f3,IMPACT_Var* E,IMPACT_Var* B,int *i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  /*
**********************************************************
INSERT f2Equation terms CODE HERE
**********************************************************
   */
  f2->set(vtemp,equation_switches::e_viscosity_on*c->idt(),i,j,k,x1,x2); 
  // set f2 n+1 element to 1.0

  //Collision part...
  double ei_cols=Initial_Conditions::Z.get(i,j)*Initial_Conditions::Z.get(i,j)*Initial_Conditions::ni.get(i,j)/c->v3(k);
  f2->inc(vtemp,3.0*equation_switches::Cei_on*ei_cols,i,j,k,x1,x2); 
  
 // Gradient of f1 terms (3 of them)
  //first v/2(dfj/dxi + dfi/dxj)
  double multiplier=c->v(k)*0.5;
  IMPACT_stencil temp_sten = (*O->ddxi(i,j,x1))*multiplier;
  f1->insert(&temp_sten,vtemp,i,j,k,x2);
  temp_sten = (*O->ddxi(i,j,x2))*multiplier;
  f1->insert(&temp_sten,vtemp,i,j,k,x1);
  //v/3 div f1 delta_ij
  multiplier=-c->v(k)*oneover3*Kronecker_Delta(x1,x2);
  for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
    {
      temp_sten=(*O->ddxi(i,j,&x3))*multiplier;
      f1->insert(&temp_sten,vtemp,i,j,k,&x3);
    }

  //_________________________________________
  //Now velocity gradient terms
  multiplier=c->v(k)*0.5;
  //v/2aid/dvfj/v + v/2ajd/dvfi/v
  E->inc(vtemp,-multiplier*f1->ddv1overv(vlagged, c,i,j,k,x1),i,j,x2);
  E->inc(vtemp,-multiplier*f1->ddv1overv(vlagged, c,i,j,k,x2),i,j,x1);
  //v/2 .2/3 .deltaij a.df1/dv
  multiplier=-oneover3*c->v(k)*Kronecker_Delta(x1,x2);
  if (c->NE()>0)
    for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
      {
	E->inc(vtemp,-multiplier*f1->ddv1overv(vlagged,c,i,j,k,&x3),i,j,&x3);
      }

  //Now B terms... complicated
  /*
    wB_x4 (eps_x1,x4,x3 f2_x2
  */
  if (c->NB()>0)
    for (IMPACT_Dim x4=3;x4>3-c->NB();--x4)
      {	
	for (IMPACT_Dim x3=1;x3<=c->N3f2();++x3)
	  {
	    // wB_k epsilon_ikn f2_jn j is x2 in this case
	    //   if (x3.get()+x2->get()<6)
	    f2->inc(vtemp,-B->get(vlagged,i,j,&x4)*CeviLevita(x1,&x4,&x3),
		    i,j,k,&x3,x2);
	    
	    // wB_k epsilon_jkn f2_in i is x1 in this case
	    //     if (x3.get()+x1->get()<6)
	    f2->inc(vtemp,-B->get(vlagged,i,j,&x4)*CeviLevita(x2,&x4,&x3)
		    ,i,j,k,&x3,x1);
	  } //end of x3 loop
	
      }
  //************************************************************
  //IB heating term
  // Polarization matrix terms:
  //-----------------------------
 //  double Pol_E=IMPACT_Heating::vosc_hat.get(x1->get(),x2->get());
 //  IMPACT_Dim xdir(1);
 //  double Pol_B=0.0;
 // if (c->NB()>1)
 //   {
 //     Pol_B=IMPACT_Heating::vosc_hat_Bx.get(x1->get(),x2->get())
 //       *B->get(vlagged,i,j,&xdir);
 //     xdir=2;
 //     Pol_B+=IMPACT_Heating::vosc_hat_By.get(x1->get(),x2->get())
 //       *B->get(vlagged,i,j,&xdir);
 //   }
 // xdir=3;
 // if (c->NB()>0)
 //   Pol_B+=IMPACT_Heating::vosc_hat_Bz.get(x1->get(),x2->get())
 //     *B->get(vlagged,i,j,&xdir);
 // Pol_B/=ei_cols;
 
 //  //-----------------------------
 //  int kminus=*k-1;
 //  double Cee0_const = -1.0*(Pol_E+Pol_B)*c->v(k)*c->idv(k);
 //  IMPACT_Vel_Sten Fk(0.0,0.0,0.0);
 //  Fk.inc(0,-Cee0_flux_store::f2_IB->get(i,j,&kminus)->get(1));
 //  Fk.inc(1,-Cee0_flux_store::f2_IB->get(i,j,&kminus)->get(2)
	//  +Cee0_flux_store::f2_IB->get(i,j,k)->get(1));
 //  Fk.inc(2,Cee0_flux_store::f2_IB->get(i,j,k)->get(2));
 //  Fk*Cee0_const;
 //  // Now add final term for B effect on IB
 //  Fk.inc(1,-3.0*Pol_B*Cee0_flux_store::f2_IB->get(i,j,k)->get(1));
 //  Fk.inc(2,-3.0*Pol_B*Cee0_flux_store::f2_IB->get(i,j,k)->get(2));

 //  f0->insert(&Fk,vtemp,i,j,k);

 //************************************************************

  /*
    Diffusion of f2 - obtained by substituting f2 into f3
    equation, then back into f2 equation
  */
  if (equation_switches::df3_by_dt_on)
    {
      IMPACT_Dim xtemp;
      GetOrthogonal(x1,&xtemp);
      IMPACT_Dim xdim(1),ydim(2);
      int iplus=*i+1,iminus=*i-1;
      double idx = 1.0/(c->dx(&iplus)+c->dx(&iminus));
      double diffusion_const=-c->v2(k)/(42.0*ei_cols);
      temp_sten = O->Laplacian(i,j);
      temp_sten = temp_sten*diffusion_const;
      f2->insert(&temp_sten,vtemp,i,j,k,x1,x2);
      if (x1->get()<3&&x2->get()<3)
	{
	  if (x1->get()==x2->get())
	    {
	      
	      temp_sten = (*O->d2dx2i(i,j,x1));
	      temp_sten = (temp_sten*diffusion_const)*0.8;
	      f2->insert(&temp_sten,vtemp,i,j,k,x1,x1);
	      temp_sten = (*O->d2dx2i(i,j,&xtemp));
	      temp_sten = (temp_sten*diffusion_const)*(-0.4);
	      f2->insert(&temp_sten,vtemp,i,j,k,&xtemp,&xtemp);
	      // Now dxdy derivative
	      temp_sten = (*O->ddy(j))*(0.4*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iplus,j,k,&xdim,&ydim);
	      temp_sten = (*O->ddy(j))*(-0.4*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iminus,j,k,&xdim,&ydim);	
	    }
	  else 
	    {
	      temp_sten = (*O->d2dx2i(i,j,&xdim))+(*O->d2dx2i(i,j,&ydim));
	      temp_sten = (temp_sten*diffusion_const)*0.6;
	      f2->insert(&temp_sten,vtemp,i,j,k,x1,x2);
	      // Now dxdy derivative
	      temp_sten = (*O->ddy(j))*(0.6*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iplus,j,k,&xdim,&xdim);
	      f2->insert_BC(&temp_sten,vtemp,&iplus,j,k,&ydim,&ydim);
	      temp_sten = (*O->ddy(j))*(-0.6*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iminus,j,k,&xdim,&xdim);
	      f2->insert_BC(&temp_sten,vtemp,&iminus,j,k,&ydim,&ydim);
	    }
	}
      else
	{
	  if (x1->get()==3) temp_sten = (*O->d2dx2i(i,j,x2));
	  else temp_sten = (*O->d2dx2i(i,j,x1));
	  temp_sten = (temp_sten*diffusion_const)*0.6;
	  f2->insert(&temp_sten,vtemp,i,j,k,x1,x2);
	}
    }
  
  //*******************************************
   if (IMPACTA_ions::ion_motion)
    { 
      // Ion motion terms
      // + C.grad f1
      for (IMPACT_Dim x4=1;x4<=c->Nf1();++x4)
	{
	  temp_sten = (*O->ddxi(i,j,&x4))*
	    Initial_Conditions::C_i[x4.get()-1].get(i,j);
	  f2->insert(&temp_sten,vtemp,i,j,k,x1,x2);
	}
    }
}
//code for outer cells of f0 equation - i.e.  boundaries
inline void f2equation_outer(IMPACT_Config *c, IMPACT_ParVec *vlagged,IMPACT_ParVec * vtemp, IMPACT_StenOps * O, IMPACT_Var* f0,IMPACT_Var* f1,IMPACT_Var* f2,IMPACT_Var* f3,IMPACT_Var* E,IMPACT_Var* B,int* i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  /*
**********************************************************
INSERT f2 Equation terms CODE HERE - BC versions
**********************************************************
*/
   f2->set(vtemp,equation_switches::e_viscosity_on*c->idt(),i,j,k,x1,x2); 
  // set f2 n+1 element to 1.0

  //Collision part...
  double ei_cols=Initial_Conditions::Z.get(i,j)
    *Initial_Conditions::Z.get(i,j)*Initial_Conditions::ni.get(i,j)/c->v3(k);
  f2->inc(vtemp,3.0*equation_switches::Cei_on*ei_cols,i,j,k,x1,x2); 
  

  // Gradient of f1 terms (3 of them)
  //first v/2(dfj/dxi + dfi/dxj)
  double multiplier=c->v(k)*0.5;
  IMPACT_stencil temp_sten = (*O->ddxi(i,j,x1))*multiplier;

  f1->insert_f1BC(&temp_sten,vtemp,i,j,k,x2);
  temp_sten = (*O->ddxi(i,j,x2))*multiplier;
  f1->insert_f1BC(&temp_sten,vtemp,i,j,k,x1);
  //v/3 div f1 delta_ij
  multiplier=-c->v(k)*oneover3*Kronecker_Delta(x1,x2);
  for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
    {
      temp_sten=(*O->ddxi(i,j,&x3))*multiplier;
      f1->insert_f1BC(&temp_sten,vtemp,i,j,k,&x3);
    }
    
  //_________________________________________
  //Now velocity gradient terms
  multiplier=c->v(k)*0.5;
  //v/2aid/dvfj/v + v/2ajd/dvfi/v
  E->inc(vtemp,-multiplier*f1->ddv1overv_BC(vlagged, c,i,j,k,x1),i,j,x2);
  E->inc(vtemp,-multiplier*f1->ddv1overv_BC(vlagged, c,i,j,k,x2),i,j,x1);
  //v/2 .2/3 .deltaij a.df1/dv
  multiplier=-c->v(k)*oneover3*oneover3*Kronecker_Delta(x1,x2);
  if (c->NE()>0)
    for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
      {
	E->inc(vtemp,-multiplier*f1->ddv1overv_BC(vlagged,c,i,j,k,&x3),i,j,&x3);
      }
  
  //Now B terms... complicated
  /*
    wB_x4 (eps_x1,x4,x3 f2_x2
  */
  if (c->NB()>0)
    for (IMPACT_Dim x4=3;x4>3-c->NB();--x4)
      {	
	for (IMPACT_Dim x3=1;x3<=c->N3f2();++x3)
	  {
	    // wB_k epsilon_ikn f2_jn j is x2 in this case
	    //   if (x3.get()+x2->get()<6)
	    f2->inc(vtemp,-B->get(vlagged,i,j,&x4)*CeviLevita(x1,&x4,&x3),
		    i,j,k,&x3,x2);
	    
	    // wB_k epsilon_jkn f2_in i is x1 in this case
	    //     if (x3.get()+x1->get()<6)
	    f2->inc(vtemp,-B->get(vlagged,i,j,&x4)*CeviLevita(x2,&x4,&x3)
		    ,i,j,k,&x3,x1);
	  } //end of x3 loop
	
      }
  //-----------------------------
 //  IB heating term
  // Polarization matrix terms:
 // -----------------------------
 //  double Pol_E=IMPACT_Heating::vosc_hat.get(x1->get(),x2->get());
 //  IMPACT_Dim xdir(1);
 //  double Pol_B=0.0;
 // if (c->NB()>1)
 //   {
 //     Pol_B=IMPACT_Heating::vosc_hat_Bx.get(x1->get(),x2->get())
 //       *B->get(vlagged,i,j,&xdir);
 //     xdir=2;
 //     Pol_B+=IMPACT_Heating::vosc_hat_By.get(x1->get(),x2->get())
 //       *B->get(vlagged,i,j,&xdir);
 //   }
 // xdir=3;
 // if (c->NB()>0)
 //   Pol_B+=IMPACT_Heating::vosc_hat_Bz.get(x1->get(),x2->get())
 //     *B->get(vlagged,i,j,&xdir);
 // Pol_B/=ei_cols;

 // //-----------------------------
 // int kminus=*k-1;
 // double Cee0_const = -1.0*(Pol_E+Pol_B)*c->v(k)*c->idv(k);
 // IMPACT_Vel_Sten Fk(0.0,0.0,0.0);
 // Fk.inc(0,-Cee0_flux_store::f2_IB->get(i,j,&kminus)->get(1));
 // Fk.inc(1,-Cee0_flux_store::f2_IB->get(i,j,&kminus)->get(2)
	// +Cee0_flux_store::f2_IB->get(i,j,k)->get(1));
 // Fk.inc(2,Cee0_flux_store::f2_IB->get(i,j,k)->get(2));
 // Fk*Cee0_const;
 // // Now add final term for B effect on IB
 // Fk.inc(1,-3.0*Pol_B*Cee0_flux_store::f2_IB->get(i,j,k)->get(1));
 // Fk.inc(2,-3.0*Pol_B*Cee0_flux_store::f2_IB->get(i,j,k)->get(2));
 
 // f0->insert(&Fk,vtemp,i,j,k);
 //-----------------------------

 //************************************************************

  /*
    Diffusion of f2 - obtained by substituting f2 into f3
    equation, then back into f2 equation
  */
  if (equation_switches::df3_by_dt_on)
    {
      IMPACT_Dim xtemp;
      GetOrthogonal(x1,&xtemp);
      IMPACT_Dim xdim(1),ydim(2);
      int iplus=*i+1,iminus=*i-1;
      double idx = 1.0/(c->dx(&iplus)+c->dx(&iminus));
      double diffusion_const=-c->v2(k)/(42.0*ei_cols);
      temp_sten = O->Laplacian(i,j);
      temp_sten = temp_sten*diffusion_const;
      f2->insert_BC(&temp_sten,vtemp,i,j,k,x1,x2);
      if (x1->get()<3&&x2->get()<3)
	{
	  if (x1->get()==x2->get())
	    {
	      temp_sten = (*O->d2dx2i(i,j,x1));
	      temp_sten = (temp_sten*diffusion_const)*0.8;
	      f2->insert_BC(&temp_sten,vtemp,i,j,k,x1,x1);
	      temp_sten = (*O->d2dx2i(i,j,&xtemp));
	      temp_sten = (temp_sten*diffusion_const)*(-0.4);
	      f2->insert_BC(&temp_sten,vtemp,i,j,k,&xtemp,&xtemp);
	      // Now dxdy derivative
	      temp_sten = (*O->ddy(j))*(0.4*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iplus,j,k,&xdim,&ydim);
	      temp_sten = (*O->ddy(j))*(-0.4*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iminus,j,k,&xdim,&ydim);	
	    }
	  else 
	    {
	      temp_sten = (*O->d2dx2i(i,j,&xdim))+(*O->d2dx2i(i,j,&ydim));
	      temp_sten = (temp_sten*diffusion_const)*0.6;
	      f2->insert_BC(&temp_sten,vtemp,i,j,k,x1,x2);
	      // Now dxdy derivative
	      temp_sten = (*O->ddy(j))*(0.6*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iplus,j,k,&xdim,&xdim);
	      f2->insert_BC(&temp_sten,vtemp,&iplus,j,k,&ydim,&ydim);
	      temp_sten = (*O->ddy(j))*(-0.6*diffusion_const*idx);
	      f2->insert_BC(&temp_sten,vtemp,&iminus,j,k,&xdim,&xdim);
	      f2->insert_BC(&temp_sten,vtemp,&iminus,j,k,&ydim,&ydim);
	    }
	}
      else
	{
	  if (x1->get()==3) temp_sten = (*O->d2dx2i(i,j,x2));
	  else temp_sten = (*O->d2dx2i(i,j,x1));
	  temp_sten = (temp_sten*diffusion_const)*0.6;
	  f2->insert_BC(&temp_sten,vtemp,i,j,k,x1,x2);
	  }
    }
  //**************************************
  if (IMPACTA_ions::ion_motion)
    { 
      // Ion motion terms
      // + C.grad f1
      for (IMPACT_Dim x4=1;x4<=c->Nf1();++x4)
	{
	  temp_sten = (*O->ddxi(i,j,&x4))*
	    Initial_Conditions::C_i[x4.get()-1].get(i,j);
	  f2->insert_BC(&temp_sten,vtemp,i,j,k,x1,x2);
	}
    }
}
