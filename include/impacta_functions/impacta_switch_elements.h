/*
**********************************************************
 Switch elements IMPACTA code
- This adjusts the RHS vector of the matrix equation
- for d/dt=0, or explicit bits

1
AGRT

1/25/2010 AGRT 
Having problems with checkboard pattern arising due to 
checkboard like pattern (decoupling of adjacent gridcells
due to center differencing?) 

To fix it I will try adding some numerical diffusion by
introducing a Lax-like averaging over adjacent grid-cells
in f0 and f1.

2/8/2010
This works but I can't diffuse f0 as density structures then melt 
away. Instead, f1 only works but a similar problem arises in the
B and E fields. I am going to diffuse away electric field 
perturbations at the grid level using the same method.

march 2011 - added contribution from ionization collisional stopping
**********************************************************

*/

#ifndef INC_IMPACTA_SWITCH_ELEMENTS_H
#define INC_IMPACTA_SWITCH_ELEMENTS_H



void IMPACT_Switch_Elements(IMPACT_Config *c, IMPACT_MPI_Config *MPIc, 
			    IMPACT_ParVec *v, IMPACT_StenOps *O,
			    IMPACT_Var* f0,IMPACT_Var* f1,
			    IMPACT_Var* f2, IMPACT_Var* f3,IMPACT_Var* E,
			    IMPACT_Var* B, int n)
{
  IMPACT_ParVec *vtemp; // temporary store of row elements
  vtemp = new IMPACT_ParVec(MPIc);
  vtemp->Duplicate(v);
  double jint_temp[3]={0.0,0.0,0.0};
  double f0temp; // temporary store of f0
  IMPACT_Dim x3d(3);
  int istart=MPIc->istart();
  int iend=MPIc->iend();
  // New -> correct istart and iend for fixed boundaries:
  if (istart==1) {
    istart=(1+IMPACT_Boundaries::fixed_any[0]);
  }
  if (iend==c->Nx()) {
    iend=(c->Nx()-IMPACT_Boundaries::fixed_any[1]);
  }
  
  double BxC[3];
  double correctE[3];
  double dv3f2bydv=0.0;

  IMPACT_Dim x2f,x3f;
  int one =1;
  int nadd=n;
  if (nadd<1) nadd=1;

  MPI_Barrier(MPI_COMM_WORLD);

  /*
    AGRT 25/1/2010
    Smoothing step - somewhat equivalent to Lax's method
   */
  // *****************************************************
  double tempvalue=0.0;
  IMPACT_stencil temp_sten;
  IMPACT_stencil dummy_sten; // This is only used as a placeholder
  int iplus,iminus,jplus,jminus;

  for (int i=istart;i<=iend;++i)
    {  
      for (int j=(1+IMPACT_Boundaries::fixed_any[2]); j<=(c->Ny()-IMPACT_Boundaries::fixed_any[3]); ++j)
	{
	  
	  //*******************************
	  // f1
	  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
	    {
	      // set elements of stencil
	      temp_sten=LAX;
	      //reset plus and minus elements
	      iplus = i+1; 
	      iminus = i-1;
	      jplus = j+1; 
	      jminus = j-1;
	  
	      // Apply boundary condition for f1
	      IMPACT_f1_E_bound(&temp_sten,c->Nx(),c->Ny(),&iplus,&iminus,&jplus,&jminus,&x1);
	      for (int k=1;k<=c->Nv();++k)
		{    
		  tempvalue = f1->get(v,&i,&j,&k,&x1)*temp_sten(0);
		  tempvalue += f1->get(v,&iplus,&j,&k,&x1)*temp_sten(1);
		  tempvalue += f1->get(v,&iminus,&j,&k,&x1)*temp_sten(2);
		  tempvalue += f1->get(v,&i,&jplus,&k,&x1)*temp_sten(3);
		  tempvalue += f1->get(v,&i,&jminus,&k,&x1)*temp_sten(4);
		  
		  f1->set(v,tempvalue,&i,&j,&k,&x1); 
		}//end of k loop
	      //I can now use the same stenicl fo rthe electric field.
	      tempvalue = E->get(v,&i,&j,&x1)*temp_sten(0);
	      tempvalue += E->get(v,&iplus,&j,&x1)*temp_sten(1);
	      tempvalue += E->get(v,&iminus,&j,&x1)*temp_sten(2);
	      tempvalue += E->get(v,&i,&jplus,&x1)*temp_sten(3);
	      tempvalue += E->get(v,&i,&jminus,&x1)*temp_sten(4);

	      E->set(v,tempvalue,&i,&j,&x1);
	    }// end of x1 loop*/
	}// end of j loop
    }// end of i loop
  
  // *****************************************************
  
  SwapGhosts(c,MPIc, v);
  double df0dt_neutrals=0.0; // df/dt due to neutral collisions

  // New 2011 - now adding collisional stopping due to ionization
  // Non-relativistic Bethe formula 1.e+05 is 1/2mc^2 / 10 eV
  
  //Here, ni(1-Z) is the unionized density
  double neutralcollisions=0.0;
 
  //**************************************************************************************************
  for (int i=istart;i<=iend;++i)
    for (int j=(1+IMPACT_Boundaries::fixed_any[2]);j<=(c->Ny()-IMPACT_Boundaries::fixed_any[3]);++j)
      {
	// New 2011 - now adding collisional stopping due to ionization
	// Non-relativistic Bethe formula 1.e+05 is 1/2mc^2 / 10 eV
	neutralcollisions=Initial_Conditions::ni.get(&i,&j)*(1.0-Initial_Conditions::Z.get(&i,&j))*ionize_coll_log*oneover_atomic_Z;
        
	// NEW for correction to E field 2011 - see thomas JCP 2011
        // E -> E+CxB-DUbyDt

        BxC[0]=0.0;  BxC[1]=0.0; BxC[2]=0.0;
	if (IMPACTA_ions::ion_motion) {
        for (IMPACT_Dim x1=1;x1<=3;++x1)
          {
            // this is +DUbyDt
            correctE[x1.get()-1] = -1.0*(IMPACTA_DUbyDt(c,MPIc,v,O,f0,f1,f2,E,B,i,j,x1,nadd));
            
            GetOrthogonal(&x1,&x2f,&x3f);
            if (x3f>3-c->NB()) {
              BxC[x1.get()-1] += Initial_Conditions::C_i[x2f.get()-1].get(&i,&j)*B->get(v,&i,&j,&x3f);
            }
            if (x2f>3-c->NB()) {
              BxC[x1.get()-1] -= Initial_Conditions::C_i[x3f.get()-1].get(&i,&j)*B->get(v,&i,&j,&x2f);
            }
            // Notice correction is already -ve (to put on RHS)
            correctE[x1.get()-1] += BxC[x1.get()-1];
          }
        }
	//for semi-implicit step
	jint_temp[0]=jint_temp[1]=jint_temp[2]=0.0;

	for (int k=1;k<=c->Nv();++k)
	  {
	    // std::cout<<df0dt_neutrals<<','<<f0temp*c->idt()<<','<<c->v2(&k)<<','<<f0->ddv(v,c,&i,&j,&k)<<'\n';
	    if (equation_switches::f0_equation_on)
	      {
	
		f0temp=f0->get(v,&i,&j,&k);
		if (IMPACTA_ions::ionization_on) {
		  df0dt_neutrals =0.0;
		  if (k>1){
		    df0dt_neutrals = f0->ddv_BC(v,c,&i,&j,&k)*neutralcollisions/c->v2(&k);	 
		    df0dt_neutrals = equation_switches::df0_by_dt_on*f0temp*c->idt()+df0dt_neutrals;
		    if (df0dt_neutrals<0.0) {
		      df0dt_neutrals=0.0;
		      f0->inc(vtemp,f0temp*c->v2(&k)*c->dv(&k)/c->v2(&one)/c->dv(&one)*c->idt(),&i,&j,&one); // this adds all the charge to the first bin
		    }
		  } else {
		    df0dt_neutrals = equation_switches::df0_by_dt_on*f0temp*c->idt();
		  }
		} else {
		  df0dt_neutrals = equation_switches::df0_by_dt_on*f0temp*c->idt();
		}
		f0->set(vtemp,df0dt_neutrals,&i,&j,&k);
		
                // ION MOTION
                //------------------------------------------
		// Now do BxC/3v^2d/dv(v^2f1)
		
		if (IMPACTA_ions::ion_motion) {
                  for (IMPACT_Dim x1=1;x1<c->Nf1();++x1) {
                    f0->inc(vtemp,oneover3/c->v2(&k)*f1->ddv_v2_BC(v,c,&i,&j,&k,&x1)*correctE[x1.get()-1],&i,&j,&k);
                  }
                }
		//------------------------------------------

		
	      }
	    //turn off electron inertia....
	    for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
	      {
		jint_temp[x1.get()-1]+=f1->get(v,&i,&j,&k,&x1)
		  *O->je(&x1)->get(k);
		if (equation_switches::f1_equation_on)
		  {
		    f1->set(vtemp,equation_switches::e_inert_on
			    *f1->get(v,&i,&j,&k,&x1)*c->idt(),&i,&j,&k,&x1);
		    
                    //------------------------------------------
                    //ION MOTION
                    // Cxwdf0/dv	
		    if (IMPACTA_ions::ion_motion) {
		      
		      f1->inc(vtemp,f0->ddv_BC(v,c,&i,&j,&k)*correctE[x1.get()-1] ,&i,&j,&k,&x1);
		    
		      if (c->Nf2()>0) {
			for (IMPACT_Dim x2=1;x2<=c->N3f2();++x2) { //N3f2 is whether z is included or not
			  dv3f2bydv=f2->ddv_v3_BC(v,c,&i,&j,&k,&x1,&x2);
			  f1->inc(vtemp,0.4/c->v3(&k)*dv3f2bydv*correctE[x2.get()-1],&i,&j,&k,&x1);
			}
		      }
		    }
		  }
	      }
	    
	    //or viscosity!
	    if (c->Nf2()>0) {
              double multiplier=1.0;
	      for (IMPACT_Dim x1=1;x1<=c->N3f2();++x1)
		for (IMPACT_Dim x2=x1.get();x2<=c->N3f2();++x2)
		  if (x1.get()+x2.get()<6)
		    {
		      f2->set(vtemp,equation_switches::e_viscosity_on
			      *f2->get(v,&i,&j,&k,&x1,&x2)*c->idt(),
			      &i,&j,&k,&x1,&x2);

                      //------------------------------------------
                      // ION MOTION
                      // v(A I - 1/3 I A) . d(f1/v)/dv
			if (IMPACTA_ions::ion_motion) {
                      //v/2aid/dvfj/v + v/2ajd/dvfi/v
                      multiplier=c->v(&k)*0.5;
                      f2->inc(vtemp,-multiplier*f1->ddv1overv_BC(v, c,&i,&j,&k,&x1)*correctE[x2.get()-1],&i,&j,&k,&x1,&x2);
                      f2->inc(vtemp,-multiplier*f1->ddv1overv_BC(v, c,&i,&j,&k,&x2)*correctE[x1.get()-1],&i,&j,&k,&x1,&x2);
                      //v/2 .2/3 .deltaij a.df1/dv
                      multiplier=-oneover3*c->v(&k)*Kronecker_Delta(&x1,&x2);
                      if (c->NE()>0)
                        for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
                          {
                            f2->inc(vtemp,-multiplier*f1->ddv1overv_BC(v,c,&i,&j,&k,&x3)*correctE[x3.get()-1],&i,&j,&k,&x1,&x2);
                          }
			}
                      //------------------------------------------
		    }
            }
	    
	  } //end of k loop
	//These next two are to make the j and B terms in the E
	//equation at n+1/2 or n+1
	for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
	  {
	    if (!n) {
	      jint_temp[x1.get()-1] = 
		B->Curl(v,c,O,&i,&j,&x1)*c_L*c_L;
	    }
	    jint_temp[x1.get()-1]*=
	      (1.0-equation_switches::jimp_in_E_equation);
	      /*
		eliminated Jan 2011 
	      /equation_switches::jimp_in_E_equation;

	      - not sure why it was there.
	      */
	  }
	//********************************
	if (equation_switches::E_equation_on)
	  for (IMPACT_Dim x1=1;x1<=c->NE();++x1) //v2
	    E->set(vtemp,E->get(v,&i,&j,&x1)*equation_switches::disp_j_on*c->idt()
		   +B->Curl(v,c,O,&i,&j,&x1)*c_L*c_L*(1.0-equation_switches::Bimp_in_E_equation+equation_switches::dEbydtfix)
		   -jint_temp[x1.get()-1]
		   ,&i,&j,&x1);

	//********************************
	if (equation_switches::B_equation_on)
	  for (IMPACT_Dim x1=3;x1>3-c->NB();--x1)
	    B->set(vtemp,B->get(v,&i,&j,&x1)*c->idt()
		   -E->Curl(v,c,O,&i,&j,&x1)
		   *(1.0-equation_switches::Eimp_in_B_equation)
		   ,&i,&j,&x1);
      }
  SwapGhosts(c,MPIc, vtemp);
  v->Duplicate(vtemp);

  delete vtemp;
}



#endif /* INC_IMPACTA_SWITCH_ELEMENTS_H  */
