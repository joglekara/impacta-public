/*
 IMPACT EXPlicit

Makes linear solution for f0, f1 E and B 
26/10/07

*/

int IMPACT_Linear_Soln(IMPACT_ParVec *v,IMPACT_Config *c, IMPACT_MPI_Config *M,
		       IMPACT_StenOps *O,IMPACT_Var *f0,IMPACT_Var *f1,
		       IMPACT_Var *f2,IMPACT_Var *f3,IMPACT_Var *E,
		       IMPACT_Var *B)
{
  // E = -1.25 grad T , f1=v/nu(4-v^2/T)*f0*gradT/T
  // f0 =f0 +heating
  
  double gradTx=0.0,gradTy=0.0,T=0.0,alpha=-1.25,ei_cols;
  IMPACT_Dim x1(1),x2(2);
  int iplus,iminus,jplus,jminus;
  IMPACT_stencil stencil;

  // First heat f0

  IMPACT_Maxwellianf0evolve(v,c,M,O,f0,f1,E); // get f0 from heating
      
  //Now get f1 and E
  for (int i=M->istart();i<=M->iend();++i)
    for (int j=1;j<=c->Ny();++j)
	{
	  iplus =i+1;
	  iminus=i-1;
	  jplus =j+1;
	  jminus=j-1;
	  T = Local_Te(v,f0,c,&i,&j);
	  
	  stencil = *O->ddx(&i);
	  IMPACT_f0_bound(&stencil,c->Nx(),c->Ny(),&iplus,&iminus,
			  &jplus,&jminus);
	  gradTx=Local_Te(v,f0,c,&iplus,&j)*stencil(1)
	    +Local_Te(v,f0,c,&iminus,&j)*stencil(2);

	  stencil = *O->ddy(&j);
	  IMPACT_f0_bound(&stencil,c->Nx(),c->Ny(),&iplus,&iminus,
			  &jplus,&jminus);
	  gradTy=Local_Te(v,f0,c,&i,&jplus)*stencil(3)
	    +Local_Te(v,f0,c,&i,&jminus)*stencil(4);
	  
	  //Now get linear solns for f1 and E
	  // E = -1.25*gradT
	  
	  E->set(v,alpha*gradTx,&i,&j,&x1);
	  E->set(v,alpha*gradTy,&i,&j,&x2);
	  ei_cols=equation_switches::Cei_on*Initial_Conditions::Z.get(&i,&j)
	    *Initial_Conditions::Z.get(&i,&j)*Initial_Conditions::ni.get(&i,&j);

	    // f1=v^4/nu_ei*(4-v^2/vth^2)*f0/l_T
	  for (int k=1;k<=c->Nv();++k)
	    {
	      f1->set(v,c->v2(&k)*c->v2(&k)/ei_cols*(4.0-c->v2(&k)/T)*
		      f0->get(v,&i,&j,&k)*gradTx/T,&i,&j,&k,&x1);
	      f1->set(v,c->v2(&k)*c->v2(&k)/ei_cols*(4.0-c->v2(&k)/T)*
		      f0->get(v,&i,&j,&k)*gradTy/T,&i,&j,&k,&x2);
	    }
	}
  return 0;
}
