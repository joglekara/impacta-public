/*
**********************************************************
 E equation IMPACTA code for inner and outer cells

Version 1.0.0
1
AGRT
20/2/07
24/5/07 = major new reformation of equation - now curl B
is at n+1/2 not n+1
**********************************************************

*/
//code for inner cells of f0 equation - i.e. not at boundaries
void Eequation_inner(IMPACT_Config *c, IMPACT_ParVec *vlagged,IMPACT_ParVec *vtemp, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,int *i,int *j,IMPACT_Dim *x1)
{
  /*
**********************************************************
INSERT E Equation terms CODE HERE
**********************************************************
   */
  // En+1
  E->inc(vtemp,equation_switches::disp_j_on*c->idt(),i,j,x1);
  //E->inc(vtemp,1.0,i,j,x1);

 // This cyclicly works out orthoganol correct coordinates for x1 andx2
  double multiplier; 
  IMPACT_Dim x2,x3;

  IMPACT_stencil temp_sten;
  //**************************************************
  // New bit first insert c^2dt^2 grad of div E * alpha 
  int iplus=*i+1,iminus=*i-1;
 
  multiplier = equation_switches::dEbydtfix*c_L*c_L
    *equation_switches::Bimp_in_E_equation*c->dt();

  if (x1->get()<3)
    {
      GetOrthogonal(x1,&x2);
      
      //d/dx_i d/dx_jE_j
      temp_sten = (*O->ddy(j))*(multiplier*c->idx(i)*0.5);
      E->insert_yonly(&temp_sten,vtemp,&iplus,j,&x2);
      iplus=*i+1,iminus=*i-1;
      temp_sten = (*O->ddy(j))*(-multiplier*c->idx(i)*0.5);
      E->insert_yonly(&temp_sten,vtemp,&iminus,j,&x2);		  
      
      //-d2/dx2_j Ei
      temp_sten = (*O->d2dx2i(i,j,&x2))*(-multiplier);
      E->insert(&temp_sten,vtemp,i,j,x1);
    }
  else  // z component
    {
      // - del^2Ez
      temp_sten = O->Laplacian(i,j);
      temp_sten=temp_sten*(-multiplier);
      E->insert(&temp_sten,vtemp,i,j,x1);
    }
  //**************************************************
  
  //curl B:
  GetOrthogonal(x1,&x2,&x3);
  multiplier=-c_L*c_L*(equation_switches::Bimp_in_E_equation
		       -equation_switches::dEbydtfix);
  //c_L - speed of light.
  
  if (x3>3-c->NB())
    {
      temp_sten = (*O->ddxi(i,j,&x2))*(multiplier);
      B->insert(&temp_sten,vtemp,i,j,&x3);
    }
  if (x2>3-c->NB())
    {
      temp_sten = (*O->ddxi(i,j,&x3))*(-multiplier);
      B->insert(&temp_sten,vtemp,i,j,&x2);
    }
  //now insert j
  f1->insert(O->je(x1),vtemp, i,j,x1);
}
//code for outer cells of f0 equation - i.e.  boundaries
void Eequation_outer(IMPACT_Config *c, IMPACT_ParVec *vlagged,IMPACT_ParVec * vtemp, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,int *i,int *j,IMPACT_Dim *x1)
{
  /*
**********************************************************
INSERT E Equation terms CODE HERE - BC versions
**********************************************************
*/
   // En+1
  E->inc(vtemp,equation_switches::disp_j_on*c->idt(),i,j,x1);
  //E->inc(vtemp,1.0,i,j,x1);
    //curl B:
  double multiplier;
  // This cyclicly works out orthoganal correct coordinates for x1 andx2
  IMPACT_Dim x2,x3;
  IMPACT_stencil temp_sten;
  //**************************************************
  // New bit first insert c^2dt^2 grad of div E * alpha 
  int iplus=*i+1,iminus=*i-1;
   multiplier = equation_switches::dEbydtfix*c_L*c_L
     *equation_switches::Bimp_in_E_equation*c->dt();

  if (x1->get()<3)
    {  
      GetOrthogonal(x1,&x2);
      
      //d/dx_i d/dx_jE_j
      temp_sten =(*O->ddy(j))*(multiplier*c->idx(i)*0.5);
      E->insert_yonly_EBC(&temp_sten,vtemp,&iplus,j,&x2);
      iplus=*i+1,iminus=*i-1;
      temp_sten = (*O->ddy(j))*(-multiplier*c->idx(i)*0.5);
      E->insert_yonly_EBC(&temp_sten,vtemp,&iminus,j,&x2);
      
      //-d2/dx2_j Ei
      temp_sten = (*O->d2dx2i(i,j,&x2))*(-multiplier);
      E->insert_EBC(&temp_sten,vtemp,i,j,x1);
    }
  else // z component
    {
      // - del^2Ez
      temp_sten = O->Laplacian(i,j);
      temp_sten=temp_sten*(-multiplier);
      E->insert_EBC(&temp_sten,vtemp,i,j,x1);
      }
  //**************************************************

  //curl B:
  GetOrthogonal(x1,&x2,&x3);
  multiplier=-c_L*c_L*(equation_switches::Bimp_in_E_equation
		       -equation_switches::dEbydtfix);
 
  if (x3>3-c->NB())
    {
      temp_sten = (*O->ddxi(i,j,&x2))*(multiplier);
      B->insert_BBC(&temp_sten,vtemp,i,j,&x3);
    }
  if (x2>3-c->NB())
    {
      temp_sten = (*O->ddxi(i,j,&x3))*(-multiplier);
      B->insert_BBC(&temp_sten,vtemp,i,j,&x2);
    }
  //now insert j
  f1->insert(O->je(x1),vtemp, i,j,x1);
}
