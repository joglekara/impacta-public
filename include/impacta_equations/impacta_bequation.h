/*
**********************************************************
 B equation IMPACT code for inner and outer cells

Version 1.0.0
1
AGRT
20/2/07

**********************************************************

*/
//code for inner cells of f0 equation - i.e. not at boundaries
void Bequation_inner(IMPACT_Config *c, IMPACT_ParVec *vlagged,IMPACT_ParVec *vtemp, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,int *i,int *j,IMPACT_Dim *x1)
{
  /*
**********************************************************
INSERT B Equation terms CODE HERE
**********************************************************
   */
  // Bn+1
  B->inc(vtemp,c->idt(),i,j,x1);
   // - curl E:
  double multiplier;
  // This cyclicly works out orthoganol correct coordinates for x1 andx2
  IMPACT_Dim x2,x3;
  GetOrthogonal(x1,&x2,&x3);

  multiplier=equation_switches::Eimp_in_B_equation; 
  IMPACT_stencil temp_sten;
  if (c->NE()>0)
    {
    if (x3<=c->NE())
      {
	temp_sten = *O->ddxi(i,j,&x2)*(multiplier);
	E->insert(&temp_sten,vtemp,i,j,&x3);
      }
    if (x2<=c->NE())
      {
	temp_sten = *O->ddxi(i,j,&x3)*(-multiplier);
	E->insert(&temp_sten,vtemp,i,j,&x2);
      }
    }
}
//code for outer cells of f0 equation - i.e.  boundaries
void Bequation_outer(IMPACT_Config *c, IMPACT_ParVec *vlagged,IMPACT_ParVec * vtemp, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,int *i,int *j,IMPACT_Dim *x1)
{
  /*
**********************************************************
INSERT B Equation terms CODE HERE - BC versions
**********************************************************
*/
// Bn+1
  B->inc(vtemp,c->idt(),i,j,x1);
 // - curl E:
  double multiplier;
   // This cyclicly works out orthoganol correct coordinates for x1 andx2
IMPACT_Dim x2,x3;
  GetOrthogonal(x1,&x2,&x3);

  multiplier=equation_switches::Eimp_in_B_equation; 
  IMPACT_stencil temp_sten;
  if (c->NE()>0)
    {
    if (x3<=c->NE())
      {
	temp_sten = *O->ddxi(i,j,&x2)*(multiplier);
	E->insert_EBC(&temp_sten,vtemp,i,j,&x3);
      }
    if (x2<=c->NE())
      {
	temp_sten = *O->ddxi(i,j,&x3)*(-multiplier);
	E->insert_EBC(&temp_sten,vtemp,i,j,&x2);
      }
    }
}
