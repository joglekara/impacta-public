/*
**********************************************************
IMPACTA Sparse counting functions

Need one for each equation - sorry!
Version 1.1
AGRT

7/3/07

8/3/07 - chk_ij is periodic as it doesn't matter in this context
even if the boundary conditions are not periodic - if they are not
periodic the values will be zero anyway.
- now moved to Operators.h

13/3/07 - slight changes due to parallelized version....
          i boundary check now eliminated.
24/2/08 - now ion motion terms included
**********************************************************
*/


//This function counts non zeros in f0 equation 

inline  void f0equation_count(IMPACT_ParVec * v, IMPACT_ParSpa *S,IMPACT_Config *c, IMPACT_Var * f0,IMPACT_Var * E,IMPACT_Var * B,IMPACT_Var * f1,int*i,int*j,int*k)
{
  int rownumber=f0->getrow(i,j,k)+1; // plus one so Resize works
  S->setrow(rownumber,0);//numinrows[rownumber-1];
  //First put in 
  for (int kk=1;kk<=c->Nv();++kk)
    S->incrow(rownumber,v->chkzero(f0->getrow(i,j,&kk)));
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(c,&jplus,&jminus);

  //Now f0 in 5 point stencil due to ion motion
  //  S->incrow(rownumber,v->chkzero(f0->getrow(i,j,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(&iplus,j,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(&iminus,j,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(i,&jplus,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(i,&jminus,k)));

  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    {
      //Now f1 in 5 point stencil
      S->incrow(rownumber,v->chkzero(f1->getrow(i,j,k,&x1)));
      S->incrow(rownumber,v->chkzero(f1->getrow(&iplus,j,k,&x1)));
      S->incrow(rownumber,v->chkzero(f1->getrow(&iminus,j,k,&x1)));
      S->incrow(rownumber,v->chkzero(f1->getrow(i,&jplus,k,&x1)));
      S->incrow(rownumber,v->chkzero(f1->getrow(i,&jminus,k,&x1)));
    }
  //Now E
  if (c->NE()>0)
    for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
      S->incrow(rownumber,v->chkzero(E->getrow(i,j,&x1)));
  for (IMPACT_Dim x5=3;x5>3-c->NB();--x5)
    S->incrow(rownumber,v->chkzero(B->getrow(i,j,&x5)));

}
inline  void f1equation_count(IMPACT_ParVec * v, IMPACT_ParSpa *S,IMPACT_Config *c, IMPACT_Var * f0,IMPACT_Var * E,IMPACT_Var *B,IMPACT_Var * f1,IMPACT_Var * f2,int*i,int*j,int*k,IMPACT_Dim * x1)
{
  int rownumber=f1->getrow(i,j,k,x1)+1;
  S->setrow(rownumber,0);//numinrows[rownumber-1];
  //First count all f1 values for k in i,j
  for (int kk=1;kk<=c->Nv();++kk)
    S->incrow(rownumber,v->chkzero(f1->getrow(i,j,&kk,x1)));
  //And orthogonal bits!
  IMPACT_Dim x2,x3;
  GetOrthogonal(x1,&x2,&x3);
  S->incrow(rownumber,v->chkzero(f1->getrow(i,j,k,&x2)));
  S->incrow(rownumber,v->chkzero(f1->getrow(i,j,k,&x3)));
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(c,&jplus,&jminus);
  
  //Now f0 and f2 in 5 point stencil
  S->incrow(rownumber,v->chkzero(f0->getrow(i,j,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(&iplus,j,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(&iminus,j,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(i,&jplus,k)));
  S->incrow(rownumber,v->chkzero(f0->getrow(i,&jminus,k)));
  for (int pos=0;pos<c->Nf2();++pos)
    {
    S->incrow(rownumber,v->chkzero(f2->getrow(i,j,k,&pos)));
  S->incrow(rownumber,v->chkzero(f2->getrow(&iplus,j,k,&pos)));
  S->incrow(rownumber,v->chkzero(f2->getrow(&iminus,j,k,&pos)));
  S->incrow(rownumber,v->chkzero(f2->getrow(i,&jplus,k,&pos)));
  S->incrow(rownumber,v->chkzero(f2->getrow(i,&jminus,k,&pos)));
    }
  
  //Now E & B
  if (c->NE()>0)
    for (IMPACT_Dim x4=1;x4<=c->NE();++x4)
      S->incrow(rownumber,v->chkzero(E->getrow(i,j,&x4)));
  for (IMPACT_Dim x5=3;x5>3-c->NB();--x5)
    S->incrow(rownumber,v->chkzero(B->getrow(i,j,&x5)));
}
inline  void f2equation_count(IMPACT_ParVec * v, IMPACT_ParSpa *S,IMPACT_Config *c, IMPACT_Var * f0, IMPACT_Var * f1,IMPACT_Var * E,IMPACT_Var *B,IMPACT_Var * f2,IMPACT_Var * f3,int*i,int*j,int*k,IMPACT_Dim * x1,IMPACT_Dim * x2)
{

  int rownumber=f2->getrow(i,j,k,x1,x2)+1;
  S->setrow(rownumber,0);//numinrows[rownumber-1];

  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(c,&jplus,&jminus);
  
  //Now f1 and f2 in 5 point stencil
  for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
    {
      S->incrow(rownumber,v->chkzero(f1->getrow(i,j,k,&x3)));
      S->incrow(rownumber,v->chkzero(f1->getrow(&iplus,j,k,&x3)));
      S->incrow(rownumber,v->chkzero(f1->getrow(&iminus,j,k,&x3)));
      S->incrow(rownumber,v->chkzero(f1->getrow(i,&jplus,k,&x3)));
      S->incrow(rownumber,v->chkzero(f1->getrow(i,&jminus,k,&x3)));
    }
  for (int pos=0;pos<c->Nf2();++pos)
    {
      S->incrow(rownumber,v->chkzero(f2->getrow(i,j,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(&iplus,j,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(&iminus,j,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(i,&jplus,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(i,&jminus,k,&pos)));
      // New additions
      S->incrow(rownumber,v->chkzero(f2->getrow(&iplus,&jplus,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(&iminus,&jplus,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(&iplus,&jminus,k,&pos)));
      S->incrow(rownumber,v->chkzero(f2->getrow(&iminus,&jminus,k,&pos)));
    }
  
  //Now E & B
  if (c->NE()>0)
    for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
      S->incrow(rownumber,v->chkzero(E->getrow(i,j,&x3)));
  for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
    S->incrow(rownumber,v->chkzero(B->getrow(i,j,&x3)));
  // last do f0 put in 
  for (int kk=1;kk<=c->Nv();++kk)
    S->incrow(rownumber,v->chkzero(f0->getrow(i,j,&kk)));

}
inline  void Eequation_count(IMPACT_ParVec * v, IMPACT_ParSpa *S,IMPACT_Config *c,IMPACT_Var * E,IMPACT_Var *B,IMPACT_Var * f1,int*i,int*j,IMPACT_Dim * x1)
{
  int rownumber=E->getrow(i,j,x1)+1;
  S->setrow(rownumber,0);//numinrows[rownumber-1];
  //First count all f1 values for k in i,j
  for (int kk=1;kk<=c->Nv();++kk)
    S->incrow(rownumber,v->chkzero(f1->getrow(i,j,&kk,x1)));
  
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(c,&jplus,&jminus);

  //Now E in 5 point stencil 
//Also have to include elements within 9 point stencil 
 if (c->NE()>0)
    for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
      {
	S->incrow(rownumber,v->chkzero(E->getrow(i,j,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iplus,j,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iminus,j,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(i,&jplus,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(i,&jminus,&x2)));
	// New additions
	S->incrow(rownumber,v->chkzero(E->getrow(&iplus,&jplus,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iminus,&jplus,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iplus,&jminus,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iminus,&jminus,&x2)));

      }
//Now B in 5 point stencil
for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
  {
    S->incrow(rownumber,v->chkzero(B->getrow(i,j,&x2)));
    S->incrow(rownumber,v->chkzero(B->getrow(&iplus,j,&x2)));
    S->incrow(rownumber,v->chkzero(B->getrow(&iminus,j,&x2)));
    S->incrow(rownumber,v->chkzero(B->getrow(i,&jplus,&x2)));
    S->incrow(rownumber,v->chkzero(B->getrow(i,&jminus,&x2)));
  }
}
inline  void Bequation_count(IMPACT_ParVec * v, IMPACT_ParSpa *S,IMPACT_Config *c, IMPACT_Var * E,IMPACT_Var *B,int*i,int*j,IMPACT_Dim * x1)
{
  int rownumber=B->getrow(i,j,x1)+1;
  S->setrow(rownumber,0);//numinrows[rownumber-1];
  
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(c,&jplus,&jminus);

  //Now E & B
  S->incrow(rownumber,v->chkzero(B->getrow(i,j,x1)));
  //Now E in 5 point stencil
  if (c->NE()>0)
    for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
      {
	S->incrow(rownumber,v->chkzero(E->getrow(i,j,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iplus,j,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(&iminus,j,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(i,&jplus,&x2)));
	S->incrow(rownumber,v->chkzero(E->getrow(i,&jminus,&x2)));
      }
}
