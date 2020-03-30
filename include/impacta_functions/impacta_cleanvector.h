/*
**********************************************************
Code for setting appropriate elements in v to zero

NOTE -WHEN f2 added need to modify f1 clean equation!!!
Version IMPACTA 1.1

AGRT
6/3/07

8/3/07   - Updated so that there exists a zero checking version 
for debugging - it checks that v really has been zeroed.

13/3/07 - Boundary code updated in line for the parallel version.
11/2/08 - updated to include IB heating in f2 equation
24/2/08 - updated to include ion motion in f0 equation
**********************************************************

*/
// EFFICIENT ZEROING OF MATRIX ELEMENTS

inline void IMPACT_Var::clean(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k) // for f0
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  for (int kk=*k-1;kk<=*k+1;++kk)
    {
      // int kk=*k;
      set(v,0.0,&iplus,j,&kk);
      set(v,0.0,&iminus,j,&kk);
      set(v,0.0,i,&jplus,&kk);
      set(v,0.0,i,&jminus,&kk);
    }
 for (int kk=1;kk<=config1->Nv();++kk)
   set(v,0.0,i,j,&kk);
}
inline void IMPACT_Var::clean(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k, IMPACT_Dim * x1) //for f1
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  set(v,0.0,&iplus,j,k,x1);
  set(v,0.0,&iminus,j,k,x1);
  set(v,0.0,i,&jplus,k,x1);
  set(v,0.0,i,&jminus,k,x1);
 for (int kk=1;kk<=config1->Nv();++kk)
   set(v,0.0,i,j,&kk,x1);
}
inline void IMPACT_Var::clean(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, IMPACT_Dim* x1) // for E/B
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  set(v,0.0,i,j,x1);
  set(v,0.0,&iplus,j,x1);
  set(v,0.0,&iminus,j,x1);
  set(v,0.0,i,&jplus,x1);
  set(v,0.0,i,&jminus,x1);

}
inline void IMPACT_Var::clean(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k, int *index) //for f2
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  set(v,0.0,&iplus,j,k,index);
  set(v,0.0,&iminus,j,k,index);
  set(v,0.0,i,&jplus,k,index);
  set(v,0.0,i,&jminus,k,index);
  set(v,0.0,i,j,k,index);
}
inline void IMPACT_Var::clean(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k, IMPACT_Dim *x1, IMPACT_Dim *x2) //for f2
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  set(v,0.0,&iplus,j,k,x1,x2);
  set(v,0.0,&iminus,j,k,x1,x2);
  set(v,0.0,i,&jplus,k,x1,x2);
  set(v,0.0,i,&jminus,k,x1,x2);
  set(v,0.0,i,j,k,x1,x2);
}
inline void IMPACT_Var::cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k) // for f0
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1,kplus=*k+1,kminus=*k-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(config1,&jplus,&jminus);
  if (kminus<1) ++kminus;
  if (kplus>config1->Nv()) --kplus;
  for (int kk=kminus;kk<=kplus;++kk)
    {
      //int kk=*k;
      set(v,0.0,&iplus,j,&kk);
      set(v,0.0,&iminus,j,&kk);
      set(v,0.0,i,&jplus,&kk);
      set(v,0.0,i,&jminus,&kk);
    }
 for (int kk=1;kk<=config1->Nv();++kk)
   set(v,0.0,i,j,&kk);
}
inline void IMPACT_Var::cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k, IMPACT_Dim * x1) //for f1
{
   int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(config1,&jplus,&jminus);
  set(v,0.0,&iplus,j,k,x1);
  set(v,0.0,&iminus,j,k,x1);
  set(v,0.0,i,&jplus,k,x1);
  set(v,0.0,i,&jminus,k,x1);
 for (int kk=1;kk<=config1->Nv();++kk)
   set(v,0.0,i,j,&kk,x1);
}
inline void IMPACT_Var::cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, IMPACT_Dim* x1) // for E/B
{
 int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
 //chk_ij(c,&iplus,&jplus);
  //chk_ij(c,&iminus,&jminus);
  chk_j(config1,&jplus,&jminus);
  set(v,0.0,i,j,x1);
  set(v,0.0,&iplus,j,x1);
  set(v,0.0,&iminus,j,x1);
  set(v,0.0,i,&jplus,x1);
  set(v,0.0,i,&jminus,x1);
  set(v,0.0,&iplus,&jplus,x1);
  set(v,0.0,&iminus,&jplus,x1);
  set(v,0.0,&iplus,&jminus,x1);
  set(v,0.0,&iminus,&jminus,x1);
}

inline void IMPACT_Var::cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1, 
				   IMPACT_StenOps * O,int *i, int *j, int *k, 
				   int *index) //for f2
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  chk_j(config1,&jplus,&jminus);
  set(v,0.0,&iplus,j,k,index);
  set(v,0.0,&iminus,j,k,index);
  set(v,0.0,i,&jplus,k,index);
  set(v,0.0,i,&jminus,k,index);
  set(v,0.0,i,j,k,index);
  // New bits for Weibel damping
  set(v,0.0,&iplus,&jplus,k,index);
  set(v,0.0,&iminus,&jplus,k,index);
  set(v,0.0,&iplus,&jminus,k,index);
  set(v,0.0,&iminus,&jminus,k,index);
}
inline void IMPACT_Var::cleanouter(IMPACT_ParVec * v,IMPACT_Config *config1, IMPACT_StenOps * O,int *i, int *j, int *k, IMPACT_Dim *x1, IMPACT_Dim *x2) //for f2
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
 chk_j(config1,&jplus,&jminus);
  set(v,0.0,&iplus,j,k,x1,x2);
  set(v,0.0,&iminus,j,k,x1,x2);
  set(v,0.0,i,&jplus,k,x1,x2);
  set(v,0.0,i,&jminus,k,x1,x2);
  set(v,0.0,i,j,k,x1,x2);
  // New bits for Weibel damping
  set(v,0.0,&iplus,&jplus,k,x1,x2);
  set(v,0.0,&iminus,&jplus,k,x1,x2);
  set(v,0.0,&iplus,&jplus,k,x1,x2);
  set(v,0.0,&iminus,&jminus,k,x1,x2);
}
inline void vzerochk(IMPACT_ParVec *v) //checks v is actually zero throughout.
{
  double check=fabs(v->VecSum());
  if (check>0.0)
    {
      std::cout<<"IMPACT: ERROR - Vector not correctly zeroed."<<std::endl;
      exit(0);
    }
}
namespace with_zero_vector_checking
{
  int withzerovectorchecking=1;
inline void f0equation_clean_outer(IMPACT_ParVec * vtemp, IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var * B,IMPACT_Var* f1,int *i,int *j,int *k)
{
  f0->cleanouter(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
      E->cleanouter(vtemp,c,O,i,j,&x1);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->cleanouter(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    f1->cleanouter(vtemp,c,O,i,j,k,&x1);
  vzerochk(vtemp);
}
inline void f0equation_clean_inner(IMPACT_ParVec * vtemp, IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var * B,IMPACT_Var* f1,int *i,int *j,int *k)
{
  f0->clean(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
      E->clean(vtemp,c,O,i,j,&x1);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->clean(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    f1->clean(vtemp,c,O,i,j,k,&x1);
  vzerochk(vtemp);
}
inline void f1equation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,IMPACT_Var* f2,int* i,int *j,int *k,IMPACT_Dim *x1)
{
  f0->cleanouter(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
      E->cleanouter(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x2=1;x2<=c->Nf1();++x2)
    f1->cleanouter(vtemp,c,O,i,j,k,&x2);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->cleanouter(vtemp,c,O,i,j,&x2);
  if (c->Nf2()>0)
    for (int pos=0;pos<c->Nf2();++pos)
      f2->cleanouter(vtemp,c,O,i,j,k,&pos);
 vzerochk(vtemp);
}
inline  void f1equation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,IMPACT_Var* f2,int* i,int *j,int *k,IMPACT_Dim *x1)
 {
   f0->clean(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
      E->clean(vtemp,c,O,i,j,&x2);
   for (IMPACT_Dim x2=1;x2<=c->Nf1();++x2)
     f1->clean(vtemp,c,O,i,j,k,&x2);
   for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
     B->clean(vtemp,c,O,i,j,&x2);
if (c->Nf2()>0)
  for (int pos=0;pos<c->Nf2();++pos)
    f2->clean(vtemp,c,O,i,j,k,&pos);
   vzerochk(vtemp);
 }
 inline void f2equation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f0,IMPACT_Var* f1,IMPACT_Var* f2,IMPACT_Var* f3,int* i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
    f0->cleanouter(vtemp,c,O,i,j,k);
  if (c->Nf2()>0)
    for (int pos=0;pos<c->Nf2();++pos)
      f2->cleanouter(vtemp,c,O,i,j,k,&pos);
  if (c->Nf1()>0)
    for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
      f1->cleanouter(vtemp,c,O,i,j,k,&x3);
  for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
    E->cleanouter(vtemp,c,O,i,j,&x3);
  for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
    B->cleanouter(vtemp,c,O,i,j,&x3);
  vzerochk(vtemp);
}
inline void f2equation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f0,IMPACT_Var* f1,IMPACT_Var* f2,IMPACT_Var* f3,int* i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  f0->clean(vtemp,c,O,i,j,k);
  if (c->Nf2()>0)
 for (int pos=0;pos<c->Nf2();++pos)
    f2->cleanouter(vtemp,c,O,i,j,k,&pos);
if (c->Nf1()>0)
  for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
    f1->clean(vtemp,c,O,i,j,k,&x3);
  for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
    E->clean(vtemp,c,O,i,j,&x3);
  for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
    B->clean(vtemp,c,O,i,j,&x3);
  vzerochk(vtemp);
}
inline  void Bequation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,int* i,int *j,IMPACT_Dim *x1)
{
    B->cleanouter(vtemp,c,O,i,j,x1);
    if (c->NE()>0)
      for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
	E->cleanouter(vtemp,c,O,i,j,&x2);
    vzerochk(vtemp);
}
inline void Bequation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,int* i,int *j,IMPACT_Dim *x1)
{
   B->clean(vtemp,c,O,i,j,x1);
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->clean(vtemp,c,O,i,j,&x2);
    vzerochk(vtemp);
}
inline void Eequation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,int* i,int *j,IMPACT_Dim *x1)
{
  int ktemp=1;
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->cleanouter(vtemp,c,O,i,j,&x2);
   f1->cleanouter(vtemp,c,O,i,j,&ktemp,x1);
   if (c->NB()>0)
     for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
       B->cleanouter(vtemp,c,O,i,j,&x2);
   vzerochk(vtemp);
}
inline void Eequation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,int* i,int *j,IMPACT_Dim *x1)
{
   int ktemp=1;
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->cleanouter(vtemp,c,O,i,j,&x2);
   f1->clean(vtemp,c,O,i,j,&ktemp,x1);
   if (c->NB()>0)
     for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
       B->clean(vtemp,c,O,i,j,&x2);
   vzerochk(vtemp);
}
}

namespace no_zero_vector_checking
{
  int withzerovectorchecking=0;
inline void f0equation_clean_outer(IMPACT_ParVec * vtemp, IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var * B,IMPACT_Var* f1,int *i,int *j,int *k)
{
  f0->cleanouter(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
      E->cleanouter(vtemp,c,O,i,j,&x1);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->cleanouter(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    f1->cleanouter(vtemp,c,O,i,j,k,&x1);
}
inline void f0equation_clean_inner(IMPACT_ParVec * vtemp, IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var * B,IMPACT_Var* f1,int *i,int *j,int *k)
{
  f0->clean(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
      E->clean(vtemp,c,O,i,j,&x1);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->clean(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    f1->clean(vtemp,c,O,i,j,k,&x1);
  
}
inline void f1equation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,IMPACT_Var* f2,int* i,int *j,int *k,IMPACT_Dim *x1)
{
  f0->cleanouter(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
      E->cleanouter(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x2=1;x2<=c->Nf1();++x2)
    f1->cleanouter(vtemp,c,O,i,j,k,&x2);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->cleanouter(vtemp,c,O,i,j,&x2);
  if (c->Nf2()>0)
    for (int pos=0;pos<c->Nf2();++pos)
      f2->cleanouter(vtemp,c,O,i,j,k,&pos);
}
inline void f1equation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* f0,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,IMPACT_Var* f2,int* i,int *j,int *k,IMPACT_Dim *x1)
{
  f0->clean(vtemp,c,O,i,j,k);
  if (c->NE()>0)
    for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
      E->clean(vtemp,c,O,i,j,&x2);
  for (IMPACT_Dim x2=1;x2<=c->Nf1();++x2)
    f1->clean(vtemp,c,O,i,j,k,&x2);
  for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
    B->clean(vtemp,c,O,i,j,&x2);
  if (c->Nf2()>0)
    for (int pos=0;pos<c->Nf2();++pos)
      f2->clean(vtemp,c,O,i,j,k,&pos);
}
inline void f2equation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f0,IMPACT_Var* f1,IMPACT_Var* f2,IMPACT_Var* f3,int* i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  f0->cleanouter(vtemp,c,O,i,j,k);
  if (c->Nf2()>0)
 for (int pos=0;pos<c->Nf2();++pos)
    f2->cleanouter(vtemp,c,O,i,j,k,&pos);
if (c->Nf1()>0)
  for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
    f1->cleanouter(vtemp,c,O,i,j,k,&x3);
  for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
    E->cleanouter(vtemp,c,O,i,j,&x3);
  for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
    B->cleanouter(vtemp,c,O,i,j,&x3);
}
inline void f2equation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O, IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f0,IMPACT_Var* f1,IMPACT_Var* f2,IMPACT_Var* f3,int* i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  f0->clean(vtemp,c,O,i,j,k);
  if (c->Nf2()>0)
 for (int pos=0;pos<c->Nf2();++pos)
    f2->cleanouter(vtemp,c,O,i,j,k,&pos);
if (c->Nf1()>0)
  for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
    f1->clean(vtemp,c,O,i,j,k,&x3);
  for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
    E->clean(vtemp,c,O,i,j,&x3);
  for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
    B->clean(vtemp,c,O,i,j,&x3);
}
inline  void Bequation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,int* i,int *j,IMPACT_Dim *x1)
 {
   B->cleanouter(vtemp,c,O,i,j,x1);
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->cleanouter(vtemp,c,O,i,j,&x2);
 }
 inline void Bequation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,int* i,int *j,IMPACT_Dim *x1)
 {
   B->clean(vtemp,c,O,i,j,x1);
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->clean(vtemp,c,O,i,j,&x2);
 }
inline  void Eequation_clean_outer(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,int* i,int *j,IMPACT_Dim *x1)
 {
   int ktemp=1;
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->cleanouter(vtemp,c,O,i,j,&x2);
    f1->cleanouter(vtemp,c,O,i,j,&ktemp,x1);
   if (c->NB()>0)
     for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
       B->cleanouter(vtemp,c,O,i,j,&x2);
 }
inline  void Eequation_clean_inner(IMPACT_ParVec * vtemp,IMPACT_Config *c, IMPACT_StenOps * O,IMPACT_Var* E,IMPACT_Var* B,IMPACT_Var* f1,int* i,int *j,IMPACT_Dim *x1)
 {
   int ktemp=1;
   if (c->NE()>0)
     for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
       E->cleanouter(vtemp,c,O,i,j,&x2);
   f1->clean(vtemp,c,O,i,j,&ktemp,x1);
   if (c->NB()>0)
     for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
       B->clean(vtemp,c,O,i,j,&x2);
 }
}
