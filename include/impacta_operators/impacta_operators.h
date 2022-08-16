/*
**********************************************************
IMPACT Version 1.9
All the operator namespaces
Basically the namespaces contain different bits of code depending
on the boundary conditions required.

Version 1.6
1
AGRT

12/2/07

19/02/07 - note can change d/dv term to make relativistic 

5/3/07

IMPORTANT CHANGE!!! IMPACT_Vector replaced with IMPACT_ParVec
-> parallel vector class.

Dimensions new class IMPACT_Dim

8/3 - chk_ij moved here.

13/3 - Boundary condition code (periodic boundary conditions)
 updated for parallel version (see _BC versions of code).

25/5/07 - Curl operator added for non matrix method (explicit)

4/6/07 - f2 functions added
**********************************************************

*/

int Open_Bound_x()
{
  return Boundaries_open_x;
}
int Open_Bound_y()
{
  return Boundaries_open_y;
}

 inline void chk_ij(IMPACT_Config *c,int *i,int *j)
{
  if (*i>c->Nx()) *i=1;
  if (*i<1) *i=c->Nx();
  if (*j>c->Ny()) *j=1;
  if (*j<1) *j=c->Ny();
}
 inline void chk_j(IMPACT_Config *c,int *jplus,int *jminus)
{
    if (*jplus>c->Ny())*jplus=1;
    if (*jminus<1)*jminus=c->Ny();
}
 // Methods for inserting stencil into ParVec v
  //NOT AT BOUNDARIES!!!
  // remeber i,j,k are location indices
  // x1,x2,x3  are the indices of f/e/b e.g. xyy
  /* so in stencil 0 - (i,j)
     1 - (i+1,j)
     2 - (i-1,j)
     3 - (i,j+1)
     4 - (i,j-1) 
   */
  inline void IMPACT_Var::insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int* i, int *j, int *k)
  {
    int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
    set(v,get(v,i,j,k)+stencil->get(0),i,j,k);
    set(v,get(v,&iplus,j,k)+stencil->get(1),&iplus,j,k);
    set(v,get(v,&iminus,j,k)+stencil->get(2),&iminus,j,k);
    set(v,get(v,i,&jplus,k)+stencil->get(3),i,&jplus,k);
    set(v,get(v,i,&jminus,k)+stencil->get(4),i,&jminus,k);
  }
  inline void IMPACT_Var::insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1)
  {
    int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
    set(v,get(v,i,j,k,x1)+stencil->get(0),i,j,k,x1);
    set(v,get(v,&iplus,j,k,x1)+stencil->get(1),&iplus,j,k,x1);
    set(v,get(v,&iminus,j,k,x1)+stencil->get(2),&iminus,j,k,x1);
    set(v,get(v,i,&jplus,k,x1)+stencil->get(3),i,&jplus,k,x1);
    set(v,get(v,i,&jminus,k,x1)+stencil->get(4),i,&jminus,k,x1);
  }
 inline void IMPACT_Var::insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1)
  {
    int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
    set(v,get(v,i,j,x1)+stencil->get(0),i,j,x1);
    set(v,get(v,&iplus,j,x1)+stencil->get(1),&iplus,j,x1);
    set(v,get(v,&iminus,j,x1)+stencil->get(2),&iminus,j,x1);
    set(v,get(v,i,&jplus,x1)+stencil->get(3),i,&jplus,x1);
    set(v,get(v,i,&jminus,x1)+stencil->get(4),i,&jminus,x1);
  }
// This is a 3 point stencil only, in y direction
 inline void IMPACT_Var::insert_yonly(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1)
  {
    int jplus=*j+1,jminus=*j-1;
    set(v,get(v,i,j,x1)+stencil->get(0),i,j,x1);
    set(v,get(v,i,&jplus,x1)+stencil->get(3),i,&jplus,x1);
    set(v,get(v,i,&jminus,x1)+stencil->get(4),i,&jminus,x1);
  }
inline void IMPACT_Var::insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  set(v,get(v,i,j,k,x1,x2)+stencil->get(0),i,j,k,x1,x2);
  set(v,get(v,&iplus,j,k,x1,x2)+stencil->get(1),&iplus,j,k,x1,x2);
  set(v,get(v,&iminus,j,k,x1,x2)+stencil->get(2),&iminus,j,k,x1,x2);
  set(v,get(v,i,&jplus,k,x1,x2)+stencil->get(3),i,&jplus,k,x1,x2);
  set(v,get(v,i,&jminus,k,x1,x2)+stencil->get(4),i,&jminus,k,x1,x2);
}
/*
  void IMPACT_Var::insert(IMPACT_stencil *stencil, IMPACT_ParVec *v, int i, int j, int k,int x1,int x2,int x3)
  {
    set(v,get(v,i,j,k,x1,x2,x3)+stencil->get(0),i,j,k,x1,x2,x3);
    set(v,get(v,i+1,j,k,x1,x2,x3)+stencil->get(1),i+1,j,k,x1,x2,x3);
    set(v,get(v,i-1,j,k,x1,x2,x3)+stencil->get(2),i-1,j,k,x1,x2,x3);
    set(v,get(v,i,j+1,k,x1,x2,x3)+stencil->get(3),i,j+1,k,x1,x2,x3);
    set(v,get(v,i,j-1,k,x1,x2,x3)+stencil->get(4),i,j-1,k,x1,x2,x3);
    }*/

//___________________________________________________________________________________________
//Now code for the boundaries

/*
  This is code for periodic x and periodic y boundaries for f0!
 */

inline void IMPACT_Var::insert_BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k)
  {
    int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;

    IMPACT_f0_bound(stencil,Nx,Ny,&iplus,&iminus,&jplus,&jminus);

    set(v,get(v,i,j,k)+stencil->get(0),i,j,k);
    set(v,get(v,&iplus,j,k)+stencil->get(1),&iplus,j,k);
    set(v,get(v,&iminus,j,k)+stencil->get(2),&iminus,j,k);
    set(v,get(v,i,&jplus,k)+stencil->get(3),i,&jplus,k);
    set(v,get(v,i,&jminus,k)+stencil->get(4),i,&jminus,k);
  }


/*
  This is code for periodic x and periodic y boundaries for f1!
 */
 inline  void IMPACT_Var::insert_f1BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1)
  {
    int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;

    IMPACT_f1_E_bound(stencil,Nx,Ny,&iplus,&iminus,&jplus,&jminus,x1);
    
     set(v,get(v,i,j,k,x1)+stencil->get(0),i,j,k,x1);
     set(v,get(v,&iplus,j,k,x1)+stencil->get(1),&iplus,j,k,x1);
     set(v,get(v,&iminus,j,k,x1)+stencil->get(2),&iminus,j,k,x1);
     set(v,get(v,i,&jplus,k,x1)+stencil->get(3),i,&jplus,k,x1);
     set(v,get(v,i,&jminus,k,x1)+stencil->get(4),i,&jminus,k,x1);
  }

/*
  This is code for periodic x and periodic y boundaries for E!
 */

  inline void IMPACT_Var::insert_EBC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, IMPACT_Dim *x1)
  {
    int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
    
    IMPACT_f1_E_bound(stencil,Nx,Ny,&iplus,&iminus,&jplus,&jminus,x1);
    
    set(v,get(v,i,j,x1)+stencil->get(0),i,j,x1);
    set(v,get(v,&iplus,j,x1)+stencil->get(1),&iplus,j,x1);
    set(v,get(v,&iminus,j,x1)+stencil->get(2),&iminus,j,x1);
    set(v,get(v,i,&jplus,x1)+stencil->get(3),i,&jplus,x1);
    set(v,get(v,i,&jminus,x1)+stencil->get(4),i,&jminus,x1);
  }

  inline void IMPACT_Var::insert_yonly_EBC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, IMPACT_Dim *x1)
  {

    int itemp=*i,jplus=*j+1,jminus=*j-1;
    
    // set temporary values - these multiply the elements 0,3,4 to take account of the fact i can be over a boundary and
    // so with reflective the sign may be flipped.
    stencil->set(1,1.0); // represents lower boundary
    stencil->set(2,1.0); // represents upper boundary
    
    IMPACT_f1_E_bound(stencil,Nx,Ny,&itemp,&itemp,&jplus,&jminus,x1);
    
    double signflip = (stencil->get(1))*(stencil->get(2));
 
    set(v,get(v,&itemp,j,x1)+signflip*(stencil->get(0)),&itemp,j,x1);
    set(v,get(v,&itemp,&jplus,x1)+signflip*(stencil->get(3)),&itemp,&jplus,x1);
    set(v,get(v,&itemp,&jminus,x1)+signflip*(stencil->get(4)),&itemp,&jminus,x1);
  }

/*
  This is code for periodic x and periodic y boundaries for B!
 */

inline void IMPACT_Var::insert_BBC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, IMPACT_Dim *x1)
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  
  IMPACT_B_bound(stencil,Nx,Ny,&iplus,&iminus,&jplus,&jminus,x1);
  
  set(v,get(v,i,j,x1)+stencil->get(0),i,j,x1);
  set(v,get(v,&iplus,j,x1)+stencil->get(1),&iplus,j,x1);
  set(v,get(v,&iminus,j,x1)+stencil->get(2),&iminus,j,x1);
  set(v,get(v,i,&jplus,x1)+stencil->get(3),i,&jplus,x1);
  set(v,get(v,i,&jminus,x1)+stencil->get(4),i,&jminus,x1);
}

inline void IMPACT_Var::insert_BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int *i, int *j, int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  IMPACT_f2_bound(stencil,Nx,Ny,&iplus,&iminus,&jplus,&jminus,x1,x2);
  set(v,get(v,i,j,k,x1,x2)+stencil->get(0),i,j,k,x1,x2);
  if (fabs(stencil->get(1))>0.0)
    set(v,get(v,&iplus,j,k,x1,x2)+stencil->get(1),&iplus,j,k,x1,x2);
  if (fabs(stencil->get(2))>0.0)
    set(v,get(v,&iminus,j,k,x1,x2)+stencil->get(2),&iminus,j,k,x1,x2);
  set(v,get(v,i,&jplus,k,x1,x2)+stencil->get(3),i,&jplus,k,x1,x2);
  set(v,get(v,i,&jminus,k,x1,x2)+stencil->get(4),i,&jminus,k,x1,x2);
}
/*
  void IMPACT_Var::insert_BC(IMPACT_stencil *stencil, IMPACT_ParVec *v, int i, int j, int k,IMPACT_Dim *x1,IMPACT_Dim *x2,IMPACT_Dim *x3)
  {
    set(v,get(v,i,j,k,x1,x2,x3)+stencil->get(0),i,j,k,x1,x2,x3);
    set(v,get(v,i+1,j,k,x1,x2,x3)+stencil->get(1),i+1,j,k,x1,x2,x3);
    set(v,get(v,i-1,j,k,x1,x2,x3)+stencil->get(2),i-1,j,k,x1,x2,x3);
    set(v,get(v,i,j+1,k,x1,x2,x3)+stencil->get(3),i,j+1,k,x1,x2,x3);
    set(v,get(v,i,j-1,k,x1,x2,x3)+stencil->get(4),i,j-1,k,x1,x2,x3);
    }*/
//******************************************************************
// Now similar methods for inserting integrals instead.
//******************************************************************
inline void IMPACT_Var::insert(IMPACT_vint *vint, IMPACT_ParVec *v, int *i, int *j)
{
  for (int k=1;k<=vint->Nv();++k)
    {
      set(v,get(v,i,j,&k)+vint->get(k),i,j,&k);
    } 
}
  inline void IMPACT_Var::insert(IMPACT_vint *vint, IMPACT_ParVec *v, int *i, int *j,IMPACT_Dim *x1)
{
  for (int k=1;k<=vint->Nv();++k)
    {
      set(v,get(v,i,j,&k,x1)+vint->get(k),i,j,&k,x1);
    } 
}
//******************************************************************
// Velocity stencil operations
//******************************************************************

// for f0
inline void IMPACT_Var::insert(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k)
{
  int kminus=*k-1,kplus=*k+1;
  inc(v,vs->get(0),i,j,&kminus);
  inc(v,vs->get(1),i,j,k);
  inc(v,vs->get(2),i,j,&kplus);
}
inline void IMPACT_Var::insert_BC(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k)
{
  int kminus=*k-1,kplus=*k+1;
  if (kminus<1) ++kminus;
  if (kplus>Nv) --kplus;
  inc(v,vs->get(0),i,j,&kminus);
  inc(v,vs->get(1),i,j,k);
  inc(v,vs->get(2),i,j,&kplus);
}


// for f1
inline void IMPACT_Var::insert(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k, IMPACT_Dim *x1)
{
  int kminus=*k-1,kplus=*k+1;
  inc(v,vs->get(0),i,j,&kminus,x1);
  inc(v,vs->get(1),i,j,k,x1);
  inc(v,vs->get(2),i,j,&kplus,x1);
}
inline void IMPACT_Var::insert_BC(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k, IMPACT_Dim *x1)
{
  int kminus=*k-1,kplus=*k+1;
  if (kminus<1) ++kminus;
  if (kplus>Nv) --kplus;
  inc(v,vs->get(0),i,j,&kminus,x1);
  inc(v,vs->get(1),i,j,k,x1);
  inc(v,vs->get(2),i,j,&kplus,x1);
}

// for f2
inline void IMPACT_Var::insert(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k, IMPACT_Dim *x1, IMPACT_Dim *x2)
{
  int kminus=*k-1,kplus=*k+1;
  inc(v,vs->get(0),i,j,&kminus,x1,x2);
  inc(v,vs->get(1),i,j,k,x1,x2);
  inc(v,vs->get(2),i,j,&kplus,x1,x2);
}
inline void IMPACT_Var::insert_BC(IMPACT_Vel_Sten *vs, IMPACT_ParVec *v, int *i, int *j, int *k, IMPACT_Dim *x1, IMPACT_Dim *x2)
{
  int kminus=*k-1,kplus=*k+1;
  if (kminus<1) ++kminus;
  if (kplus>Nv) --kplus;
  inc(v,vs->get(0),i,j,&kminus,x1,x2);
  inc(v,vs->get(1),i,j,k,x1,x2);
  inc(v,vs->get(2),i,j,&kplus,x1,x2);
}








//This is for the non-stencil differentials wrt v
//namespace Diff_df0dv_centred //center differenced df0/dv term

  inline double IMPACT_Var::ddv(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k)// returns the local value of df0/dv
  {
    double answer;
    int kplus=*k+1,kminus=*k-1;
    answer = (get(v,i,j,&kplus)-get(v,i,j,&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember idv is the inverse of dv[k];
    return answer;
  }
 inline  double IMPACT_Var::ddv_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k)// returns the local value of df0/dv
  {
    int kplus=*k+1,kminus=*k-1;
    if (kplus>Nv) kplus=Nv;
    if (kminus<1) kminus=1;
    double answer;
    answer = (get(v,i,j,&kplus)-get(v,i,j,&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember idv is the inverse of dv[k];
    return answer;
  }
 inline  double IMPACT_Var::ddv_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k, IMPACT_Dim *x1)// returns the local value of df0/dv
  {
    int kplus=*k+1,kminus=*k-1;
    if (kplus>Nv) kplus=Nv;
    if (kminus<1) kminus=1;
    double answer;
    answer = (get(v,i,j,&kplus,x1)-get(v,i,j,&kminus,x1))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember idv is the inverse of dv[k];
    return answer;
  }

// Now ddv1/v for f2 equation - of f1
inline double IMPACT_Var::ddv1overv(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k,IMPACT_Dim *x1)// returns the local value of df1/dv
  {
    double answer;
    int kplus=*k+1,kminus=*k-1;
    answer = (get(v,i,j,&kplus,x1)/config1->v(&kplus)-get(v,i,j,&kminus,x1)/config1->v(&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember idv is the inverse of dv[k];
    return answer;
  }
inline  double IMPACT_Var::ddv1overv_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k,IMPACT_Dim *x1)// returns the local value of df1/dv
  {
    int kplus=*k+1,kminus=*k-1;
    if (kplus>Nv) kplus=Nv;
    if (kminus<1) kminus=1;
    double answer;
    answer = (get(v,i,j,&kplus,x1)/config1->v(&kplus)-get(v,i,j,&kminus,x1)/config1->v(&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember idv is the inverse of dv[k];
    return answer;
  }

 inline  double IMPACT_Var::ddv_v2(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k,IMPACT_Dim *x1)// returns the local value of d/dv(f1v^2)
  {
    double answer;
    int kplus=*k+1,kminus=*k-1;
    answer = (get(v,i,j,&kplus,x1)*config1->v2(&kplus)-get(v,i,j,&kminus,x1)*config1->v2(&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember to fix divide by at some point
    return answer;
  }
  inline double IMPACT_Var::ddv_v2_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k, IMPACT_Dim *x1)// returns the local value of d/dv(f1v^2) at boundary
  {
    int kplus=*k+1,kminus=*k-1;
    if (kplus>Nv) kplus=Nv;
    if (kminus<1) kminus=1;
    double answer;
    answer = (get(v,i,j,&kplus,x1)*config1->v2(&kplus)-get(v,i,j,&kminus,x1)*config1->v2(&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
  // remember to fix divide by at some point
    return answer;
  }

 inline  double IMPACT_Var::ddv_v3(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)// returns the local value of d/dv(f1v^2)
  {
    double answer;
    int kplus=*k+1,kminus=*k-1;
    answer = (get(v,i,j,&kplus,x1,x2)*config1->v3(&kplus)-get(v,i,j,&kminus,x1,x2)*config1->v3(&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
    // remember to fix divide by at some point
    return answer;
  }
  inline double IMPACT_Var::ddv_v3_BC(IMPACT_ParVec * v, IMPACT_Config *config1, int *i, int *j, int *k, IMPACT_Dim *x1,IMPACT_Dim *x2)// returns the local value of d/dv(f1v^2) at boundary
  {
    int kplus=*k+1,kminus=*k-1;
    if (kplus>Nv) kplus=Nv;
    if (kminus<1) kminus=1;
    double answer;
    answer = (get(v,i,j,&kplus,x1,x2)*config1->v3(&kplus)-get(v,i,j,&kminus,x1,x2)*config1->v3(&kminus))/(config1->dv(&kplus)+config1->dv(&kminus));
  // remember to fix divide by at some point
    return answer;
  }
inline int CeviLevita(IMPACT_Dim *i, IMPACT_Dim *j, IMPACT_Dim *k) //i.e. epsilon_ijk
{
  int answer=0;
  IMPACT_Dim total = *i+*j+*k;

  switch (total.get())
    {
    case 6: //clever eh!? 
      if (*i==2&&*j==2&&*k==2) break;
      answer=(*i-*k).get();
      if (*i-*k>1||*i-*k<-1) {answer/=-2;break;} //doesn't seem like it would work
      break;                                 //......but it does!
    default:
      break;
    }
  return answer;
}
inline int Kronecker_Delta(IMPACT_Dim *i,IMPACT_Dim *j)
{
  int answer=(i->get()==j->get());
  return answer;
}
//This curl is for E or B only
inline double IMPACT_Var::Curl(IMPACT_ParVec * v,IMPACT_Config *c,IMPACT_StenOps * O, int *i, int *j,  IMPACT_Dim * x1)
{
  // This cyclicly works out orthoganol correct coordinates for x1 andx2
  IMPACT_Dim x2;
  IMPACT_Dim x3;
  GetOrthogonal(x1,&x2,&x3);

  IMPACT_stencil temp_sten;
  double answer=0.0;
  int iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;

  if ((objecttype==5&&x2>3-c->NB())||(objecttype==4&&x2<=c->NE()))
    {
      temp_sten =(*O->ddxi(i,j,&x3))*(-1.0);
      if (objecttype==5)
	IMPACT_B_bound(&temp_sten,Nx,Ny,&iplus,&iminus,&jplus,&jminus,&x2);
      if (objecttype==4)
	IMPACT_f1_E_bound(&temp_sten,Nx,Ny,&iplus,&iminus,&jplus,&jminus,&x2);
      
      answer+=get(v,i,j,&x2)*temp_sten(0);
      answer+=get(v,&iplus,j,&x2)*temp_sten(1);
      answer+=get(v,&iminus,j,&x2)*temp_sten(2);
      answer+=get(v,i,&jplus,&x2)*temp_sten(3);
      answer+=get(v,i,&jminus,&x2)*temp_sten(4);
    }
  iplus=*i+1,iminus=*i-1,jplus=*j+1,jminus=*j-1;
  if ((objecttype==5&&x3>3-c->NB())||(objecttype==4&&x3<=c->NE()))
    {
      temp_sten = (*O->ddxi(i,j,&x2));
      if (objecttype==5)
	IMPACT_B_bound(&temp_sten,Nx,Ny,&iplus,&iminus,&jplus,&jminus,&x3);
      if (objecttype==4)
	IMPACT_f1_E_bound(&temp_sten,Nx,Ny,&iplus,&iminus,&jplus,&jminus,&x3);
      
      answer+=get(v,i,j,&x3)*temp_sten(0);
      answer+=get(v,&iplus,j,&x3)*temp_sten(1);
      answer+=get(v,&iminus,j,&x3)*temp_sten(2);
      answer+=get(v,i,&jplus,&x3)*temp_sten(3);
      answer+=get(v,i,&jminus,&x3)*temp_sten(4);
    }
  
  
  return answer;
}
