/*
**********************************************************
Functions for packing sparse matrix quickly

Version 1.0

AGRT

8/3 
This is painful.
**********************************************************
*/
 //_________________________________________________________________
  //Counts non zeros in a row and packs them:


inline void PackSparse_int(IMPACT_ParSpa *S,IMPACT_ParVec *v,int packstart,int packend,int *runningtotal)
{
   double val=0.0;
  for (int i=packstart;i<=packend;++i) //Now pack sparse matrix
      {
	val=v->Get(i);
	if(fabs(val)>zerotolerance::zerothreshold)
	  {
	    S->setval(*runningtotal+1,val);
	    S->setcol(*runningtotal+1,i);
	    *runningtotal+=1;
	  }
      }
}
inline void PackSparse(IMPACT_ParSpa *S,IMPACT_ParVec *v,int i,int *runningtotal)
{
  double val;
  val=v->Get(i);
  if(fabs(val)>zerotolerance::zerothreshold)
	  {
	    S->setval(*runningtotal+1,val);
	    S->setcol(*runningtotal+1,i);
	    *runningtotal+=1;
	  }
}
inline void PackSparse_outer(IMPACT_ParSpa *S,IMPACT_ParVec *v,int i,int *i_old,int *runningtotal,int *rownumber)
{
  double val;
  int check=0;
  val=v->Get(*i_old); //get the shifted i...
  if(fabs(val)>zerotolerance::zerothreshold)
    {// (first cycle over sparse to check there are no repeats)
      for (int index=S->getrow(*rownumber);index<=*runningtotal;++index)
	if (S->getcol(index)==i)
	  {
	    val+=S->getval(index);
	    if(fabs(val)>zerotolerance::zerothreshold)
	      S->setval(index,val); //increment existing
	    else 
	      {
		for (int innerindex=index;innerindex<=*runningtotal;++innerindex)
		  { 		    
		    S->setval(innerindex,S->getval(innerindex+1));
		    S->setcol(innerindex,S->getcol(innerindex+1));
		    }
		S->setval(index,0.0);
		S->setcol(index,0);
		*runningtotal-=1;
	      }
	    check=1;
	    break;
	  }
      if (check==0)
      {
	  S->setval(*runningtotal+1,val);
	  S->setcol(*runningtotal+1,i); //but put in the periodic value
	  *runningtotal+=1;
      }
      
    }
}
inline void f0equation_Pack_inner(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *E,IMPACT_Var * B,int *i,int *j,int *k)
  {
    int rownumber=f0->getrow(i,j,k);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
    int runningtotal=index-1; //to keep tabs of the position in sparse.
    int Nv=c->Nv();
    int one=1;
    //First, appropriate f0 components - integral over k
    PackSparse_int(S,v,f0->getrow(i,j,&one),f0->getrow(i,j,&Nv),&runningtotal);

    // Now f0 components - 5point stencil: - due to ion motion
    // PackSparse(S,v,f0->getrow(i,j,k),&runningtotal); - no because will be 
    //counted twice.
     PackSparse(S,v,f0->getrow(&iplus,j,k),&runningtotal);
      PackSparse(S,v,f0->getrow(&iminus,j,k),&runningtotal);
      PackSparse(S,v,f0->getrow(i,&jplus,k),&runningtotal);
      PackSparse(S,v,f0->getrow(i,&jminus,k),&runningtotal);

    // Now f1 components - 5point stencil:
    for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    {
      PackSparse(S,v,f1->getrow(i,j,k,&x1),&runningtotal);
      PackSparse(S,v,f1->getrow(&iplus,j,k,&x1),&runningtotal);
      PackSparse(S,v,f1->getrow(&iminus,j,k,&x1),&runningtotal);
      PackSparse(S,v,f1->getrow(i,&jplus,k,&x1),&runningtotal);
      PackSparse(S,v,f1->getrow(i,&jminus,k,&x1),&runningtotal);
    }
    //lastly E:
    if (c->NE()>0)
      for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
	PackSparse(S,v,E->getrow(i,j,&x1),&runningtotal);
    for (IMPACT_Dim x4=3;x4>3-c->NB();--x4)
      PackSparse(S,v,B->getrow(i,j,&x4),&runningtotal);
    
    chksparse(S,&rownumber,&runningtotal);
  }
inline void f0equation_Pack_outer(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *E,IMPACT_Var * B,int *i,int *j,int *k)
  {
    int rownumber=f0->getrow(i,j,k);
    int index = S->getrow(rownumber); 
    
    int runningtotal=index-1; //to keep tabs of the position in sparse.
    int Nv=c->Nv();
    int one=1;
    //First, appropriate f0 components - integral over k
    PackSparse_int(S,v,f0->getrow(i,j,&one),f0->getrow(i,j,&Nv),&runningtotal);
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
    int f0_iplus=f0->getrow(&iplus,j,k); 
    // These are to pass to packsparse outer
    int f0_iminus=f0->getrow(&iminus,j,k);
    chk_ij(c,&iplus,&jplus); //because if shifted is already counted
    chk_ij(c,&iminus,&jminus);
    // Now f0 components - 5point stencil:
    // PackSparse(S,v,f0->getrow(i,j,k),&runningtotal); //no because 2 count
    if (jplus!=*j) PackSparse(S,v,f0->getrow(i,&jplus,k),&runningtotal);
    if (jminus!=*j) PackSparse(S,v,f0->getrow(i,&jminus,k),&runningtotal);
    if (iplus!=*i) PackSparse_outer(S,v,f0->getrow(&iplus,j,k),&f0_iplus,&runningtotal,&rownumber);
    if (iminus!=*i) PackSparse_outer(S,v,f0->getrow(&iminus,j,k),&f0_iminus,&runningtotal,&rownumber);

    // Now f1 components - 5point stencil:
    for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
    {
      PackSparse(S,v,f1->getrow(i,j,k,&x1),&runningtotal);
      iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
      int f1_iplus=f1->getrow(&iplus,j,k,&x1);
      int f1_iminus=f1->getrow(&iminus,j,k,&x1);
      chk_ij(c,&iplus,&jplus);
      chk_ij(c,&iminus,&jminus);
      PackSparse(S,v,f1->getrow(i,&jplus,k,&x1),&runningtotal);
      PackSparse(S,v,f1->getrow(i,&jminus,k,&x1),&runningtotal);
      PackSparse_outer(S,v,f1->getrow(&iplus,j,k,&x1),&f1_iplus,&runningtotal,&rownumber);
      PackSparse_outer(S,v,f1->getrow(&iminus,j,k,&x1),&f1_iminus,&runningtotal,&rownumber);    
    }
    //lastly E:
    if (c->NE()>0)
      for (IMPACT_Dim x1=1;x1<=c->NE();++x1)
	PackSparse(S,v,E->getrow(i,j,&x1),&runningtotal);
    for (IMPACT_Dim x4=3;x4>3-c->NB();--x4)
      PackSparse(S,v,B->getrow(i,j,&x4),&runningtotal);
    chksparse(S,&rownumber,&runningtotal);
  }

inline void f1equation_Pack_inner(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var* f2,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,int *k,IMPACT_Dim *x1)
  {
   
    int rownumber= f1->getrow(i,j,k,x1);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
    int runningtotal=index-1; //to keep tabs of the position in sparse.
    int Nv=c->Nv();
    int one=1;
    //First, appropriate f1 components - integral over k
    PackSparse_int(S,v,f1->getrow(i,j,&one,x1),f1->getrow(i,j,&Nv,x1),&runningtotal);
    IMPACT_Dim x2,x3;
    GetOrthogonal(x1,&x2,&x3);
    if (x2<=c->Nf1()) PackSparse(S,v,f1->getrow(i,j,k,&x2),&runningtotal);
    if (x3<=c->Nf1()) PackSparse(S,v,f1->getrow(i,j,k,&x3),&runningtotal);
   
    // Now f0 components - 5point stencil:
      PackSparse(S,v,f0->getrow(i,j,k),&runningtotal);
      PackSparse(S,v,f0->getrow(&iplus,j,k),&runningtotal);
      PackSparse(S,v,f0->getrow(&iminus,j,k),&runningtotal);
      PackSparse(S,v,f0->getrow(i,&jplus,k),&runningtotal);
      PackSparse(S,v,f0->getrow(i,&jminus,k),&runningtotal);
      // and f2 components in 5point stencil also
      
  for (int pos=0;pos<c->Nf2();++pos)
    {
      PackSparse(S,v,f2->getrow(i,j,k,&pos),&runningtotal);
      PackSparse(S,v,f2->getrow(&iplus,j,k,&pos),&runningtotal);
      PackSparse(S,v,f2->getrow(&iminus,j,k,&pos),&runningtotal);
      PackSparse(S,v,f2->getrow(i,&jplus,k,&pos),&runningtotal);
      PackSparse(S,v,f2->getrow(i,&jminus,k,&pos),&runningtotal);
    }
    //Now E & B
    if (c->NE()>0)
      for (IMPACT_Dim x5=1;x5<=c->NE();++x5)
	PackSparse(S,v,E->getrow(i,j,&x5),&runningtotal);
    for (IMPACT_Dim x4=3;x4>3-c->NB();--x4)
      PackSparse(S,v,B->getrow(i,j,&x4),&runningtotal);
    
    chksparse(S,&rownumber,&runningtotal);
  }
inline void f1equation_Pack_outer(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var* f2,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,int *k,IMPACT_Dim *x1)
  {
    int rownumber= f1->getrow(i,j,k,x1);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
    int f0_iplus=f0->getrow(&iplus,j,k); 
// These are to pass to packsparse outer
    int f0_iminus=f0->getrow(&iminus,j,k);

    int *f2_iplus;
    f2_iplus = new int[c->Nf2()];
    // These are to pass to packsparse outer
    int *f2_iminus;
    f2_iminus = new int[c->Nf2()];
    
    for (int pos=0;pos<c->Nf2();++pos)
      {
	f2_iplus[pos]=f2->getrow(&iplus,j,k,&pos);
	f2_iminus[pos]=f2->getrow(&iminus,j,k,&pos);
      }
    chk_ij(c,&iplus,&jplus);
    chk_ij(c,&iminus,&jminus);
    int runningtotal=index-1; //to keep tabs of the position in sparse.
    int Nv=c->Nv();
    int one=1;
    //First, appropriate f1 components - integral over k
    PackSparse_int(S,v,f1->getrow(i,j,&one,x1),f1->getrow(i,j,&Nv,x1),&runningtotal);
    //And orthogonal bits!
    IMPACT_Dim x2,x3;
    GetOrthogonal(x1,&x2,&x3);
    if (x2<=c->Nf1()) PackSparse(S,v,f1->getrow(i,j,k,&x2),&runningtotal);
    if (x3<=c->Nf1()) PackSparse(S,v,f1->getrow(i,j,k,&x3),&runningtotal);
    
    // Now f0 components - 5point stencil:
    PackSparse(S,v,f0->getrow(i,j,k),&runningtotal);
    PackSparse(S,v,f0->getrow(i,&jplus,k),&runningtotal);
    PackSparse(S,v,f0->getrow(i,&jminus,k),&runningtotal);
    PackSparse_outer(S,v,f0->getrow(&iplus,j,k),&f0_iplus,&runningtotal,&rownumber);
    PackSparse_outer(S,v,f0->getrow(&iminus,j,k),&f0_iminus,&runningtotal,&rownumber);
    //and f2...
      for (int pos=0;pos<c->Nf2();++pos)
    {
      PackSparse(S,v,f2->getrow(i,j,k,&pos),&runningtotal);
      PackSparse(S,v,f2->getrow(i,&jplus,k,&pos),&runningtotal);
      PackSparse(S,v,f2->getrow(i,&jminus,k,&pos),&runningtotal);
      PackSparse_outer(S,v,f2->getrow(&iplus,j,k,&pos),&f2_iplus[pos],&runningtotal,&rownumber);
      PackSparse_outer(S,v,f2->getrow(&iminus,j,k,&pos),&f2_iminus[pos],&runningtotal,&rownumber);
    }
    //Now E & B
    if (c->NE()>0)
      for (IMPACT_Dim x5=1;x5<=c->NE();++x5)
	PackSparse(S,v,E->getrow(i,j,&x5),&runningtotal);
    for (IMPACT_Dim x4=3;x4>3-c->NB();--x4)
      PackSparse(S,v,B->getrow(i,j,&x4),&runningtotal);
    
    chksparse(S,&rownumber,&runningtotal);
    delete[] f2_iplus;
    delete[] f2_iminus;    
  }
inline void f2equation_Pack_inner(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,IMPACT_Var* f3,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
  {
    int Nv=c->Nv();
    int one=1;
    int rownumber= f2->getrow(i,j,k,x1,x2);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;

    chk_ij(c,&iplus,&jplus);
    chk_ij(c,&iminus,&jminus);
    int runningtotal=index-1; //to keep tabs of the position in sparse.

    // Now f1 components - 5point stencil:
    if (c->Nf1()>0)
      for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
	{
	  PackSparse(S,v,f1->getrow(i,j,k,&x3),&runningtotal);
	  PackSparse(S,v,f1->getrow(&iplus,j,k,&x3),&runningtotal);
	  PackSparse(S,v,f1->getrow(&iminus,j,k,&x3),&runningtotal);
	  PackSparse(S,v,f1->getrow(i,&jplus,k,&x3),&runningtotal);
	  PackSparse(S,v,f1->getrow(i,&jminus,k,&x3),&runningtotal);
	  
	  
	}
    // and f2 components
    for (int pos=0;pos<c->Nf2();++pos)
      {
	PackSparse(S,v,f2->getrow(i,j,k,&pos),&runningtotal);
	PackSparse(S,v,f2->getrow(&iplus,j,k,&pos),&runningtotal);
	PackSparse(S,v,f2->getrow(&iminus,j,k,&pos),&runningtotal);   
	PackSparse(S,v,f2->getrow(i,&jplus,k,&pos),&runningtotal);
	PackSparse(S,v,f2->getrow(i,&jminus,k,&pos),&runningtotal);
	//New bits
	PackSparse(S,v,f2->getrow(&iplus,&jplus,k,&pos),&runningtotal);
	PackSparse(S,v,f2->getrow(&iminus,&jplus,k,&pos),&runningtotal);   
	PackSparse(S,v,f2->getrow(&iplus,&jminus,k,&pos),&runningtotal);
	PackSparse(S,v,f2->getrow(&iminus,&jminus,k,&pos),&runningtotal);
      }
 
    //Now E & B
    if (c->NE()>0)
      for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
	PackSparse(S,v,E->getrow(i,j,&x3),&runningtotal);
    for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
      PackSparse(S,v,B->getrow(i,j,&x3),&runningtotal);
    
    chksparse(S,&rownumber,&runningtotal);

    //Finally, appropriate f0 components for IB heating - integral over k
    PackSparse_int(S,v,f0->getrow(i,j,&one),f0->getrow(i,j,&Nv),&runningtotal);

  }
inline void f2equation_Pack_outer(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,IMPACT_Var* f3,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,int *k,IMPACT_Dim *x1,IMPACT_Dim *x2)
  {
    int Nv=c->Nv();
    int one=1;
    int rownumber= f2->getrow(i,j,k,x1,x2);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;

    int *f2_iplus;
    f2_iplus = new int[c->Nf2()];
    // These are to pass to packsparse outer
    int *f2_iminus;
    f2_iminus = new int[c->Nf2()];

    for (int pos=0;pos<c->Nf2();++pos)
      {
	f2_iplus[pos]=f2->getrow(&iplus,j,k,&pos);
	f2_iminus[pos]=f2->getrow(&iminus,j,k,&pos);
      }
    
    int *f1_iplus;
    f1_iplus = new int[c->Nf1()];
    int *f1_iminus;
    f1_iminus = new int[c->Nf1()];
    
    for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
      {
	f1_iplus[x3.get()-1]=f1->getrow(&iplus,j,k,&x3);
	f1_iminus[x3.get()-1]=f1->getrow(&iminus,j,k,&x3);
      }
    chk_ij(c,&iplus,&jplus);
    chk_ij(c,&iminus,&jminus);
    int runningtotal=index-1; //to keep tabs of the position in sparse.

    // Now f1 components - 5point stencil:
if (c->Nf1()>0)
      for (IMPACT_Dim x3=1;x3<=c->Nf1();++x3)
	{
	  PackSparse(S,v,f1->getrow(i,j,k,&x3),&runningtotal);
	  PackSparse(S,v,f1->getrow(i,&jplus,k,&x3),&runningtotal);
	  PackSparse(S,v,f1->getrow(i,&jminus,k,&x3),&runningtotal);
	  PackSparse_outer(S,v,f1->getrow(&iplus,j,k,&x3),&f1_iplus[x3.get()-1],
			   &runningtotal,&rownumber);
	  PackSparse_outer(S,v,f1->getrow(&iminus,j,k,&x3),
			   &f1_iminus[x3.get()-1],&runningtotal,&rownumber);
	}
// and f2 components
 for (int pos=0;pos<c->Nf2();++pos)
   {
     PackSparse(S,v,f2->getrow(i,j,k,&pos),&runningtotal);
     PackSparse(S,v,f2->getrow(i,&jplus,k,&pos),&runningtotal);
     PackSparse(S,v,f2->getrow(i,&jminus,k,&pos),&runningtotal);
     PackSparse_outer(S,v,f2->getrow(&iplus,j,k,&pos),&f2_iplus[pos],
		      &runningtotal, &rownumber);
     PackSparse_outer(S,v,f2->getrow(&iminus,j,k,&pos),&f2_iminus[pos],
		      &runningtotal,&rownumber);
     //new bits
     /*PackSparse_outer(S,v,f2->getrow(&iplus,&jplus,k,&pos),&f2_iplus[pos],
		      &runningtotal, &rownumber);
     PackSparse_outer(S,v,f2->getrow(&iminus,&jplus,k,&pos),&f2_iminus[pos],
		      &runningtotal,&rownumber);
     PackSparse_outer(S,v,f2->getrow(&iplus,&jminus,k,&pos),&f2_iplus[pos],
		      &runningtotal, &rownumber);
     PackSparse_outer(S,v,f2->getrow(&iminus,&jminus,k,&pos),&f2_iminus[pos],
     &runningtotal,&rownumber);*/
     if (*i<c->Nx()&&*j<c->Ny()) PackSparse(S,v,f2->getrow(&iplus,&jplus,k,&pos),&runningtotal);
     if (*i<c->Nx()&&*j>1) PackSparse(S,v,f2->getrow(&iplus,&jminus,k,&pos),&runningtotal);
     if (*i>1&&*j<c->Ny()) PackSparse(S,v,f2->getrow(&iminus,&jplus,k,&pos),&runningtotal);
     if (*i>1&&*j>1) PackSparse(S,v,f2->getrow(&iminus,&jminus,k,&pos),&runningtotal);
   }
 
    //Now E & B
    if (c->NE()>0)
      for (IMPACT_Dim x3=1;x3<=c->NE();++x3)
	PackSparse(S,v,E->getrow(i,j,&x3),&runningtotal);
    for (IMPACT_Dim x3=3;x3>3-c->NB();--x3)
      PackSparse(S,v,B->getrow(i,j,&x3),&runningtotal);
    
    chksparse(S,&rownumber,&runningtotal);
    //Finally, appropriate f0 components for IB heating - integral over k
    PackSparse_int(S,v,f0->getrow(i,j,&one),f0->getrow(i,j,&Nv),&runningtotal);
    delete[] f2_iplus;
    delete[] f2_iminus;    
  }
inline void Eequation_Pack_inner(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f1,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,IMPACT_Dim *x1)
  {
    int rownumber= E->getrow(i,j,x1);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
    //chk_ij(c,&iplus,&jplus);
    //chk_ij(c,&iminus,&jminus);
    int runningtotal=index-1; //to keep tabs of the position in sparse.
    int Nv=c->Nv();
    int one=1;
    //First, appropriate f1 components - integral over k
    PackSparse_int(S,v,f1->getrow(i,j,&one,x1),f1->getrow(i,j,&Nv,x1),&runningtotal);

    //Now E
    if (c->NE()>0)
      for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
	{
	  PackSparse(S,v,E->getrow(i,j,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(&iplus,j,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(&iminus,j,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(i,&jplus,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(i,&jminus,&x2),&runningtotal);
	  //New bits
	  PackSparse(S,v,E->getrow(&iplus,&jplus,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(&iminus,&jplus,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(&iplus,&jminus,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(&iminus,&jminus,&x2),&runningtotal);
	}
    // finally B  - 5point stencil:
for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
  {
    PackSparse(S,v,B->getrow(i,j,&x2),&runningtotal);
    PackSparse(S,v,B->getrow(&iplus,j,&x2),&runningtotal);
    PackSparse(S,v,B->getrow(&iminus,j,&x2),&runningtotal);
    PackSparse(S,v,B->getrow(i,&jplus,&x2),&runningtotal);
    PackSparse(S,v,B->getrow(i,&jminus,&x2),&runningtotal);
  }
    
    
    chksparse(S,&rownumber,&runningtotal);
  }
inline void Eequation_Pack_outer(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *f1,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,IMPACT_Dim *x1)
  {
    int rownumber= E->getrow(i,j,x1);
    int index = S->getrow(rownumber); 
    
    int runningtotal=index-1; //to keep tabs of the position in sparse.
    int Nv=c->Nv();
    int one=1;
    //First, appropriate f1 components - integral over k
    PackSparse_int(S,v,f1->getrow(i,j,&one,x1),f1->getrow(i,j,&Nv,x1),&runningtotal);

    //Now E
    if (c->NE()>0)
      for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
	{
	  PackSparse(S,v,E->getrow(i,j,&x2),&runningtotal);
	  int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
	  int E_iplus=E->getrow(&iplus,j,&x2);
	  int E_iminus=E->getrow(&iminus,j,&x2);
	  chk_j(c,&jplus,&jminus);
	  int E_iplusjplus=E->getrow(&iplus,&jplus,&x2);
	  int E_iminusjplus=E->getrow(&iminus,&jplus,&x2);
	  int E_iplusjminus=E->getrow(&iplus,&jminus,&x2);
	  int E_iminusjminus=E->getrow(&iminus,&jminus,&x2);

	  chk_ij(c,&iplus,&jplus);
	  chk_ij(c,&iminus,&jminus);
	  PackSparse(S,v,E->getrow(i,&jplus,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(i,&jminus,&x2),&runningtotal);
	  PackSparse_outer(S,v,E->getrow(&iplus,j,&x2),&E_iplus,&runningtotal,&rownumber);
	  PackSparse_outer(S,v,E->getrow(&iminus,j,&x2),&E_iminus,&runningtotal,&rownumber);
	  //New Bits

	    PackSparse_outer(S,v,E->getrow(&iplus,&jplus,&x2),&E_iplusjplus,&runningtotal,&rownumber);
	    PackSparse_outer(S,v,E->getrow(&iminus,&jplus,&x2),&E_iminusjplus,&runningtotal,&rownumber);

	    PackSparse_outer(S,v,E->getrow(&iplus,&jminus,&x2),&E_iplusjminus,&runningtotal,&rownumber);
	    PackSparse_outer(S,v,E->getrow(&iminus,&jminus,&x2),&E_iminusjminus,&runningtotal,&rownumber);
	  	  
	  /*if (*i<c->Nx()&&*j<c->Ny())
	    PackSparse_outer(S,v,E->getrow(&iplus,&jplus,&x2),&E_iplus,&runningtotal,&rownumber);
	  if (*i<c->Nx()&&*j>1)
	    PackSparse_outer(S,v,E->getrow(&iplus,&jminus,&x2),&E_iplus,&runningtotal,&rownumber);
	  if (*i>1&&*j<c->Ny())
	    PackSparse_outer(S,v,E->getrow(&iminus,&jplus,&x2),&E_iminus,&runningtotal,&rownumber);
	  if (*i>1&&*j>1)
	  PackSparse(S,v,E->getrow(&iminus,&jminus,&x2),&E_iminus,&runningtotal,&rownumber);*/
	  /*	  if (*i<c->Nx()&&*j<c->Ny()) PackSparse(S,v,E->getrow(&iplus,&jplus,&x2),&runningtotal);
	  if (*i<c->Nx()&&*j>1) PackSparse(S,v,E->getrow(&iplus,&jminus,&x2),&runningtotal);
	  if (*i>1&&*j<c->Ny()) PackSparse(S,v,E->getrow(&iminus,&jplus,&x2),&runningtotal);
	  if (*i>1&&*j>1) PackSparse(S,v,E->getrow(&iminus,&jminus,&x2),&runningtotal);*/

	  
	}
    // finally B  - 5point stencil:
    for (IMPACT_Dim x2=3;x2>3-c->NB();--x2)
      {
	PackSparse(S,v,B->getrow(i,j,&x2),&runningtotal);
	int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
	int B_iplus=B->getrow(&iplus,j,&x2);
	int B_iminus=B->getrow(&iminus,j,&x2);
	chk_ij(c,&iplus,&jplus);
	chk_ij(c,&iminus,&jminus);
	PackSparse(S,v,B->getrow(i,&jplus,&x2),&runningtotal);
	PackSparse(S,v,B->getrow(i,&jminus,&x2),&runningtotal);
	PackSparse_outer(S,v,B->getrow(&iplus,j,&x2),&B_iplus,&runningtotal,&rownumber);
	PackSparse_outer(S,v,B->getrow(&iminus,j,&x2),&B_iminus,&runningtotal,&rownumber);
	
      }
    
    
    chksparse(S,&rownumber,&runningtotal);
  }
inline void Bequation_Pack_inner(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,IMPACT_Dim *x1)
  {
    int rownumber= B->getrow(i,j,x1);
    int index = S->getrow(rownumber); 
    int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
    int runningtotal=index-1; //to keep tabs of the position in sparse.

    //Now that single B element
    PackSparse(S,v,B->getrow(i,j,x1),&runningtotal);

    // finally E  - 5point stencil:
if (c->NE()>0)
  for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
    {
      PackSparse(S,v,E->getrow(i,j,&x2),&runningtotal);
      PackSparse(S,v,E->getrow(&iplus,j,&x2),&runningtotal);
      PackSparse(S,v,E->getrow(&iminus,j,&x2),&runningtotal);
      PackSparse(S,v,E->getrow(i,&jplus,&x2),&runningtotal);
      PackSparse(S,v,E->getrow(i,&jminus,&x2),&runningtotal);
    }
 
    
    chksparse(S,&rownumber,&runningtotal);
  }
inline void Bequation_Pack_outer(IMPACT_ParVec *v, IMPACT_ParSpa *S, IMPACT_Config *c,IMPACT_Var *E,IMPACT_Var *B,int *i,int *j,IMPACT_Dim *x1)
  {
    int rownumber= B->getrow(i,j,x1);
    int index = S->getrow(rownumber); 
    int runningtotal=index-1; //to keep tabs of the position in sparse.

    //Now that single B element
    PackSparse(S,v,B->getrow(i,j,x1),&runningtotal);

    // finally E  - 5point stencil:
    if (c->NE()>0)
      for (IMPACT_Dim x2=1;x2<=c->NE();++x2)
	{
	  PackSparse(S,v,E->getrow(i,j,&x2),&runningtotal);
	  int iminus=*i-1,iplus=*i+1,jplus=*j+1,jminus=*j-1;
	  int E_iplus=E->getrow(&iplus,j,&x2);
	  int E_iminus=E->getrow(&iminus,j,&x2);
	  chk_ij(c,&iplus,&jplus);
	  chk_ij(c,&iminus,&jminus);
	  PackSparse(S,v,E->getrow(i,&jplus,&x2),&runningtotal);
	  PackSparse(S,v,E->getrow(i,&jminus,&x2),&runningtotal);
	  PackSparse_outer(S,v,E->getrow(&iplus,j,&x2),&E_iplus,&runningtotal,
			   &rownumber);
	  PackSparse_outer(S,v,E->getrow(&iminus,j,&x2),&E_iminus,
			   &runningtotal,&rownumber);
	  
	}
    
    
    chksparse(S,&rownumber,&runningtotal);
  }
