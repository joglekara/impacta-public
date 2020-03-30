//Choose the correct moment
void IMPACT_Select_Moment(struct IMPACT_Moment * answer,int choice,
			  IMPACT_Config *c, 
			  IMPACT_MPI_Config *MPIc, IMPACT_ParVec *v,
			  IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,
			  IMPACT_Var *f3,
			  IMPACT_Var *E,IMPACT_Var *B)
{
  answer->values.ChangeSize(c->Nx(),c->Ny());
  switch (choice)
    {
    case VISIT_PARAMS::data_id0:
      for (int j=1;j<=c->Ny();++j)
	for (int i=1;i<=c->Nx();++i)
	  answer->values.set(i,j,0.0);
      break;
    case VISIT_PARAMS::data_id1:     //ne    
      new_ne(v,f0,c,MPIc,answer);
      break;
    case VISIT_PARAMS::data_id2:     //ni    
      new_Constant(&Initial_Conditions::ni,c,MPIc,answer,"ni");
      break;
    case VISIT_PARAMS::data_id3:     //Te  
      new_Te(v,f0,c,MPIc,answer);
      break;
    case VISIT_PARAMS::data_id4:     //Bz    
      IMPACT_Dim x3(3);
      new_B(v,B,c,MPIc,answer,&x3);
      break;
    }
}
void IMPACT_Normalize_Vector(struct IMPACT_Moment * answerx,
			     struct IMPACT_Moment * answery,
			     struct IMPACT_Moment * answerz,
			     IMPACT_Config *c)
{
  double val=0.0,max=0.0;
  for (int j=1;j<=c->Ny();++j)

    for (int i=1;i<=c->Nx();++i)
      {
	val=answerx->values.get(i,j)*answerx->values.get(i,j)
	  +answery->values.get(i,j)*answery->values.get(i,j)
	  +answerz->values.get(i,j)*answerz->values.get(i,j);
	if (val>max) max=val;
      }
  if (max==0.0) max=1;
  int ii=1;
  max=1/sqrt(max)*c->dx(&ii);
  for (int j=1;j<=c->Ny();++j)
    for (int i=1;i<=c->Nx();++i){
      answerx->values.set(i,j,answerx->values.get(i,j)*max);
      answery->values.set(i,j,answery->values.get(i,j)*max);
      answerz->values.set(i,j,answerz->values.get(i,j)*max);
    }
}
//Choose the correct vector
void IMPACT_Select_Vector(struct IMPACT_Moment * answerx,
			  struct IMPACT_Moment * answery,
			  struct IMPACT_Moment * answerz,
			  int choice,
			  IMPACT_Config *c, 
			  IMPACT_MPI_Config *MPIc, IMPACT_ParVec *v,
			  IMPACT_Var *f0,IMPACT_Var *f1,IMPACT_Var *f2,
			  IMPACT_Var *f3,
			  IMPACT_Var *E,IMPACT_Var *B)
{
  //answer->values.ChangeSize(c->Nx(),c->Ny());
  struct IMPACT_Moment* answerpointer[3]={answerx,answery,answerz};
  for (int t=0;t<3;++t){
    answerpointer[t]->values.ChangeSize(c->Nx(),c->Ny());
    for (int j=1;j<=c->Ny();++j)
      for (int i=1;i<=c->Nx();++i)
	answerpointer[t]->values.set(i,j,0.0);}
  switch (choice)
    {
      /*case VISIT_PARAMS::vec_id0:
	
	break;*/
    case VISIT_PARAMS::vec_id1:     //E    
      if(c->NE()>0)
	for (IMPACT_Dim x1=1;x1<=c->NE();++x1){
	  MPI_Barrier(MPI_COMM_WORLD);
	  new_E(v,E,c,MPIc,answerpointer[x1.get()-1],&x1);}
      break;
    case VISIT_PARAMS::vec_id2:     //B    
      if(c->NB()>0)
	for (IMPACT_Dim x1=3;x1>3-c->NB();--x1)
	  {
	    MPI_Barrier(MPI_COMM_WORLD);
	    new_B(v,B,c,MPIc,answerpointer[x1.get()-1],&x1);}
      break;
    case VISIT_PARAMS::vec_id3:     //j    
      for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
	  {	
	    MPI_Barrier(MPI_COMM_WORLD);
	    new_je(v,f1,c,MPIc,answerpointer[x1.get()-1],&x1);}
      break;
    case VISIT_PARAMS::vec_id4:     //q    
      for (IMPACT_Dim x1=1;x1<=c->Nf1();++x1)
	  {	
	    MPI_Barrier(MPI_COMM_WORLD);
	    new_qT(v,f1,c,MPIc,answerpointer[x1.get()-1],&x1);}
      break;
      }
}
