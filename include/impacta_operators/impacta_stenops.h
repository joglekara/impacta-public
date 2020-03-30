/*
**********************************************************
IMPACT Version 2.3
Stencil operator "box" - contains series of operator
stencils to be passed from function to funciton

Version 2.7
1
AGRT

21/2/07

16/4/07 - since j operator components are actually all
the same(!) I have eliminated two of them,

25/1/2010 AGRT
I am including a stencil which introduces some numerical
diffusion. See IMPACTA_Switch Elements for details

**********************************************************

*/


#ifndef INC_IMPACTA_STENOPS_H
#define INC_IMPACTA_STENOPS_H

//system

//user

using namespace globalconsts;


// Center differenced operators
namespace differentials_x_c_y_c
{
  char xdiff[10]="center";
  char ydiff[10]="center";

  

  class IMPACT_StenOps // stencil operators
  {
   private:
    IMPACT_stencil * null_;
    IMPACT_stencil * ddx_; //i.e. d/d|x| ^x (x hat)
    IMPACT_stencil * ddy_; //ditto for y
    IMPACT_stencil * ddxfw_; //i.e. d/d|x| ^x (x hat) forward differenced 
    IMPACT_stencil * ddyfw_; //ditto for y
    IMPACT_stencil * ddxbw_; //i.e. d/d|x| ^x (x hat) backward differenced 
    IMPACT_stencil * ddybw_; //ditto for y
    IMPACT_stencil ** div_; // divergence stencil
    IMPACT_stencil *d2dx_;
    IMPACT_stencil *d2dy_;
    IMPACT_vint je_; // current stencil for E equation
   public:
    IMPACT_StenOps(IMPACT_Config *c);
    ~IMPACT_StenOps();
    //Access methods for differential stencils
    IMPACT_stencil *ddx(int *i);
    IMPACT_stencil *ddy(int *j);
    IMPACT_stencil *div(int *i,int *j);
    IMPACT_stencil *ddxi(int *i,int *j, IMPACT_Dim *x1);
    IMPACT_stencil *ddxi_uw(int *i,int *j, IMPACT_Dim *x1, bool uw);
    IMPACT_stencil *d2dx2i(int *i, int *j, IMPACT_Dim *x1);
    IMPACT_stencil Laplacian(int *i, int *j);
    IMPACT_vint * je(IMPACT_Dim *x1);
    IMPACT_stencil * null();
  };

  IMPACT_StenOps::IMPACT_StenOps(IMPACT_Config *c)
  {
    //This constructor creates a templateobject with
    // variaous cartesian differential and integral operators.
    ddx_=new IMPACT_stencil[c->Nx()+2];
    ddy_=new IMPACT_stencil[c->Ny()+2];
	ddxfw_=new IMPACT_stencil[c->Nx()+2];
    ddyfw_=new IMPACT_stencil[c->Ny()+2];
	ddxbw_=new IMPACT_stencil[c->Nx()+2];
    ddybw_=new IMPACT_stencil[c->Ny()+2];
    d2dx_=new IMPACT_stencil[c->Nx()+2];
    d2dy_=new IMPACT_stencil[c->Ny()+2];
    null_ = new IMPACT_stencil(0.0,0.0,0.0,0.0,0.0);
    int iplus,iminus,jplus,jminus;
    for (int i=1;i<c->Nx()+1;++i) // (1/(dx[i+1]+dx[i-1]) and its negative in i+1 and i-1
    {
      iplus=i+1;
      iminus=i-1;
      ddx_[i]=IMPACT_stencil(0.0,1.0/(c->dx(&iminus)+c->dx(&iplus)),-1.0/(c->dx(&iminus)+c->dx(&iplus)),0.0,0.0);
      ddxfw_[i]=IMPACT_stencil(-1.0/(c->dx(&i)),1.0/(c->dx(&i)),0.0,0.0,0.0);
      ddxbw_[i]=IMPACT_stencil(1.0/(c->dx(&i)),0.0,-1.0/(c->dx(&i)),0.0,0.0);
      d2dx_[i]=IMPACT_stencil(-2.0/c->dx(&i)/c->dx(&i),1.0/c->dx(&i)/c->dx(&i),1.0/c->dx(&i)/c->dx(&i),0.0,0.0);
    }
    for (int j=1;j<c->Ny()+1;++j)
    {
      jplus=j+1;
      jminus=j-1;
      ddy_[j]=IMPACT_stencil(0.0,0.0,0.0,1.0/(c->dy(&jminus)+c->dy(&jplus)),-1.0/(c->dy(&jminus)+c->dy(&jplus)));
	  ddyfw_[j]=IMPACT_stencil(-1.0/(c->dy(&j)),0.0,0.0,1.0/(c->dy(&j)),0.0);
	  ddybw_[j]=IMPACT_stencil(1.0/(c->dy(&j)),0.0,0.0,0.0,-1.0/(c->dy(&j)));
      d2dy_[j]=IMPACT_stencil(-2.0/c->dy(&j)/c->dy(&j),0.0,0.0,1.0/c->dy(&j)/c->dy(&j),1.0/c->dy(&j)/c->dy(&j));
    }
    div_=new IMPACT_stencil * [c->Nx()+2];

    div_[0]=new IMPACT_stencil [c->Ny()+2];
    div_[c->Nx()+1]=new IMPACT_stencil [c->Ny()+2];

    for (int i=1;i<c->Nx()+1;++i)
    {
      div_[i]=new IMPACT_stencil [c->Ny()+2];
      for (int j=1;j<c->Ny()+1;++j)
      {
        div_[i][j]=(ddx_[i]+ddy_[j]);
      }
    }

    /*
    This next stencil is for the current density (int f1 v^3 dv)
    but actually it is j/epsilon0 as it is in the equation:
    E_n+1+dtj_n+1/epsilon0-c^2 dt curlB_n+1 = E_n
     */
    je_.resize(c);
    double val;
    for (int k=1;k<=c->Nv();++k)
    {
      val = -fourthirdspi*c->v3(&k)*c->dv(&k)*c_L*c_L
      *equation_switches::jimp_in_E_equation/deltasquared;
      // /epsilon0*e_charge*c->dt(); 
      // v^3dv/epsilon0*-e *dt
      je_.set(k,val);
    }
  }

  IMPACT_StenOps::~IMPACT_StenOps()
  {
    delete[] ddx_;
    delete[] ddy_;
    delete[] ddxfw_;
    delete[] ddyfw_;
    delete[] ddxbw_;
    delete[] ddybw_;
    delete[] div_;
    delete null_;
    delete[] d2dx_;
    delete[] d2dy_;
   }
  IMPACT_stencil * IMPACT_StenOps::ddx(int *i)
  {
    return &ddx_[*i];
  }
  IMPACT_stencil * IMPACT_StenOps::ddy(int *j)
  {
    return &ddy_[*j];
  }
  IMPACT_stencil * IMPACT_StenOps::div(int *i,int *j)
  {
    return &div_[*i][*j];
  }
  IMPACT_stencil * IMPACT_StenOps::d2dx2i(int *i, int *j, IMPACT_Dim *x1)
  {
    IMPACT_stencil * temp;
    temp=null_;
    switch (x1->get())
      {
      case 1:
        temp = &d2dx_[*i];
        break;
      case 2:
         temp = &d2dy_[*j];
        break;
      case 3:
        temp = null_;
      }
    return temp;
  }
  IMPACT_stencil IMPACT_StenOps::Laplacian(int *i, int *j)
  {
    IMPACT_stencil temp;
    temp=d2dx_[*i]+d2dy_[*j];
    return temp;
  }
  IMPACT_stencil * IMPACT_StenOps::ddxi(int *i,int *j, IMPACT_Dim *x1)
  {
    IMPACT_stencil * temp;
    temp=null_;
    switch (x1->get())
      {
      case 1:
        temp = &ddx_[*i];
        break;
      case 2:
        temp = &ddy_[*j];
        break;
      case 3:
        temp = null_;
      }
    return temp;
  }
  IMPACT_stencil * IMPACT_StenOps::ddxi_uw(int *i,int *j, IMPACT_Dim *x1, bool uw)
  {
    IMPACT_stencil * temp;
    temp=null_;
    switch (x1->get())
      {
      case 1:
        if (uw)			temp = &ddxfw_[*i];
        else			temp = &ddxbw_[*i];
//                temp = &ddx_[*i];
        break;
      case 2:
        if (uw)			temp = &ddyfw_[*j];
        else			temp = &ddybw_[*j]; 
//                temp = &ddy_[*j];
        break;
      case 3:
        temp = null_;
      }
    return temp;
  }
  IMPACT_vint * IMPACT_StenOps::je(IMPACT_Dim *x1)
  {
    return &je_;
  } 
  IMPACT_stencil * IMPACT_StenOps::null()
  {
    return null_;
  }
  // Lax stencil - averaging over adjacent cells
  IMPACT_stencil LAX(1.0,0.0,0.0,0.0,0.0);

} // end of namespace





#endif /* INC_IMPACTA_STENOPS_H  */
