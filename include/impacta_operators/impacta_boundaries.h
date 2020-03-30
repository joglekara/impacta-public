/*
 * *********************************************************
 * IMPACTA 
 * 
 * // BOundary condition code 
 * Boundary conditions added - open boundaries
 * - vector gradients go to 0 at bounds
 * - not vectors (as in  reflecting)
 * - all others as reflecting
 * 
 * - to account for new boundary conditions, rearrangemet of 
 * namespace c.f. earlier versions
 * 
 * NB NB NB !!!! int Open_Bound_x() & int Open_Bound_y()
 * in IMPACT_Operatros
 * 
 * 1/7/08 - New boundary condition for B field - reflective
 * but with j crossing boundary (& E), B doesn't flip over bound
 * bound is _b_b
 * 
 * Also - fixed boundaries changed to the lazy way 
 * - f0 not updated on boundary cells.
 * 
 * july 2011 - AGRT - added functino pointers to boundary conditions!
 */


void IMPACT_fixed_bound(IMPACT_stencil *stencil,int index,int *i)
{
  /* IMPACT_Boundaries::fixed_f0_grid[index][*i]=stencil->get(index);
   *     stencil->set(index,0.0);*/
}
int Open_Bound_x();
int Open_Bound_y();

//Periodic for both axes
namespace BC_x_p_y_p
{
  char boundaries[2][11]={"periodic","periodic"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  //no fixed for periodic!
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    if (*jplus>Ny) *jplus=1;
    if (*jminus<1) *jminus=Ny;
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_p_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_p_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_p_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_p_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);  
    if (*iplus>Nx) *iplus=1; //for periodic bounds
    if (*iminus<1) *iminus=Nx;
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_p_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
    if (*iplus>Nx) *iplus=1; //for periodic bounds
    if (*iminus<1) *iminus=Nx;
  }
}
//**************************************************************************
// reflective in y
namespace BC_x_p_y_rg
{
  // char boundaries[2][11]={"periodic","reflecting"};
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f0[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]); 
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f0[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
    }
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_ni[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]); 
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_ni[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
    }
    if (*iplus>Nx) *iplus=1; //for periodic bounds
    if (*iminus<1) *iminus=Nx;
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f1_E[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(3);
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f1_E[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(4);
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_Ci[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(3);
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f1_E[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(4);
    }
    if (*iplus>Nx) *iplus=1; //for periodic bounds
    if (*iminus<1) *iminus=Nx;
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    if (x1->get()==x2->get())
    {
      if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f2[3]) *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f2[2]) *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
    }
    else
    {
      if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f2[3]) {
	*jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
	if(!Open_Bound_y())
	  if (x1->get()==2||x2->get()==2) stencil->invert(3);}
	  if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f2[2]) {
	    *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
	    if(!Open_Bound_y())
	      if (x1->get()==2||x2->get()==2) stencil->invert(4);}
    }
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_B[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if(!Open_Bound_y())
	if (x1->get()==3||x1->get()==1) stencil->invert(3);
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_B[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
      if(!Open_Bound_y())
	if (x1->get()==1||x1->get()==3) stencil->invert(4);
    }
  }
}
//**************************************************************************

//reflective in x
namespace BC_x_rg_y_p
{
  //  char boundaries[2][11]={"reflecting","periodic"};
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f0[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f0[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
    }
    if (*jplus>Ny) *jplus=1;
    if (*jminus<1) *jminus=Ny;
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_ni[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_ni[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
    }
    if (*jplus>Ny) *jplus=1;
    if (*jminus<1) *jminus=Ny;
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f1_E[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(1);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f1_E[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(2);
    }
    
    if (*jplus>Ny) *jplus=1;
    if (*jminus<1) *jminus=Ny;
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_Ci[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(1);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_Ci[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(2);
    }
    
    if (*jplus>Ny) *jplus=1;
    if (*jminus<1) *jminus=Ny;
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    if (x1->get()==x2->get())
    {
      if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f2[1]) *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f2[0]) *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if (*jplus>Ny) *jplus=1;
      if (*jminus<1) *jminus=Ny;
    }
    else
    {
      if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f2[1]) {
	*iplus=(Nx-IMPACT_Boundaries::fixed_any[1]); 
	if(!Open_Bound_x())
	  if (x1->get()==1||x2->get()==1) stencil->invert(1);}
	  if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f2[0]) {
	    *iminus=(1+IMPACT_Boundaries::fixed_any[0]); 
	    if(!Open_Bound_x())
	      if (x1->get()==1||x2->get()==1) stencil->invert(2);}
	      if (*jplus>Ny) *jplus=1;
	      if (*jminus<1) *jminus=Ny;
    }
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_B[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if(!Open_Bound_x())
	if (x1->get()==2||x1->get()==3) stencil->invert(1);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_B[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if(!Open_Bound_x())
	if (x1->get()==2||x1->get()==3) stencil->invert(2);
    }
    if (*jplus>Ny) *jplus=1;
    if (*jminus<1) *jminus=Ny;
  }
}
//**************************************************************************
//reflective in both
namespace BC_x_rg_y_rg
{
  // char boundaries[2][11]={"reflecting","reflecting"};
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f0[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f0[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
    }
    BC_x_p_y_rg::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_ni[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_ni[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
    }
    BC_x_p_y_rg::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f1_E[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(1);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f1_E[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(2);
    }
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f1_E[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(3);
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f1_E[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(4);
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_Ci[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(1);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_Ci[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if(!Open_Bound_x())
	if (x1->get()==1) stencil->invert(2);
    }
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_Ci[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(3);
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_Ci[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
      if(!Open_Bound_y())
	if (x1->get()==2) stencil->invert(4);
    }
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    if (x1->get()==x2->get())
    {
      if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f2[1]) *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f2[0]) *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f2[3]) *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f2[2]) *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
    }
    else
    {
      if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_f2[1]) {
	*iplus=(Nx-IMPACT_Boundaries::fixed_any[1]); 
	if(!Open_Bound_x())
	  if (x1->get()==1||x2->get()==1) stencil->invert(1);}
	  if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_f2[0]) {
	    *iminus=(1+IMPACT_Boundaries::fixed_any[0]); 
	    if(!Open_Bound_x())
	      if (x1->get()==1||x2->get()==1) stencil->invert(2);}
	      if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_f2[3]) {
		*jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
		if(!Open_Bound_y())
		  if (x1->get()==2||x2->get()==2) stencil->invert(3);}
		  if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_f2[2]) {
		    *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
		    if(!Open_Bound_y())
		      if (x1->get()==2||x2->get()==2) stencil->invert(4);}
    }
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*iplus>(Nx-IMPACT_Boundaries::fixed_any[1])&&!IMPACT_Boundaries::fix_B[1]) 
    {
      *iplus=(Nx-IMPACT_Boundaries::fixed_any[1]);
      if(!Open_Bound_x())
	if (x1->get()==2||x1->get()==3) stencil->invert(1);
    }
    if (*iminus<(1+IMPACT_Boundaries::fixed_any[0])&&!IMPACT_Boundaries::fix_B[0]) 
    {
      *iminus=(1+IMPACT_Boundaries::fixed_any[0]);
      if(!Open_Bound_x())
	if (x1->get()==2||x1->get()==3) stencil->invert(2);
    }
    
    if (*jplus>(Ny-IMPACT_Boundaries::fixed_any[3])&&!IMPACT_Boundaries::fix_B[3]) 
    {
      *jplus=(Ny-IMPACT_Boundaries::fixed_any[3]);
      if(!Open_Bound_y())
	if (x1->get()==3||x1->get()==1) stencil->invert(3);
    }
    if (*jminus<(1+IMPACT_Boundaries::fixed_any[2])&&!IMPACT_Boundaries::fix_B[2]) 
    {
      *jminus=(1+IMPACT_Boundaries::fixed_any[2]);
      if(!Open_Bound_y())
	if (x1->get()==1||x1->get()==3) stencil->invert(4);
    }
  }
}
namespace BC_x_r_y_r
{
  char boundaries[2][11]={"reflecting","reflecting"};
  using namespace BC_x_rg_y_rg;
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
}
namespace BC_x_r_y_p
{
  char boundaries[2][11]={"reflecting","periodic"};
  using namespace BC_x_rg_y_p;
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
}
namespace BC_x_p_y_r
{
  char boundaries[2][11]={"periodic","reflecting"};
  using namespace BC_x_p_y_rg;
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
}
//---------------------------------------------------------------------------

//Reflective with magnetic field continuous across boundary
namespace BC_x_b_y_b
{
  char boundaries[2][11]={"reflect_Bz","reflect_Bz"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_r::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_r::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_r::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_r_y_r::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}
//---------------------------------------------------------------------------

//Reflective with magnetic field continuous across boundary
namespace BC_x_b_y_p
{
  char boundaries[2][11]={"reflect_Bz","periodic"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_p::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_p::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_p::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_p::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_r_y_p::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}

//Reflective with magnetic field continuous across right boundary boundary but going to 0 at the left
namespace BC_xl_r_xr_b_y_p
{
  char boundaries[2][11]={"reflect_Bz","periodic"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_p::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
      //left side
    if (*iminus<=1)
    {
      BC_x_r_y_p::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
    else
    {
      //right side
      BC_x_r_y_p::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);    
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_p::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
      //left side
    if (*iminus<=1)
    { 
      BC_x_r_y_p::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);    
    }
    else
    {  
      BC_x_r_y_p::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_r_y_p::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}
namespace BC_xl_r_xr_b_y_r
{
  char boundaries[2][11]={"reflect_Bz","periodic"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_p::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    //left side
    if (*iminus<=1)
    {
      BC_x_r_y_p::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
    else
    {
      //right side
      BC_x_r_y_p::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);    
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_p::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
      //left side
    if (*iminus<=1)
    { 
      BC_x_r_y_p::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);    
    }
    else
    {  
      BC_x_r_y_p::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_r_y_p::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}
namespace BC_xl_r_xr_p_y_r
{
  char boundaries[2][11]={"reflect_Bz","reflecting"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
    
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    //left side
    if (*iminus<=1)
    {
      BC_x_r_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
    }
    else
    {
      //right side
      BC_x_p_y_p::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
    }
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
     //left side
    if (*iminus<=1)
    {
      BC_x_r_y_r::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
    }
    else
    {
      //right side
      BC_x_p_y_r::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
    }
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
      //left side
    if (*iminus<=1)
    {
      BC_x_r_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
    else
    {
      //right side
      BC_x_p_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
     //left side
    if (*iminus<=1)
    {
      BC_x_r_y_r::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
    else
    {
      //right side
      BC_x_p_y_r::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    
      //left side
    if (*iminus<=1)
    { 
      BC_x_r_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);    
    }
    else
    {  
      BC_x_p_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);    
    }
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
     //left side
    if (*iminus<=1)
    {
      BC_x_r_y_r::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
    }
    else
    {
      //right side
      BC_x_p_y_r::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
    }
  }
}
//Reflective with magnetic field continuous across boundary
namespace BC_x_p_y_b
{
  char boundaries[2][11]={"periodic","reflect_Bz"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_p_y_r::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_p_y_r::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_p_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  } 
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_p_y_r::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_p_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_p_y_r::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}
//---------------------------------------------------------------------------

//Reflective with magnetic field continuous across boundary in one, and normal reflective the other way
namespace BC_x_b_y_r
{
  char boundaries[2][11]={"reflect_Bz","reflecting"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_r::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_r::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*iminus<=1||*iplus>=Nx) {
      BC_x_r_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    } else {
      BC_x_r_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_r::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*iminus<=1||*iplus>=Nx) {
      BC_x_r_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    } else {
      BC_x_r_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    } 
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_r_y_r::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}
namespace BC_x_r_y_b
{
  char boundaries[2][11]={"reflecting","reflect_Bz"};
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=0;
  
  inline void IMPACT_f0_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_r::IMPACT_f0_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline void IMPACT_ni_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus)
  {
    BC_x_r_y_r::IMPACT_ni_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus);
  }
  inline  void IMPACT_f1_E_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*jminus<=1||*jplus>=Ny) {
      BC_x_r_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    } else {
      BC_x_r_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    }
  }
  inline  void IMPACT_Ci_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    BC_x_r_y_r::IMPACT_Ci_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
  }
  inline  void IMPACT_B_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1)
  {
    if (*jminus<=1||*jplus>=Ny) {
      BC_x_r_y_r::IMPACT_f1_E_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    } else {
      BC_x_r_y_r::IMPACT_B_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1);
    } 
  }
  inline  void IMPACT_f2_bound(IMPACT_stencil *stencil,int Nx, int Ny, int *iplus,int *iminus, int *jplus, int *jminus, IMPACT_Dim *x1, IMPACT_Dim *x2)
  {
    BC_x_r_y_r::IMPACT_f2_bound(stencil,Nx,Ny,iplus,iminus,jplus,jminus,x1,x2);
  }
}
//---------------------------------------------------------------------------
// Open boundary conditions

namespace BC_x_o_y_o
{
  char boundaries[2][11]={"open","open"};
  using namespace BC_x_rg_y_rg;
  const int Boundaries_open_x=1;
  const int Boundaries_open_y=1;
}
namespace BC_x_r_y_o
{
  char boundaries[2][11]={"reflecting","open"};
  using namespace BC_x_rg_y_rg;
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=1;
}
namespace BC_x_o_y_r
{
  char boundaries[2][11]={"open","reflecting"};
  using namespace BC_x_rg_y_rg;
  const int Boundaries_open_x=1;
  const int Boundaries_open_y=0;
}
namespace BC_x_o_y_p
{
  char boundaries[2][11]={"open","periodic"};
  using namespace BC_x_rg_y_p;
  const int Boundaries_open_x=1;
  const int Boundaries_open_y=0;
}
namespace BC_x_p_y_o
{
  char boundaries[2][11]={"periodic","open"};
  using namespace BC_x_p_y_rg;
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=1;
}

namespace BC_x_b_y_o
{
  char boundaries[2][11]={"reflect_Bz","open"};
  using namespace BC_x_b_y_b;
  const int Boundaries_open_x=0;
  const int Boundaries_open_y=1;
}
namespace BC_x_o_y_b
{
  char boundaries[2][11]={"open","reflect_Bz"};
  using namespace BC_x_b_y_b;
  const int Boundaries_open_x=1;
  const int Boundaries_open_y=0;
}


void IMPACTA_initiate_bounds()
{
  
  switch (IMPACT_Boundaries::boundary_type)
  {
    case 0:
      IMPACT_f0_bound=&BC_x_p_y_p::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_p_y_p::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_p_y_p::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_p_y_p::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_p_y_p::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_p_y_p::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_p_y_p::boundaries[0][ii];
	boundaries[1][ii]=BC_x_p_y_p::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_p_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_p_y_p::Boundaries_open_y;
      break;
    case 1:
      IMPACT_f0_bound=&BC_x_r_y_p::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_r_y_p::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_r_y_p::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_r_y_p::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_r_y_p::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_r_y_p::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_r_y_p::boundaries[0][ii];
	boundaries[1][ii]=BC_x_r_y_p::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_r_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_r_y_p::Boundaries_open_y;
      break;
    case 2:
      IMPACT_f0_bound=&BC_x_b_y_p::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_b_y_p::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_b_y_p::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_b_y_p::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_b_y_p::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_b_y_p::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_b_y_p::boundaries[0][ii];
	boundaries[1][ii]=BC_x_b_y_p::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_p::Boundaries_open_y;
      break;
    case 3:
      IMPACT_f0_bound=&BC_x_o_y_p::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_o_y_p::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_o_y_p::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_o_y_p::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_o_y_p::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_o_y_p::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_o_y_p::boundaries[0][ii];
	boundaries[1][ii]=BC_x_o_y_p::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_o_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_o_y_p::Boundaries_open_y;
      break;
    case 4:
      IMPACT_f0_bound=&BC_xl_r_xr_b_y_p::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_xl_r_xr_b_y_p::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_xl_r_xr_b_y_p::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_xl_r_xr_b_y_p::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_xl_r_xr_b_y_p::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_xl_r_xr_b_y_p::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_b_y_p::boundaries[0][ii];
	boundaries[1][ii]=BC_x_b_y_p::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_p::Boundaries_open_y;
      break;
    case 10:
      IMPACT_f0_bound=&BC_x_p_y_r::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_p_y_r::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_p_y_r::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_p_y_r::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_p_y_r::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_p_y_r::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_p_y_r::boundaries[0][ii];
	boundaries[1][ii]=BC_x_p_y_r::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_p_y_r::Boundaries_open_x;
      Boundaries_open_y=BC_x_p_y_r::Boundaries_open_y;
      break;
    case 11:
      IMPACT_f0_bound=&BC_x_r_y_r::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_r_y_r::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_r_y_r::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_r_y_r::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_r_y_r::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_r_y_r::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_r_y_r::boundaries[0][ii];
	boundaries[1][ii]=BC_x_r_y_r::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_r_y_r::Boundaries_open_x;
      Boundaries_open_y=BC_x_r_y_r::Boundaries_open_y;
      break;
    case 12:
      IMPACT_f0_bound=&BC_x_b_y_r::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_b_y_r::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_b_y_r::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_b_y_r::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_b_y_r::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_b_y_r::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_b_y_r::boundaries[0][ii];
	boundaries[1][ii]=BC_x_b_y_r::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_r::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_r::Boundaries_open_y;
      break;
    case 13:
      IMPACT_f0_bound=&BC_x_o_y_r::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_o_y_r::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_o_y_r::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_o_y_r::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_o_y_r::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_o_y_r::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_o_y_r::boundaries[0][ii];
	boundaries[1][ii]=BC_x_o_y_r::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_o_y_r::Boundaries_open_x;
      Boundaries_open_y=BC_x_o_y_r::Boundaries_open_y;
      break;
  case 14:
      IMPACT_f0_bound=&BC_xl_r_xr_b_y_r::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_xl_r_xr_b_y_r::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_xl_r_xr_b_y_r::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_xl_r_xr_b_y_r::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_xl_r_xr_b_y_r::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_xl_r_xr_b_y_r::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
  boundaries[0][ii]=BC_x_b_y_r::boundaries[0][ii];
  boundaries[1][ii]=BC_x_b_y_r::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_p::Boundaries_open_y;
      break;
   case 15:
      IMPACT_f0_bound=&BC_xl_r_xr_p_y_r::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_xl_r_xr_p_y_r::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_xl_r_xr_p_y_r::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_xl_r_xr_p_y_r::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_xl_r_xr_p_y_r::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_xl_r_xr_p_y_r::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
  boundaries[0][ii]=BC_x_p_y_r::boundaries[0][ii];
  boundaries[1][ii]=BC_x_p_y_r::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_p::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_p::Boundaries_open_y;
      break;
    case 20:
      IMPACT_f0_bound=&BC_x_p_y_b::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_p_y_b::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_p_y_b::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_p_y_b::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_p_y_b::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_p_y_b::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_p_y_b::boundaries[0][ii];
	boundaries[1][ii]=BC_x_p_y_b::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_p_y_b::Boundaries_open_x;
      Boundaries_open_y=BC_x_p_y_b::Boundaries_open_y;
      break;
    case 21:
      IMPACT_f0_bound=&BC_x_r_y_b::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_r_y_b::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_r_y_b::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_r_y_b::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_r_y_b::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_r_y_b::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_r_y_b::boundaries[0][ii];
	boundaries[1][ii]=BC_x_r_y_b::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_r_y_b::Boundaries_open_x;
      Boundaries_open_y=BC_x_r_y_b::Boundaries_open_y;
      break;
    case 22:
      IMPACT_f0_bound=&BC_x_b_y_b::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_b_y_b::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_b_y_b::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_b_y_b::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_b_y_b::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_b_y_b::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_b_y_b::boundaries[0][ii];
	boundaries[1][ii]=BC_x_b_y_b::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_b::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_b::Boundaries_open_y;
      break;
    case 23:
      IMPACT_f0_bound=&BC_x_o_y_b::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_o_y_b::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_o_y_b::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_o_y_b::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_o_y_b::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_o_y_b::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_o_y_b::boundaries[0][ii];
	boundaries[1][ii]=BC_x_o_y_b::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_o_y_b::Boundaries_open_x;
      Boundaries_open_y=BC_x_o_y_b::Boundaries_open_y;
      break;
    case 30:
      IMPACT_f0_bound=&BC_x_p_y_o::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_p_y_o::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_p_y_o::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_p_y_o::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_p_y_o::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_p_y_o::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_p_y_o::boundaries[0][ii];
	boundaries[1][ii]=BC_x_p_y_o::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_p_y_o::Boundaries_open_x;
      Boundaries_open_y=BC_x_p_y_o::Boundaries_open_y;
      break;
    case 31:
      IMPACT_f0_bound=&BC_x_r_y_o::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_r_y_o::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_r_y_o::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_r_y_o::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_r_y_o::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_r_y_o::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_r_y_o::boundaries[0][ii];
	boundaries[1][ii]=BC_x_r_y_o::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_r_y_o::Boundaries_open_x;
      Boundaries_open_y=BC_x_r_y_o::Boundaries_open_y;
      break;
    case 32:
      IMPACT_f0_bound=&BC_x_b_y_o::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_b_y_o::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_b_y_o::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_b_y_o::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_b_y_o::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_b_y_o::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_b_y_o::boundaries[0][ii];
	boundaries[1][ii]=BC_x_b_y_o::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_b_y_o::Boundaries_open_x;
      Boundaries_open_y=BC_x_b_y_o::Boundaries_open_y;
      break;
    case 33:
      IMPACT_f0_bound=&BC_x_o_y_o::IMPACT_f0_bound;
      IMPACT_f1_E_bound=&BC_x_o_y_o::IMPACT_f1_E_bound;
      IMPACT_f2_bound=&BC_x_o_y_o::IMPACT_f2_bound;
      IMPACT_B_bound=&BC_x_o_y_o::IMPACT_B_bound;
      IMPACT_Ci_bound=&BC_x_o_y_o::IMPACT_Ci_bound;
      IMPACT_ni_bound=&BC_x_o_y_o::IMPACT_ni_bound;
      for (int ii=0;ii<11;++ii) {
	boundaries[0][ii]=BC_x_o_y_o::boundaries[0][ii];
	boundaries[1][ii]=BC_x_o_y_o::boundaries[1][ii];
      }
      Boundaries_open_x=BC_x_o_y_o::Boundaries_open_x;
      Boundaries_open_y=BC_x_o_y_o::Boundaries_open_y;
      break;
    default:
      std::cout<<"Boundary conditions not valid\n";
      exit(0);
  } 
  
}
