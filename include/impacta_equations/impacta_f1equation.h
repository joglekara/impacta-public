/*
**********************************************************
 f1 equation IMPACT code for inner and outer cells

Version 1.2
1
AGRT
20/2/07
27/2 - Need to put in Bxf1 term!!! done.

15/04/07 - Changed div for grad
AGRT 2022 - now on github, updates there.

**********************************************************

*/

// code for inner cells of f0 equation - i.e. not at boundaries
inline void f1equation_inner(IMPACT_Config *c, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *E, IMPACT_Var *B, IMPACT_Var *f1, IMPACT_Var *f2, int *i, int *j, int *k, IMPACT_Dim *x1)
{
    /*
  **********************************************************
  INSERT f1 Equation terms CODE HERE
  **********************************************************
     */

    f1->set(vtemp, equation_switches::e_inert_on * c->idt(), i, j, k, x1);
    // set f1 n+1 element to 1.0

    // vgrad f0
    double v = c->v(k) * equation_switches::inf1_vgradf0_on;
    IMPACT_stencil temp_sten = (*O->ddxi(i, j, x1)) * v;

    f0->insert(&temp_sten, vtemp, i, j, k);

    // Gradient in z term
    if (x1->get() == 3)
    {
        f0->inc(vtemp, v * IMPACT_Heating::Dnbydzgrid.get(i, j), i, j, k);
        f0->inc(vtemp, (v * v / Initial_Conditions::Te.get(i, j) - 1.5) * v * IMPACT_Heating::DTbydzgrid.get(i, j), i, j, k);
    }

    //-Edf0/dv
    if (!if_use_hybrid_nonlinear)
    {
        if (c->NE() > 0)
            E->set(vtemp, -equation_switches::inf1_Edf0dv_on * f0->ddv(vlagged, c, i, j, k), i, j, x1);
    }
    else
    {
        if (c->NE() > 0)
        {
            E->set(vtemp, -0.5 * equation_switches::inf1_Edf0dv_on * f0->ddv(vlagged, c, i, j, k), i, j, x1);
            IMPACT_Vel_Sten ddv = O->ddv(c, k);
            ddv *(-0.5 * equation_switches::inf1_Edf0dv_on *E->get(vlagged, i, j, x1));
            f0->insert(&ddv, vtemp, i, j, k);
        }
    }
    // double erfv=erf(v);
    // double expmv2oversqrtpi = exp(-1.0*c->v2(k))/sqrt(globalconsts::pi);

    // double ei_cols=(equation_switches::Cei_on*Initial_Conditions::Z.get(i,j)*Initial_Conditions::Z.get(i,j)*Initial_Conditions::ni.get(i,j)
    //    -(erfv*(0.5/c->v2(k)-1.0)-expmv2oversqrtpi/c->v(k))
    //    )/c->v3(k)-2.0*f0->get(vlagged,i,j,k);

    double ei_cols = equation_switches::Cei_on * Initial_Conditions::Z.get(i, j) * Initial_Conditions::Z.get(i, j) * Initial_Conditions::ni.get(i, j) / c->v3(k) / ((Initial_Conditions::Z.get(i, j) / globalconsts::oneover_atomic_Z + 0.236) / (1.0 + 0.236 * Initial_Conditions::Z.get(i, j) / globalconsts::oneover_atomic_Z) * .236);

    //      +
    // 0*(expmv2oversqrtpi-globalconsts::fourthirds*
    //  (expmv2oversqrtpi*c->v2(k)+c->v3(k)*(1.0-erfv))+erfv*(c->v(k)-0.5/c->v(k)))
    // *((Initial_Conditions::Z.get(i,j)/globalconsts::oneover_atomic_Z+0.236)/(1.0+0.236*Initial_Conditions::Z.get(i,j)/globalconsts::oneover_atomic_Z)
    // *.236)

    // collision term - only ion collisions first iteration
    f1->inc(vtemp, ei_cols, i, j, k, x1);

    // B x f1
    double multiplier = -equation_switches::inf1_f1xB_on;

    IMPACT_Dim x2, x3;
    GetOrthogonal(x1, &x2, &x3);
    if (!if_use_hybrid_nonlinear)
    {
        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            f1->inc(vtemp, B->get(vlagged, i, j, &x2) * multiplier * Bimp, i, j, k, &x3);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            f1->inc(vtemp, -(B->get(vlagged, i, j, &x3) * multiplier * Bimp), i, j, k, &x2);

        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            B->inc(vtemp, f1->get(vlagged, i, j, k, &x3) * multiplier * (1 - Bimp), i, j, &x2);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            B->inc(vtemp, -f1->get(vlagged, i, j, k, &x2) * multiplier * (1 - Bimp), i, j, &x3);
    }
    else
    {
        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            f1->inc(vtemp, B->get(vlagged, i, j, &x2) * multiplier * 0.5, i, j, k, &x3);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            f1->inc(vtemp, -(B->get(vlagged, i, j, &x3) * multiplier * 0.5), i, j, k, &x2);

        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            B->inc(vtemp, f1->get(vlagged, i, j, k, &x3) * multiplier * 0.5, i, j, &x2);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            B->inc(vtemp, -f1->get(vlagged, i, j, k, &x2) * multiplier * 0.5, i, j, &x3);
    }
    // f2 parts
    // a_j d/dvv^3f2_ij
    multiplier = -0.4 / c->v3(k);
    double dv3f2bydv = 0.0;
    if (c->Nf2() > 0)
    {
        for (x2 = 1; x2 <= c->N3f2(); ++x2) // N3f2 is whether z is included or not
        {
            dv3f2bydv = f2->ddv_v3(vlagged, c, i, j, k, x1, &x2);
            E->inc(vtemp, dv3f2bydv * multiplier, i, j, &x2);
        }
        // 2/5vdiv f2
        multiplier = 0.4 * c->v(k);

        for (x2 = 1; x2 <= c->N3f2(); ++x2) // N3f2 is whether z is included or not
        {
            temp_sten = (*O->ddxi(i, j, &x2)) * multiplier;
            f2->insert(&temp_sten, vtemp, i, j, k, x1, &x2);
        }
    }

    // bool uw; int dir;

    if (IMPACTA_ions::ion_motion)
    {
        // Ion motion terms
        // + C.grad f1
        for (IMPACT_Dim x4 = 1; x4 <= c->Nf1(); ++x4)
        {

            temp_sten = (*O->ddxi(i, j, &x4)) *
                        Initial_Conditions::C_i[x4.get() - 1].get(i, j);

            // Upwind
            /*dir = x4.get()-1;
            uw = Initial_Conditions::C_i[dir].sign(i,j);
            temp_sten = (*O->ddxi_uw(i,j,&x4,uw))*
            Initial_Conditions::C_i[dir].get(i,j);*/
            f1->insert(&temp_sten, vtemp, i, j, k, x1);
        }

        /*  // Cxwdf0/dv - now explicit - agrt 2011
        GetOrthogonal(x1,&x2,&x3);
        if (x3>3-c->NB())
      B->inc(vtemp,-f0->ddv(vlagged,c,i,j,k)
             *Initial_Conditions::C_i[x2.get()-1].get(i,j),i,j,&x3);
        if (x2>3-c->NB())
      B->inc(vtemp,f0->ddv(vlagged,c,i,j,k)
          *Initial_Conditions::C_i[x3.get()-1].get(i,j),i,j,&x2);*/
    }
}
// code for outer cells of f0 equation - i.e.  boundaries
inline void f1equation_outer(IMPACT_Config *c, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *E, IMPACT_Var *B, IMPACT_Var *f1, IMPACT_Var *f2, int *i, int *j, int *k, IMPACT_Dim *x1)
{
    /*
  **********************************************************
  INSERT f1 Equation terms CODE HERE - BC versions
  **********************************************************
  */

    f1->set(vtemp, equation_switches::e_inert_on * c->idt(), i, j, k, x1);
    // set f1 n+1 element to 1.0

    // vgrad f0
    double v = c->v(k) * equation_switches::inf1_vgradf0_on;
    IMPACT_stencil temp_sten = (*O->ddxi(i, j, x1)) * v;
    f0->insert_BC(&temp_sten, vtemp, i, j, k);

    if (x1->get() == 3)
    {
        f0->inc(vtemp, v * IMPACT_Heating::Dnbydzgrid.get(i, j), i, j, k);
        f0->inc(vtemp, (v * v / Initial_Conditions::Te.get(i, j) - 1.5) * v * IMPACT_Heating::DTbydzgrid.get(i, j), i, j, k);
    }

    //-Edf0/dv
    if (!if_use_hybrid_nonlinear)
    {
        if (c->NE() > 0)
            E->set(vtemp, -equation_switches::inf1_Edf0dv_on * f0->ddv(vlagged, c, i, j, k), i, j, x1);
    }
    else
    {
        if (c->NE() > 0)
        {
            E->set(vtemp, -0.5 * equation_switches::inf1_Edf0dv_on * f0->ddv(vlagged, c, i, j, k), i, j, x1);
            IMPACT_Vel_Sten ddv = O->ddv(c, k);
            ddv *(-0.5 * equation_switches::inf1_Edf0dv_on *E->get(vlagged, i, j, x1));
            f0->insert_BC(&ddv, vtemp, i, j, k);
        }
    }
    // double erfv=erf(v);
    // double expmv2oversqrtpi = exp(-1.0*c->v2(k))/sqrt(globalconsts::pi);

    double ei_cols = equation_switches::Cei_on * Initial_Conditions::Z.get(i, j) * Initial_Conditions::Z.get(i, j) * Initial_Conditions::ni.get(i, j) / c->v3(k) / ((Initial_Conditions::Z.get(i, j) / globalconsts::oneover_atomic_Z + 0.236) / (1.0 + 0.236 * Initial_Conditions::Z.get(i, j) / globalconsts::oneover_atomic_Z) * .236);

    // double ei_cols=(equation_switches::Cei_on*Initial_Conditions::Z.get(i,j)*Initial_Conditions::Z.get(i,j)*Initial_Conditions::ni.get(i,j)
    //    -(erfv*(0.5/c->v2(k)-1.0)-expmv2oversqrtpi/c->v(k))
    //    )/c->v3(k)-2.0*f0->get(vlagged,i,j,k);

    // 0*(expmv2oversqrtpi-globalconsts::fourthirds*
    //    (expmv2oversqrtpi*c->v2(k)+c->v3(k)*(1.0-erfv))+erfv*(c->v(k)-0.5/c->v(k)))

    // collision term - only ion collisions first iteration
    f1->inc(vtemp, ei_cols, i, j, k, x1);

    // B x f1
    double multiplier = -equation_switches::inf1_f1xB_on;
    // This cyclicly works out orthogonal correct coordinates for x1 andx2

    IMPACT_Dim x2, x3;
    GetOrthogonal(x1, &x2, &x3);
    if (!if_use_hybrid_nonlinear)
    {
        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            f1->inc(vtemp, B->get(vlagged, i, j, &x2) * multiplier * Bimp, i, j, k, &x3);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            f1->inc(vtemp, -(B->get(vlagged, i, j, &x3) * multiplier * Bimp), i, j, k, &x2);

        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            B->inc(vtemp, f1->get(vlagged, i, j, k, &x3) * multiplier * (1 - Bimp), i, j, &x2);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            B->inc(vtemp, -f1->get(vlagged, i, j, k, &x2) * multiplier * (1 - Bimp), i, j, &x3);
    }
    else
    {
        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            f1->inc(vtemp, B->get(vlagged, i, j, &x2) * multiplier * 0.5, i, j, k, &x3);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            f1->inc(vtemp, -(B->get(vlagged, i, j, &x3) * multiplier * 0.5), i, j, k, &x2);

        if (x2 > 3 - c->NB() && x3 <= c->Nf1())
            B->inc(vtemp, f1->get(vlagged, i, j, k, &x3) * multiplier * 0.5, i, j, &x2);
        if (x3 > 3 - c->NB() && x2 <= c->Nf1())
            B->inc(vtemp, -f1->get(vlagged, i, j, k, &x2) * multiplier * 0.5, i, j, &x3);
    } 

    // f2 parts
    // a_j d/dvv^3f2_ij
    multiplier = -0.4 / c->v3(k);
    double dv3f2bydv = 0.0;
    if (c->Nf2() > 0)
    {
        for (x2 = 1; x2 <= c->N3f2(); ++x2) // N3f2 is whether z is included or not
        {
            dv3f2bydv = f2->ddv_v3_BC(vlagged, c, i, j, k, x1, &x2);
            E->inc(vtemp, dv3f2bydv * multiplier, i, j, &x2);
        }
        // 2/5vdiv f2
        multiplier = 0.4 * c->v(k);

        for (x2 = 1; x2 <= c->N3f2(); ++x2) // N3f2 is whether z is included or not
        {
            temp_sten = (*O->ddxi(i, j, &x2)) * multiplier;
            f2->insert_BC(&temp_sten, vtemp, i, j, k, x1, &x2);
        }
    }

    // bool uw; int dir;

    if (IMPACTA_ions::ion_motion)
    {
        // Ion motion terms
        // + C.grad f1
        for (IMPACT_Dim x4 = 1; x4 <= c->Nf1(); ++x4)
        {
            temp_sten = (*O->ddxi(i, j, &x4)) *
                        Initial_Conditions::C_i[x4.get() - 1].get(i, j);

            // Upwind
            /*dir = x4.get()-1;
                uw = Initial_Conditions::C_i[dir].sign(i,j);
                temp_sten = (*O->ddxi_uw(i,j,&x4,uw))*
                Initial_Conditions::C_i[dir].get(i,j);*/

            f1->insert_f1BC(&temp_sten, vtemp, i, j, k, x1);
        }

        /* // Cxwdf0/dv - now explicit agrt 2011
        GetOrthogonal(x1,&x2,&x3);
        if (x3>3-c->NB())
      B->inc(vtemp,-f0->ddv_BC(vlagged,c,i,j,k)
             *Initial_Conditions::C_i[x2.get()-1].get(i,j),i,j,&x3);
        if (x2>3-c->NB())
      B->inc(vtemp,f0->ddv_BC(vlagged,c,i,j,k)
          *Initial_Conditions::C_i[x3.get()-1].get(i,j),i,j,&x2); */
    }
}
