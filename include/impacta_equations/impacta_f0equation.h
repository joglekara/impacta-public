/*
**********************************************************
 f0 equation IMPACTA code for inner and outer cells

Version 1.2
1
AGRT
20/2/07
24/2/08 - updated to include simple ion model

**********************************************************

*/
// code for inner cells of f0 equation - i.e. not at boundaries
void f0equation_inner(IMPACT_Config *c, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *E, IMPACT_Var *B, IMPACT_Var *f1, int *i, int *j, int *k)
{

    /*
  **********************************************************
  INSERT f0 Equation terms CODE HERE
  **********************************************************
     */
    f0->set(vtemp, equation_switches::df0_by_dt_on * c->idt(), i, j, k);

    // set f0 n+1 element to 1.0
    IMPACT_Dim coord_1(1);
    IMPACT_Dim coord_2(2);
    // v/3 div.f1
    double vover3 = c->v(k) * oneover3 * equation_switches::inf0_grad_f1_on;
    // i.e. v/3

    IMPACT_stencil temp_sten = (*O->ddx(i)) * vover3;
    f1->insert(&temp_sten, vtemp, i, j, k, &coord_1); // x part of div
    temp_sten = (*O->ddy(j)) * vover3;
    f1->insert(&temp_sten, vtemp, i, j, k, &coord_2); // y part of div

    // ai/v^2 dv^2f1i/fv
    double oneover3v2 = -1.0 * oneover3 / c->v2(k); // 1/v squared
    if (!if_use_hybrid_nonlinear)
    {
        for (IMPACT_Dim x1 = 1; x1 <= c->NE(); ++x1)
            E->inc(vtemp, oneover3v2 * f1->ddv_v2(vlagged, c, i, j, k, &x1) * equation_switches::inf0_Edf1dv_on, i, j, &x1); // E element
    }
    else
    {
        IMPACT_Vel_Sten ddv_v2;
        for (IMPACT_Dim x1 = 1; x1 <= c->NE(); ++x1)
        {
            E->inc(vtemp, 0.5 * oneover3v2 * f1->ddv_v2(vlagged, c, i, j, k, &x1) * equation_switches::inf0_Edf1dv_on, i, j, &x1); // E element
            ddv_v2 = O->ddv_v2(c,k);
            ddv_v2 *(0.5 * oneover3v2 * E->get(vlagged, i, j, &x1) * equation_switches::inf0_Edf1dv_on);
            f1->insert(&ddv_v2, vtemp, i, j, k, &x1);
        }
    }
    // e-e collision term Cee0
    int kminus = *k - 1;
    double Cee0_const = -1.0 / c->v2(k) * c->idv(k);
    IMPACT_Vel_Sten Fk(0.0, 0.0, 0.0);
    Fk.inc(0, -Cee0_flux_store::flux->get(i, j, &kminus)->get(1));
    Fk.inc(1, -Cee0_flux_store::flux->get(i, j, &kminus)->get(2) + Cee0_flux_store::flux->get(i, j, k)->get(1));
    Fk.inc(2, Cee0_flux_store::flux->get(i, j, k)->get(2));
    Fk *Cee0_const;
    f0->insert(&Fk, vtemp, i, j, k);

    if (IMPACTA_ions::ion_motion)
    {
        // ION Motion terms...
        // + C.grad f0
        for (IMPACT_Dim x1 = 1; x1 <= c->Nf1(); ++x1)
        {
            temp_sten = (*O->ddxi(i, j, &x1)) *
                        Initial_Conditions::C_i[x1.get() - 1].get(i, j);

            f0->insert(&temp_sten, vtemp, i, j, k);
        }

        // Removed as not included in f2 and f1 agrt 2011  -- readded May 2011
        // - v/3df0/dvdivCi
        kminus = *k - 1;
        int kplus = *k + 1;
        double vover3deltav = -c->v(k) * oneover3 / (c->dv(&kminus) + c->dv(&kplus));
        f0->inc(vtemp, vover3deltav * Initial_Conditions::DivC_i.get(i, j), i, j, &kplus);
        f0->inc(vtemp, -vover3deltav * Initial_Conditions::DivC_i.get(i, j), i, j, &kminus);
    }
}

// code for outer cells of f0 equation - i.e.  boundaries
void f0equation_outer(IMPACT_Config *c, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *E, IMPACT_Var *B, IMPACT_Var *f1, int *i, int *j, int *k)
{
    /*
  **********************************************************
  INSERT f0 Equation terms CODE HERE - BC versions
  **********************************************************
  */
    f0->set(vtemp, equation_switches::df0_by_dt_on * c->idt(), i, j, k);

    IMPACT_Dim coord_1(1);
    IMPACT_Dim coord_2(2);

    // v/3 div.f1
    double vover3 = c->v(k) * oneover3 * equation_switches::inf0_grad_f1_on;
    // i.e. v/3

    IMPACT_stencil temp_sten = (*O->ddx(i)) * vover3;
    f1->insert_f1BC(&temp_sten, vtemp, i, j, k, &coord_1); // x part of div
    temp_sten = (*O->ddy(j)) * vover3;
    f1->insert_f1BC(&temp_sten, vtemp, i, j, k, &coord_2); // y part of div

    // ai/v^2 dv^2f1i/fv
    double oneover3v2 = -1.0 * oneover3 / c->v2(k); // 1/v squared

    if (!if_use_hybrid_nonlinear)
    {
        for (IMPACT_Dim x1 = 1; x1 <= c->NE(); ++x1)
            E->inc(vtemp, oneover3v2 * f1->ddv_v2_BC(vlagged, c, i, j, k, &x1) * equation_switches::inf0_Edf1dv_on, i, j, &x1);
    }
    else
    {
        IMPACT_Vel_Sten ddv_v2;
        for (IMPACT_Dim x1 = 1; x1 <= c->NE(); ++x1)
        {
            E->inc(vtemp, 0.5 * oneover3v2 * f1->ddv_v2_BC(vlagged, c, i, j, k, &x1) * equation_switches::inf0_Edf1dv_on, i, j, &x1); // E element
            ddv_v2 = O->ddv_v2(c,k);
            ddv_v2 * (0.5 * oneover3v2 * E->get(vlagged, i, j, &x1) * equation_switches::inf0_Edf1dv_on);
            f1->insert_BC(&ddv_v2, vtemp, i, j, k, &x1);
        }
    }

    // e-e collision term Cee0
    int kminus = *k - 1;
    double Cee0_const = -1.0 / c->v2(k) * c->idv(k);
    IMPACT_Vel_Sten Fk(0.0, 0.0, 0.0);
    Fk.inc(0, -Cee0_flux_store::flux->get(i, j, &kminus)->get(1));
    Fk.inc(1, -Cee0_flux_store::flux->get(i, j, &kminus)->get(2) + Cee0_flux_store::flux->get(i, j, k)->get(1));
    Fk.inc(2, Cee0_flux_store::flux->get(i, j, k)->get(2));
    Fk *Cee0_const;
    f0->insert(&Fk, vtemp, i, j, k);

    // bool uw; int dir;

    if (IMPACTA_ions::ion_motion)
    {
        // ION Motion terms...
        // C.grad f0
        for (IMPACT_Dim x1 = 1; x1 <= c->Nf1(); ++x1)
        {

            temp_sten = (*O->ddxi(i, j, &x1)) *
                        Initial_Conditions::C_i[x1.get() - 1].get(i, j);

            f0->insert_BC(&temp_sten, vtemp, i, j, k);
        }

        // Removed as not included in f2 and f1 agrt 2011  -- readded May 2011
        // v/3df0/dvdivCi
        kminus = *k - 1;
        int kplus = *k + 1;
        if (kminus < 1)
            kminus = 1;
        if (kminus > c->Nv())
            kplus = c->Nv();
        double vover3deltav = -c->v(k) * oneover3 / (c->dv(&kminus) + c->dv(&kplus));
        f0->inc(vtemp, -vover3deltav * Initial_Conditions::DivC_i.get(i, j), i, j, &kminus);
        f0->inc(vtemp, vover3deltav * Initial_Conditions::DivC_i.get(i, j), i, j, &kplus);
    }
}
