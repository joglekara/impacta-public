/*
**********************************************************
Ion motion package for IMPACT

Version 1.1
AGRT

2/7/08 - Fixed boundaries for ions added
22/2/08
jan 2011 - Added DCbyDt function - RHS of
DC/Dt in hydrodynamic equations, also used
in RHS of main kinetic equations (refer to
Thomas JCP 2011)

march 2011 - Added dZbydt - Z has to move too!

july 2011 - Added ion pressure to equation of motion!
**********************************************************
*/
void IMPACTA_Set_Ci(IMPACT_Config *c, double value)
{
    if (IMPACTA_ions::ion_motion)
    {
        for (int i = (1 + IMPACT_Boundaries::fixed_any[0]); i <= (c->Nx() - IMPACT_Boundaries::fixed_any[1]); ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {
                for (int dir = 0; dir < 3; ++dir)
                {
                    Initial_Conditions::C_i[dir].set(i, j, value);
                }
                Initial_Conditions::DivC_i.set(i, j, value);
            }
    }
}
void IMPACTA_Multiply_Ci(IMPACT_Config *c, double value)
{
    if (IMPACTA_ions::ion_motion)
    {
        for (int i = (1 + IMPACT_Boundaries::fixed_any[0]); i <= (c->Nx() - IMPACT_Boundaries::fixed_any[1]); ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
                for (int dir = 0; dir < 3; ++dir)
                    Initial_Conditions::C_i[dir].set(i, j, value * Initial_Conditions::C_i[dir].get(&i, &j));
    }
}
void IMPACTA_Update_DivCi(IMPACT_Config *c, IMPACT_MPI_Config *M,
                          IMPACT_StenOps *O)
{
    MPI_Barrier(MPI_COMM_WORLD);
    // Now update DivC
    int direc = 0, iplus = 0, iminus = 0, jplus = 0, jminus = 0;
    IMPACT_stencil temp_sten;
    int istart = M->istart();
    int iend = M->iend();
    // New -> correct istart and iend for fixed boundaries:
    if (istart == 1)
    {
        istart = (1 + IMPACT_Boundaries::fixed_any[0]);
    }
    if (iend == c->Nx())
    {
        iend = (c->Nx() - IMPACT_Boundaries::fixed_any[1]);
    }

    //bool //uw = 0;

    if (IMPACTA_ions::ion_motion)
        for (int i = istart; i <= iend; ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {
                double divC = 0.0;
                for (IMPACT_Dim x1 = 1; x1 <= 3; ++x1)
                {
                    iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;
                    direc = x1.get() - 1;
                    //uw = Initial_Conditions::C_i[direc].sign(&i, &j);
                    //	    temp_sten = (*O->ddxi_uw(&i,&j,&x1,uw));

                    temp_sten = (*O->ddxi(&i, &j, &x1));
                    IMPACT_Ci_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                                    &iminus, &jplus, &jminus, &x1);

                    divC += Initial_Conditions::C_i[direc].get(&i, &j) * temp_sten(0);
                    divC += Initial_Conditions::C_i[direc].get(&iplus, &j) * temp_sten(1);
                    divC += Initial_Conditions::C_i[direc].get(&iminus, &j) * temp_sten(2);
                    divC += Initial_Conditions::C_i[direc].get(&i, &jplus) * temp_sten(3);
                    divC += Initial_Conditions::C_i[direc].get(&i, &jminus) * temp_sten(4);
                }
                Initial_Conditions::DivC_i.set(i, j, divC);
            }
    MPI_Barrier(MPI_COMM_WORLD);
    IMPACTA_Share_Moment(c, M, &Initial_Conditions::DivC_i);
}

// this is RHS of DU/Dt equation - refer to Thomas JCP 2011
double IMPACTA_DUbyDt(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_ParVec *v,
                      IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1,
                      IMPACT_Var *f2, IMPACT_Var *E, IMPACT_Var *B,
                      int i, int j, IMPACT_Dim x1, int n_adjusted)
{
    IMPACT_stencil temp_sten;
    double gradp = 0.0, jxb = 0.0, EZnimne = 0.0;
    int dir = 0;
    double loc_ni = Initial_Conditions::ni.get(&i, &j);
    double loc_Z = Initial_Conditions::Z.get(&i, &j);
    int iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;
    int iplusni = iplus, iminusni = iminus; // these added because periodic conditions different for f and e.g. n,T,Z grids!!!
    IMPACT_Dim x2, x3;
    //bool uw;
    //-------------------------------------------------------
    // +Grad p_e

    //uw = Initial_Conditions::C_i[x1.get() - 1].sign(&i, &j);
    // temp_sten = (*O->ddxi_uw(&i,&j,&x1,uw));
    temp_sten = (*O->ddxi(&i, &j, &x1));

    IMPACT_f0_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                    &iminus, &jplus, &jminus);
    IMPACT_ni_bound(&temp_sten, c->Nx(), c->Ny(), &iplusni,
                    &iminusni, &jplus, &jminus);
    // here - added ion pressure also!
    /*gradp=(Local_pe(v,f0,c,&i,&j)+IMPACTA_ions::ion_temperature*Initial_Conditions::ni.get(&i,&j))*temp_sten(0);
    gradp+=(Local_pe(v,f0,c,&iplus,&j)+IMPACTA_ions::ion_temperature*Initial_Conditions::ni.get(&iplusni,&j))*temp_sten(1);
    gradp+=(Local_pe(v,f0,c,&iminus,&j)+IMPACTA_ions::ion_temperature*Initial_Conditions::ni.get(&iminusni,&j))*temp_sten(2);
    gradp+=(Local_pe(v,f0,c,&i,&jplus)+IMPACTA_ions::ion_temperature*Initial_Conditions::ni.get(&i,&jplus))*temp_sten(3);
    gradp+=(Local_pe(v,f0,c,&i,&jminus)+IMPACTA_ions::ion_temperature*Initial_Conditions::ni.get(&i,&jminus))*temp_sten(4);
    */
    gradp = (Local_pe(v, f0, c, &i, &j)) * temp_sten(0);
    gradp += (Local_pe(v, f0, c, &iplus, &j)) * temp_sten(1);
    gradp += (Local_pe(v, f0, c, &iminus, &j)) * temp_sten(2);
    gradp += (Local_pe(v, f0, c, &i, &jplus)) * temp_sten(3);
    gradp += (Local_pe(v, f0, c, &i, &jminus)) * temp_sten(4);

    //-------------------------------------------------------
    // Divergence of anisotropic pressure
    if (c->Nf2() > 0)
        for (IMPACT_Dim xa = 1; xa <= 3; ++xa)
        {
            iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;
            //   //uw = Initial_Conditions::C_i[xa.get()-1].sign(&i,&j);
            //	temp_sten = (*O->ddxi_uw(&i,&j,&xa,uw));

            temp_sten = (*O->ddxi(&i, &j, &xa));
            IMPACT_f2_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                            &iminus, &jplus, &jminus, &x1, &xa);
            gradp += Local_Pi(v, f2, c, &i, &j, &x1, &xa) * temp_sten(0);
            gradp += Local_Pi(v, f2, c, &iplus, &j, &x1, &xa) * temp_sten(1);
            gradp += Local_Pi(v, f2, c, &iminus, &j, &x1, &xa) * temp_sten(2);
            gradp += Local_Pi(v, f2, c, &i, &jplus, &x1, &xa) * temp_sten(3);
            gradp += Local_Pi(v, f2, c, &i, &jminus, &x1, &xa) * temp_sten(4);
        }
    //-------------------------------------------------------

    if (dir == 2)
    {
        gradp = (Local_Te(v, f0, c, &i, &j) * IMPACT_Heating::Dnz_xy.get(&i, &j) * IMPACT_Heating::Dnz_t.Get(n_adjusted) + Local_ne(v, f0, c, &i, &j) * IMPACT_Heating::DTz_xy.get(&i, &j) * IMPACT_Heating::DTz_t.Get(n_adjusted));
    }

    //-------------------------------------------------------
    //// +J x B
    /// NEW 2011 AGRT -> +jxB instead of curl BxB
    GetOrthogonal(&x1, &x2, &x3);
    if (x3 > 3 - c->NB())
    {
        jxb += B->get(v, &i, &j, &x3) * Local_je(v, f1, c, &i, &j, &x2);
    }
    if (x2 > 3 - c->NB())
    {
        jxb -= B->get(v, &i, &j, &x2) * Local_je(v, f1, c, &i, &j, &x3);
    }

    //-------------------------------------------------------
    // +E(Zni-ne)
    if (x1 <= c->NE())
    {
        EZnimne = E->get(v, &i, &j, &x1) * (loc_Z * loc_ni - Local_ne(v, f0, c, &i, &j));
    }
    //-------------------------------------------------------

    return (IMPACTA_ions::alpha_ion / loc_ni) * (EZnimne + jxb - gradp);

    /*
      old (pre 2011)
      (((gradp+EZnimne)
           +gradBsquared*IMPACTA_ions::a2bar)
           *Initial_Conditions::Z.get(&i,&j));*/
}
void IMPACTA_Ci_Smooth(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_StenOps *O, int direc);
int IMPACTA_Update_Ci(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_ParVec *v,
                      IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1,
                      IMPACT_Var *f2, IMPACT_Var *E, IMPACT_Var *B, int *n)
{
    int istart = M->istart();
    int iend = M->iend();
    // New -> correct istart and iend for fixed boundaries:
    if (istart == 1)
    {
        istart = (1 + IMPACT_Boundaries::fixed_any[0]);
    }
    if (iend == c->Nx())
    {
        iend = (c->Nx() - IMPACT_Boundaries::fixed_any[1]);
    }

    IMPACT_stencil temp_sten;
    int dir = 0, fixcheck = 0;
    int n_adjusted = *n;
    if (n_adjusted < 1)
        n_adjusted = 1;
    double DUbyDt = 0.0;
    int iplus = 2, iminus = 0, jplus = 2, jminus = 0;
    IMPACT_Dim x2, x3;

    // First update C due to effects of thermal and magnetic pressure:
    //________________________________________________________________________
    // Contains smoothing effect!!!
    if (!IMPACTA_ions::ion_motion)
        return 0;

    for (dir = 0; dir < 3; ++dir)
    {
        Initial_Conditions::C_istar[dir].copy(&Initial_Conditions::C_i[dir]);
    }

    for (IMPACT_Dim x1 = 1; x1 <= 3; ++x1)
    {

        dir = x1.get() - 1;
        for (int i = istart; i <= iend; ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {

                fixcheck = 0;
                // check fixed bounds
                if (i == 1)
                    fixcheck += (IMPACT_Boundaries::fix_Ci[0] != 0.0);
                if (i == c->Nx())
                    fixcheck += (IMPACT_Boundaries::fix_Ci[1] != 0.0);
                if (j == 1)
                    fixcheck += (IMPACT_Boundaries::fix_Ci[2] != 0.0);
                if (j == c->Ny())
                    fixcheck += (IMPACT_Boundaries::fix_Ci[3] != 0.0);

                DUbyDt = IMPACTA_DUbyDt(c, M, v, O, f0, f1, f2, E, B, i, j, x1, n_adjusted);

                if (!fixcheck)
                    Initial_Conditions::C_istar[dir].set(i, j, Initial_Conditions::C_i[dir].get(&i, &j) + DUbyDt * c->dt());
                else
                    Initial_Conditions::C_istar[dir].set(i, j, Initial_Conditions::C_i[dir].get(&i, &j));
            }
    }
    for (dir = 0; dir < 3; ++dir)
    {
        Initial_Conditions::C_i[dir].copy(&Initial_Conditions::C_istar[dir]);
        MPI_Barrier(MPI_COMM_WORLD);
        IMPACTA_Share_Moment(c, M, &Initial_Conditions::C_i[dir]);
    }
    //________________________________________________________________________
    // Next - advection at the bulk flow velocity
    double deltaCistar = 0.0;
    int dirj = 0;
    //bool uw;
    for (IMPACT_Dim x1i = 1; x1i <= 3; ++x1i)
    {
        dir = x1i.get() - 1;
        for (int i = istart; i <= iend; ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {
                fixcheck = 0;
                // check fixed bounds
                if (i == 1)
                    fixcheck += (IMPACT_Boundaries::fix_Ci[0] != 0.0);
                if (i == c->Nx())
                    fixcheck += (IMPACT_Boundaries::fix_Ci[1] != 0.0);
                if (j == 1)
                    fixcheck += (IMPACT_Boundaries::fix_Ci[2] != 0.0);
                if (j == c->Ny())
                    fixcheck += (IMPACT_Boundaries::fix_Ci[3] != 0.0);
                deltaCistar = 0.0;
                for (IMPACT_Dim x1j = 1; x1j <= 3; ++x1j)
                {
                    iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;
                    dirj = x1j.get() - 1;

                    //uw = Initial_Conditions::C_i[dir].sign(&i, &j);
                    // temp_sten = (*O->ddxi_uw(&i,&j,&x1j,uw));

                    temp_sten = (*O->ddxi(&i, &j, &x1j));

                    IMPACT_Ci_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                                    &iminus, &jplus, &jminus, &x1j);

                    deltaCistar += Initial_Conditions::C_i[dir].get(&i, &j) * Initial_Conditions::C_i[dirj].get(&i, &j) * temp_sten(0);
                    deltaCistar += Initial_Conditions::C_i[dir].get(&i, &j) * Initial_Conditions::C_i[dirj].get(&iplus, &j) * temp_sten(1);
                    deltaCistar += Initial_Conditions::C_i[dir].get(&i, &j) * Initial_Conditions::C_i[dirj].get(&iminus, &j) * temp_sten(2);
                    deltaCistar += Initial_Conditions::C_i[dir].get(&i, &j) * Initial_Conditions::C_i[dirj].get(&i, &jplus) * temp_sten(3);
                    deltaCistar += Initial_Conditions::C_i[dir].get(&i, &j) * Initial_Conditions::C_i[dirj].get(&i, &jminus) * temp_sten(4);
                } // end of xj loop
                if (!fixcheck)
                    Initial_Conditions::C_istar[dir].set(i, j, Initial_Conditions::C_i[dir].get(&i, &j) - deltaCistar * c->dt());
                else
                    Initial_Conditions::C_istar[dir].set(i, j, Initial_Conditions::C_i[dir].get(&i, &j));
            } // end of j loop
    }         // end of xi loop
    for (dir = 0; dir < 3; ++dir)
    {
        Initial_Conditions::C_i[dir].copy(&Initial_Conditions::C_istar[dir]);
        MPI_Barrier(MPI_COMM_WORLD);
        IMPACTA_Share_Moment(c, M, &Initial_Conditions::C_i[dir]);
    }
    /*// smooth Ci
    for (int dir=0;dir<3;++dir)
    {
        IMPACTA_Ci_Smooth(c,M,O,dir);
    }*/
    // Finally update div Ci
    IMPACTA_Update_DivCi(c, M, O);
    return 0;
} // end of function

int IMPACTA_Update_ni(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_ParVec *v,
                      IMPACT_StenOps *O, IMPACT_Var *f0, int *n)
{
    int istart = M->istart();
    int iend = M->iend();
    // New -> correct istart and iend for fixed boundaries:
    if (istart == 1)
    {
        istart = (1 + IMPACT_Boundaries::fixed_any[0]);
    }
    if (iend == c->Nx())
    {
        iend = (c->Nx() - IMPACT_Boundaries::fixed_any[1]);
    }

    int iplus = 2, iminus = 0, jplus = 2, jminus = 0, direc = 0;
    double gradni = 0.0, divC = 0.0, Citemp = 0.0;
    IMPACT_stencil temp_sten;
    int fixcheck = 0;
    int n_adjusted = *n;
    if (n_adjusted < 1)
        n_adjusted = 1;
    if (!IMPACTA_ions::ion_motion)
        return 0;
    Initial_Conditions::C_istar[0].copy(&Initial_Conditions::ni);

    if (IMPACTA_ions::quasin_on)
    {
        for (int i = istart; i <= iend; ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {
                //-------------------------------------------------------
                // C.Grad ni
                gradni = 0.0;
                fixcheck = 0;
                // check fixed bounds
                if (i == 1)
                    fixcheck += IMPACT_Boundaries::fix_ni[0];
                if (i == c->Nx())
                    fixcheck += IMPACT_Boundaries::fix_ni[1];
                if (j == 1)
                    fixcheck += IMPACT_Boundaries::fix_ni[2];
                if (j == c->Ny())
                    fixcheck += IMPACT_Boundaries::fix_ni[3];
                Initial_Conditions::ni.set(i, j, Local_ne(v, f0, c, &i, &j) / Initial_Conditions::Z.get(&i, &j));
            }
    }
    else
    {
        for (int i = istart; i <= iend; ++i)
        {
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {
                //-------------------------------------------------------
                // C.Grad ni
                gradni = 0.0;
                fixcheck = 0;
                // check fixed bounds
                if (i == 1)
                    fixcheck += IMPACT_Boundaries::fix_ni[0];
                if (i == c->Nx())
                    fixcheck += IMPACT_Boundaries::fix_ni[1];
                if (j == 1)
                    fixcheck += IMPACT_Boundaries::fix_ni[2];
                if (j == c->Ny())
                    fixcheck += IMPACT_Boundaries::fix_ni[3];

                for (IMPACT_Dim x1 = 1; x1 <= 3; ++x1)
                {

                    direc = x1.get() - 1;
                    iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;

                    temp_sten = (*O->ddxi(&i, &j, &x1));

                    IMPACT_ni_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                                    &iminus, &jplus, &jminus);

                    Citemp = Initial_Conditions::C_i[direc].get(&i, &j);

                    gradni += Citemp * Initial_Conditions::ni.get(&i, &j) * temp_sten(0);
                    gradni += Citemp * Initial_Conditions::ni.get(&iplus, &j) * temp_sten(1);
                    gradni += Citemp * Initial_Conditions::ni.get(&iminus, &j) * temp_sten(2);
                    gradni += Citemp * Initial_Conditions::ni.get(&i, &jplus) * temp_sten(3);
                    gradni += Citemp * Initial_Conditions::ni.get(&i, &jminus) * temp_sten(4);

                    if (direc == 2)
                    {
                        gradni = Citemp * IMPACT_Heating::Dnz_xy.get(&i, &j) * IMPACT_Heating::Dnz_t.Get(n_adjusted);
                    }
                }
                //-------------------------------------------------------
                // ni div C
                divC = Initial_Conditions::DivC_i.get(&i, &j) * Initial_Conditions::ni.get(&i, &j);
                //-------------------------------------------------------

                // using C_istar as a temporary store
                if (!fixcheck)
                    Initial_Conditions::C_istar[0].set(i, j, Initial_Conditions::ni.get(&i, &j) - (gradni + divC) * c->dt());
                else
                    Initial_Conditions::C_istar[0].set(i, j, Initial_Conditions::ni.get(&i, &j));
            }
        }
        Initial_Conditions::ni.copy(&Initial_Conditions::C_istar[0]);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    IMPACTA_Share_Moment(c, M, &Initial_Conditions::ni);

    return 0;
}
// This advects the charge state (but assumes Z is incompressible quantity - dZ/dt=0
int IMPACTA_Update_Z(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_ParVec *v,
                     IMPACT_StenOps *O, IMPACT_Var *f0, int *n)
{
    int istart = M->istart();
    int iend = M->iend();
    // New -> correct istart and iend for fixed boundaries:
    if (istart == 1)
    {
        istart = (1 + IMPACT_Boundaries::fixed_any[0]);
    }
    if (iend == c->Nx())
    {
        iend = (c->Nx() - IMPACT_Boundaries::fixed_any[1]);
    }

    int iplus = 2, iminus = 0, jplus = 2, jminus = 0, direc = 0;
    double gradni = 0.0, Citemp = 0.0;
    IMPACT_stencil temp_sten;
    int fixcheck = 0;
    int n_adjusted = *n;
    if (n_adjusted < 1)
        n_adjusted = 1;
    if (!IMPACTA_ions::ion_motion)
        return 0;
    Initial_Conditions::C_istar[0].copy(&Initial_Conditions::Z);
    for (int i = istart; i <= iend; ++i)
        for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
        {
            //-------------------------------------------------------
            // C.Grad ni
            gradni = 0.0;
            fixcheck = 0;
            // check fixed bounds
            if (i == 1)
                fixcheck += IMPACT_Boundaries::fix_ni[0];
            if (i == c->Nx())
                fixcheck += IMPACT_Boundaries::fix_ni[1];
            if (j == 1)
                fixcheck += IMPACT_Boundaries::fix_ni[2];
            if (j == c->Ny())
                fixcheck += IMPACT_Boundaries::fix_ni[3];

            //bool uw;
            for (IMPACT_Dim x1 = 1; x1 <= 3; ++x1)
            {

                direc = x1.get() - 1;
                iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;

                //uw = Initial_Conditions::C_i[direc].sign(&i, &j);
                //	    temp_sten = (*O->ddxi_uw(&i,&j,&x1,uw));
                temp_sten = (*O->ddxi(&i, &j, &x1));

                IMPACT_ni_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                                &iminus, &jplus, &jminus);

                Citemp = Initial_Conditions::C_i[direc].get(&i, &j);

                gradni += Citemp * Initial_Conditions::Z.get(&i, &j) * temp_sten(0);
                gradni += Citemp * Initial_Conditions::Z.get(&iplus, &j) * temp_sten(1);
                gradni += Citemp * Initial_Conditions::Z.get(&iminus, &j) * temp_sten(2);
                gradni += Citemp * Initial_Conditions::Z.get(&i, &jplus) * temp_sten(3);
                gradni += Citemp * Initial_Conditions::Z.get(&i, &jminus) * temp_sten(4);

                if (direc == 2)
                {
                    gradni = Citemp * IMPACT_Heating::Dnz_xy.get(&i, &j) * IMPACT_Heating::Dnz_t.Get(n_adjusted);
                }
            }

            //-------------------------------------------------------

            // using C_istar as a temporary store
            if (!fixcheck)
                Initial_Conditions::C_istar[0].set(i, j, Initial_Conditions::Z.get(&i, &j) - (gradni)*c->dt());
            else
                Initial_Conditions::C_istar[0].set(i, j, Initial_Conditions::Z.get(&i, &j));
        }

    Initial_Conditions::Z.copy(&Initial_Conditions::C_istar[0]);
    MPI_Barrier(MPI_COMM_WORLD);
    IMPACTA_Share_Moment(c, M, &Initial_Conditions::Z);
    return 0;
}
void IMPACTA_Ci_Smooth(IMPACT_Config *c, IMPACT_MPI_Config *M,
                       IMPACT_StenOps *O, int direc)
{
    MPI_Barrier(MPI_COMM_WORLD);
    // Now update DivC
    int iplus = 0, iminus = 0, jplus = 0, jminus = 0;
    IMPACT_stencil temp_sten;
    int istart = M->istart();
    int iend = M->iend();
    // New -> correct istart and iend for fixed boundaries:
    if (istart == 1)
    {
        istart = (1 + IMPACT_Boundaries::fixed_any[0]);
    }
    if (iend == c->Nx())
    {
        iend = (c->Nx() - IMPACT_Boundaries::fixed_any[1]);
    }
    IMPACT_Dim x1(direc);
    if (IMPACTA_ions::ion_motion)
    {
        for (int i = istart; i <= iend; ++i)
            for (int j = (1 + IMPACT_Boundaries::fixed_any[2]); j <= (c->Ny() - IMPACT_Boundaries::fixed_any[3]); ++j)
            {
                double divC = 0.0;

                iplus = i + 1, iminus = i - 1, jplus = j + 1, jminus = j - 1;
                temp_sten = LAX;
                IMPACT_Ci_bound(&temp_sten, c->Nx(), c->Ny(), &iplus,
                                &iminus, &jplus, &jminus, &x1);

                divC += Initial_Conditions::C_i[direc].get(&i, &j) * temp_sten(0);
                divC += Initial_Conditions::C_i[direc].get(&iplus, &j) * temp_sten(1);
                divC += Initial_Conditions::C_i[direc].get(&iminus, &j) * temp_sten(2);
                divC += Initial_Conditions::C_i[direc].get(&i, &jplus) * temp_sten(3);
                divC += Initial_Conditions::C_i[direc].get(&i, &jminus) * temp_sten(4);
                Initial_Conditions::C_istar[direc].set(i, j, divC);
            }
        MPI_Barrier(MPI_COMM_WORLD);
        IMPACTA_Share_Moment(c, M, &Initial_Conditions::C_istar[direc]);
        Initial_Conditions::C_i[direc].copy(&Initial_Conditions::C_istar[direc]);
    }
}
