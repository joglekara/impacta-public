/*
**********************************************************
Obtain moments of the equations

Version 2.4
AGRT

23/3/07 - the first set of functions turn a ParVec into
a Moment, the second set do the reverse
4/6/07 - pressure tensor (f2) term added

21/8/07 - Average value member added
5/5/08 - Anisotrpic Pressure now 2* origninal -> because
of normalization
30/6/08 - Local anisotropic pressure function added
**********************************************************
*/
inline double Local_pe(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c,
                       int *i, int *j);
inline double Local_je(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c,
                       int *i, int *j, IMPACT_Dim *x1);
inline double Local_ne(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c,
                       int *i, int *j);

// First get moments from parallel vectors
// This is a moment...
struct IMPACT_Moment
{
    char name[LenName]; // the name of the moment
    int Nx, Ny;
    char coords[Ncoords]; // what is the coordinate system
    IMPACT_Vector xaxis;  // xaxis labels
    IMPACT_Vector yaxis;  // yaxis labels
    IMPACT_Matrix values; // the Matrix of moment values
    double average;
};
/*
Output format is like:
IMPACT
<Moment name>
<coords>  - either xy, rz or r0
<Nx> <TAB> <Ny>
<xaxis values>
<yaxis values>
<matrix of data>
*/
struct IMPACT_Moment Duplicate(struct IMPACT_Moment *M1)
{
    struct IMPACT_Moment M2;
    for (int i = 0; i < 20; ++i)
        M2.name[i] = M1->name[i];
    M2.Nx = M1->Nx;
    M2.Ny = M1->Ny;
    M2.coords[0] = M1->coords[0];
    M2.coords[1] = M1->coords[1];

    M2.xaxis.ChangeSize(M1->Nx);
    M2.yaxis.ChangeSize(M1->Ny);
    M2.values.ChangeSize(M1->Nx, M1->Ny);
    M2.average = M1->average;
    for (int i = 1; i <= M1->Nx; ++i)
        M2.xaxis.Set(i, M1->xaxis.Get(i));
    for (int i = 1; i <= M1->Ny; ++i)
        M2.yaxis.Set(i, M1->yaxis.Get(i));
    for (int i = 1; i <= M1->Nx; ++i)
        for (int j = 1; j <= M1->Ny; ++j)
            M2.values.set(i, j, M1->values.get(i, j));
    return M2;
}
struct IMPACT_Moment Divide(struct IMPACT_Moment *M1, struct IMPACT_Moment *M3)
{
    struct IMPACT_Moment M2;
    M2.Nx = M1->Nx;
    M2.Ny = M1->Ny;
    M2.coords[0] = M1->coords[0];
    M2.coords[1] = M1->coords[1];

    M2.xaxis.ChangeSize(M1->Nx);
    M2.yaxis.ChangeSize(M1->Ny);
    M2.values.ChangeSize(M1->Nx, M1->Ny);

    for (int i = 1; i <= M1->Nx; ++i)
        M2.xaxis.Set(i, M1->xaxis.Get(i));
    for (int i = 1; i <= M1->Ny; ++i)
        M2.yaxis.Set(i, M1->yaxis.Get(i));
    for (int i = 1; i <= M1->Nx; ++i)
        for (int j = 1; j <= M1->Ny; ++j)
            M2.values.set(i, j, M1->values.get(i, j) / M3->values.get(i, j));
    return M2;
}
inline void GetGrid(IMPACT_Config *c, struct IMPACT_Moment *M)
{
    M->xaxis.ChangeSize(c->Nx());
    M->yaxis.ChangeSize(c->Ny());
    for (int i = 1; i <= c->Nx(); ++i)
        M->xaxis.Set(i, c->xpos(i));
    for (int i = 1; i <= c->Ny(); ++i)
        M->yaxis.Set(i, c->ypos(i));
}
inline void new_Constant(IMPACT_IMatrix *Const, IMPACT_Config *c,
                         IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer,
                         std::string name1)
{

    int Nx = c->Nx();
    int Ny = c->Ny();

    strncpy(answer->name, name1.c_str(), LenName);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);

    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    double totalsum = 0.0;

    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            answer->values.set(i, j, Const->get(&i, &j));
            totalsum += Const->get(&i, &j);
        }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
}
inline void new_ne(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    strcpy(answer->name, "n_e");
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];

    double rowsum;
    double totalsum = 0.0;

    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            //  for (int k=1;k<=Nv;++k)
            //	    rowsum+=fourpi*f0->get(&gatheredvector,&i,&j,&k)*c->v2(&k)*c->dv(&k);
            for (int k = 1; k <= Nv; ++k)
                rowsum += fourpi * data[k - 1] * c->dv(&k) * c->v2(&k);
            answer->values.set(i, j, rowsum);
            totalsum += rowsum;
        }
    }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
// the next is for quickeness - no update of grid, just values.
inline void get_ne(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            //  for (int k=1;k<=Nv;++k)
            //	    rowsum+=fourpi*f0->get(&gatheredvector,&i,&j,&k)*c->v2(&k)*c->dv(&k);
            for (int k = 1; k <= Nv; ++k)
                rowsum += fourpi * data[k - 1] * c->dv(&k) * c->v2(&k);
            answer->values.set(i, j, rowsum);
        }
    }
    delete[] data;
}
inline void new_Ue(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    strcpy(answer->name, "U_e");
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum, totalsum = 0.0;
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += twopi * e_mass * data[k - 1] * c->v2(&k) * c->v2(&k) * c->dv(&k);
            // for (int k=1;k<=Nv;++k)
            // rowsum+= twopi*e_mass*f0->get(&gatheredvector,&i,&j,&k)*c->v2(&k)*c->v2(&k)*c->dv(&k); // 4pi /2
            answer->values.set(i, j, rowsum);
            totalsum += rowsum;
        }
    }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
inline void new_m0five(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    strcpy(answer->name, "m05");
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum, totalsum = 0.0;
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += twopi * e_mass * data[k - 1] * c->v3(&k) * c->v2(&k) * c->v2(&k) * c->dv(&k);
            // for (int k=1;k<=Nv;++k)
            // rowsum+= twopi*e_mass*f0->get(&gatheredvector,&i,&j,&k)*c->v2(&k)*c->v2(&k)*c->dv(&k); // 4pi /2
            answer->values.set(i, j, rowsum);
            totalsum += rowsum;
        }
    }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
inline void new_m0three(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    strcpy(answer->name, "m03");
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum, totalsum = 0.0;
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += twopi * e_mass * data[k - 1] * c->v2(&k) * c->v3(&k) * c->dv(&k);
            // for (int k=1;k<=Nv;++k)
            // rowsum+= twopi*e_mass*f0->get(&gatheredvector,&i,&j,&k)*c->v2(&k)*c->v2(&k)*c->dv(&k); // 4pi /2
            answer->values.set(i, j, rowsum);
            totalsum += rowsum;
        }
    }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
// the next is for quickeness - no update of grid, just values.
inline void get_Ue(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += twopi * e_mass * data[k - 1] * c->v2(&k) * c->v2(&k) * c->dv(&k);
            // for (int k=1;k<=Nv;++k)
            // rowsum+= twopi*e_mass*f0->get(&gatheredvector,&i,&j,&k)*c->v2(&k)*c->v2(&k)*c->dv(&k);// 4pi /2
            answer->values.set(i, j, rowsum);
        }
    }
    delete[] data;
}

inline void new_je(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    char dimension = 'x';
    switch (x1_int)
    {
    case 1:
        dimension = 'x';
        break;
    case 2:
        dimension = 'y';
        break;
    case 3:
        dimension = 'z';
        break;
    }
    char nametemp[10] = {'j', '_', dimension};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    int kstart = x1_int * Nv - 1; // if 1 - just skip f0, if 2 skip f0 and f1x1 etc.
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += fourthirdspi * e_charge * data[k + kstart] * c->v3(&k) * c->dv(&k);

            /*	for (int k=1;k<=Nv;++k)
            {
            rowsum+=fourthirdspi*e_charge*f1->get(&gatheredvector,&i,&j,&k,x1)*c->v3(&k)*c->dv(&k); // 4pi/3 e int f1v^3 dv*/
            answer->values.set(i, j, rowsum);
            //}
        }
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
inline void new_Te(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    strcpy(answer->name, "T_e");
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    struct IMPACT_Moment Ue, ne;
    new_Ue(v, f0, c, MPIc, &Ue);
    new_ne(v, f0, c, MPIc, &ne);
    struct IMPACT_Moment Te = Divide(&Ue, &ne);
    for (int i = 1; i <= c->Nx(); ++i)
        for (int j = 1; j <= c->Ny(); ++j)
            answer->values.set(i, j, Te.values.get(i, j) * 2.0 / 3.0);
    answer->average = 2.0 / 3.0 * Ue.average / ne.average;
}

inline void new_eta(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c, IMPACT_MPI_Config *MPIc,
                    struct IMPACT_Moment *answer)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int Nv = c->Nv();

    double *data;
    data = new double[c->cellpoints()];
    double rowsum, totalsum = 0.0;
    strcpy(answer->name, "eta");
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    double netemp;
    // f0 MOMENT
    //************************************************************
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            netemp = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += data[k - 1] * c->v2(&k) * c->v3(&k) * c->dv(&k);

            for (int k = 1; k <= Nv; ++k)
                netemp += data[k - 1] * c->dv(&k) * c->v2(&k); // get ne

            answer->values.set(i, j, 0.5 * Initial_Conditions::Z.get(&i, &j) * netemp / rowsum);
            totalsum += rowsum;
        }
    }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis

    delete[] data;
}

// - vector quants

inline void get_je(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    int x1_int = x1->get();
    int kstart = x1_int * Nv - 1; // if 1 - just skip f0, if 2 skip f0 and f1x1 etc.
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += fourthirdspi * e_charge * data[k + kstart] * c->v3(&k) * c->dv(&k);
            /*	for (int k=1;k<=Nv;++k)
            {
            rowsum+=fourthirdspi*e_charge*f1->get(&gatheredvector,&i,&j,&k,x1)*c->v3(&k)*c->dv(&k); // 4pi/3 e int f1v^3 dv*/
            answer->values.set(i, j, rowsum);
            //}
        }
    delete[] data;
}
inline void new_qT(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    char dimension = 'x';
    switch (x1_int)
    {
    case 1:
        dimension = 'x';
        break;
    case 2:
        dimension = 'y';
        break;
    case 3:
        dimension = 'z';
        break;
    }
    char nametemp[10] = {'q', '_', dimension};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();

    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    int kstart = x1_int * Nv - 1; // if 1 - just skip f0, if 2 skip f0 and f1x1 etc.
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += fourthirdspi * data[k + kstart] * c->v3(&k) * c->v2(&k) * c->dv(&k);
            answer->values.set(i, j, rowsum / 2.0);
        }
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
inline void new_VN(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Config *c, IMPACT_MPI_Config *MPIc,
                   struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int Nv = c->Nv();
    int x1_int = x1->get();
    char dimension = 'x';
    double *data;
    data = new double[c->cellpoints()];
    double rowsum, totalsum = 0.0;
    char nametemp[10] = {'V', '_', 'N', dimension};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    //double jtemp, ptemp;
    // f0 MOMENT
    //************************************************************
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += data[k - 1] * c->v2(&k) * c->v3(&k) * c->dv(&k);

            answer->values.set(i, j, rowsum);
            totalsum += rowsum;
        }
    }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
                        // f1 MOMENT -> M3
                        //************************************************************
    switch (x1_int)
    {
    case 1:
        dimension = 'x';
        break;
    case 2:
        dimension = 'y';
        break;
    case 3:
        dimension = 'z';
        break;
    }

    // double data2[c->cellpoints()];

    int kstart = x1_int * Nv - 1; // if 1 - just skip f0, if 2 skip f0 and f1x1 etc.
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            // jtemp = 0.0;
            // ptemp=1.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += data[k + kstart] * c->v3(&k) * c->v3(&k) * c->dv(&k);
            //-------------------------------------------------------------
            // get current
            // for (int k=1;k<=Nv;++k)
            // jtemp+=fourthirdspi*e_charge*data[k+kstart]*c->v3(&k)*c->dv(&k);
            //-------------------------------------------------------------

            // ptemp=Local_pe(v,f0,c,&i,&j);

            answer->values.set(i, j, (oneover6 * rowsum / (answer->values.get(i, j)))); //+jtemp/ptemp);
        }
    // GetGrid(c,answer); // gets x and y axis
    delete[] data;
}

inline void get_q(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int Nv = c->Nv();
    int x1_int = x1->get();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    int kstart = x1_int * Nv - 1; // if 1 - just skip f0, if 2 skip f0 and f1x1 etc.
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += fourthirdspi * data[k + kstart] * c->v3(&k) * c->v2(&k) * c->dv(&k);
            answer->values.set(i, j, rowsum / 2.0);
        }

    delete[] data;
}
// anisotropic pressure
inline void new_P(IMPACT_ParVec *v, IMPACT_Var *f2, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1, IMPACT_Dim *x2)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    int x2_int = x2->get();
    if (x1_int < 3)
        x1_int--;
    if (x2_int < 3)
        x2_int--; // similar to the system in IMPACT_Var
    char dimension[2] = {'e', 'e'};
    switch (x1_int + x2_int)
    {
    case 0:
        dimension[0] = 'x';
        dimension[1] = 'x';
        break;
    case 1:
        dimension[0] = 'x';
        dimension[1] = 'y';
        break;
    case 2:
        dimension[0] = 'y';
        dimension[1] = 'y';
        break;
    case 3:
        dimension[0] = 'x';
        dimension[1] = 'z';
        break;
    case 4:
        dimension[0] = 'y';
        dimension[1] = 'z';
        break;
    }
    char nametemp[10] = {'P', '_', dimension[0], '_', dimension[1]};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    int kstart = (x1_int + x2_int + c->Nf1() + 1) * Nv - 1; // skip f0 and f1
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += data[k + kstart] * c->v2(&k) * c->v2(&k) * c->dv(&k);
            answer->values.set(i, j, 2.0 / 5.0 * fourthirdspi * rowsum);
        }

    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}

// anisotropic pressure term in Ohm's law
inline void new_divP(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Var *f2, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1, IMPACT_Dim *x2)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    int x2_int = x2->get();
    if (x1_int < 3)
        x1_int--;
    if (x2_int < 3)
        x2_int--; // similar to the system in IMPACT_Var
    char dimension[2] = {'e', 'e'};
    switch (x1_int + x2_int)
    {
    case 0:
        dimension[0] = 'x';
        dimension[1] = 'x';
        break;
    case 1:
        dimension[0] = 'x';
        dimension[1] = 'y';
        break;
    case 2:
        dimension[0] = 'y';
        dimension[1] = 'y';
        break;
    case 3:
        dimension[0] = 'x';
        dimension[1] = 'z';
        break;
    case 4:
        dimension[0] = 'y';
        dimension[1] = 'z';
        break;
    }
    char nametemp[10] = {'D', '_', dimension[0], '_', dimension[1]};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Nv = c->Nv();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;

    int kstart = (x1_int + x2_int + c->Nf1() + 1) * Nv - 1; // skip f0 and f1
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += data[k + kstart] * c->v2(&k) * c->v2(&k) * c->dv(&k) * (c->v2(&k) * c->v(&k)); // M2(3)
            answer->values.set(i, j, (2.0 / 15.0) * rowsum);
        }

    GetGrid(c, answer); // gets x and y axis
                        //************************************************************
    // Now take gradient
    int iplus, iminus, jplus, jminus;

    double **answermethis;
    answermethis = new double *[Nx];
    for (int i = 0; i < Nx; ++i)
    {
        answermethis[i] = new double[Ny];
    }
    if ((*x1) < 3)
    {
        for (int i = 1; i <= Nx; ++i)
        {
            for (int j = 1; j <= Ny; ++j)
            {
                iplus = i + 1;
                iminus = i - 1;
                jplus = j + 1;
                jminus = j - 1;
                // temp_sten =(*O->ddxi(i,j,x1))
                IMPACT_stencil temp_sten(0.0, (0.5 / c->dx(&i)), -(0.5 / c->dx(&i)), (0.5 / c->dy(&j)), -(0.5 / c->dy(&j)));
                IMPACT_f2_bound(&temp_sten, Nx, Ny, &iplus, &iminus, &jplus, &jminus, x1, x2);

                // fix for periodic bounds
                if (iplus > Nx)
                {
                    iplus = 1;
                }
                if (iminus < 1)
                {
                    iminus = Nx;
                }
                if (jplus > Ny)
                {
                    jplus = 1;
                }
                if (jminus < 1)
                {
                    jminus = Ny;
                }

                if (*x1 == 1)
                {
                    answermethis[i - 1][j - 1] = answer->values.get(iplus, j) * temp_sten(1);
                    answermethis[i - 1][j - 1] += answer->values.get(iminus, j) * temp_sten(2);
                }
                if (*x1 == 2)
                {
                    answermethis[i - 1][j - 1] = answer->values.get(i, jplus) * temp_sten(3);
                    answermethis[i - 1][j - 1] += answer->values.get(i, jminus) * temp_sten(4);
                }
            } // i
        }     // j
    }
    // f0 MOMENT
    //************************************************************
    for (int i = 1; i <= Nx; ++i)
    {
        for (int j = 1; j <= Ny; ++j)
        {
            rowsum = 0.0;
            Gather_kstring(v, c, MPIc, &i, &j, data);
            for (int k = 1; k <= Nv; ++k)
                rowsum += data[k - 1] * c->v2(&k) * c->v3(&k) * c->dv(&k);

            answer->values.set(i, j, answermethis[i - 1][j - 1] / (2.0 * rowsum));
            // totalsum+=rowsum;
        }
    }
    // answer->average=totalsum/Nx/Ny;
    GetGrid(c, answer); // gets x and y axis
                        //************************************************************

    delete[] data;
}

inline void new_E(IMPACT_ParVec *v, IMPACT_Var *E, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    char dimension = 'x';
    switch (x1_int)
    {
    case 1:
        dimension = 'x';
        break;
    case 2:
        dimension = 'y';
        break;
    case 3:
        dimension = 'z';
        break;
    }
    char nametemp[10] = {'E', '_', dimension};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Epos = c->cellpoints() - c->NB() - c->NE() + x1_int - 1;

    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];

    double rowsum;
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            Gather_kstring(v, c, MPIc, &i, &j, data);
            rowsum = data[Epos]; // E->get(&gatheredvector,&i,&j,x1);
            answer->values.set(i, j, rowsum);
        }

    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
inline void get_E(IMPACT_ParVec *v, IMPACT_Var *E, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    int Epos = c->cellpoints() - c->NB() - c->NE() + x1_int - 1;
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            Gather_kstring(v, c, MPIc, &i, &j, data);
            rowsum = data[Epos]; // E->get(&gatheredvector,&i,&j,x1);
            answer->values.set(i, j, rowsum);
        }
    delete[] data;
}
inline void new_B(IMPACT_ParVec *v, IMPACT_Var *B, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    char dimension = 'x';
    switch (x1_int)
    {
    case 1:
        dimension = 'x';
        break;
    case 2:
        dimension = 'y';
        break;
    case 3:
        dimension = 'z';
        break;
    }
    char nametemp[10] = {'B', '_', dimension};
    strcpy(answer->name, nametemp);
    strncpy(answer->coords, IMPACT_Coords, Ncoords);
    answer->Nx = Nx;
    answer->Ny = Ny;
    answer->values.ChangeSize(Nx, Ny);
    int Bpos = c->cellpoints() - 4 + x1_int;
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum, totalsum = 0.0;
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            Gather_kstring(v, c, MPIc, &i, &j, data);
            rowsum = data[Bpos]; // E->get(&gatheredvector,&i,&j,x1);
            answer->values.set(i, j, rowsum);
            totalsum += rowsum;
        }
    answer->average = totalsum / Nx / Ny;
    GetGrid(c, answer); // gets x and y axis
    delete[] data;
}
inline void get_B(IMPACT_ParVec *v, IMPACT_Var *B, IMPACT_Config *c, IMPACT_MPI_Config *MPIc, struct IMPACT_Moment *answer, IMPACT_Dim *x1)
{
    int Nx = c->Nx();
    int Ny = c->Ny();
    int x1_int = x1->get();
    int Bpos = c->cellpoints() - 4 + x1_int;
    // IMPACT_ParVec gatheredvector=Gather(v,c,MPIc);
    double *data;
    data = new double[c->cellpoints()];
    double rowsum;
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            Gather_kstring(v, c, MPIc, &i, &j, data);
            rowsum = data[Bpos]; // E->get(&gatheredvector,&i,&j,x1);
            answer->values.set(i, j, rowsum);
        }
    delete[] data;
}

//________________________________________________________________________

// Local moments;
inline double Local_Te(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c,
                       int *i, int *j)
{
    int Nv = c->Nv();
    double rowsum = 0.0;
    double answer = 0.0;
    for (int k = 1; k <= Nv; ++k)
        rowsum += f0->get(v, i, j, &k) * c->v2(&k) * c->v2(&k) * c->dv(&k);

    answer = rowsum;
    rowsum = 0.0;

    for (int k = 1; k <= Nv; ++k)
        rowsum += f0->get(v, i, j, &k) * c->v2(&k) * c->dv(&k);

    answer = 2.0 / 3.0 * answer / rowsum;
    return answer;
}
inline double Local_ne(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c,
                       int *i, int *j)
{
    int Nv = c->Nv();
    double rowsum = 0.0;
    double answer = 0.0;

    for (int k = 1; k <= Nv; ++k)
        rowsum += f0->get(v, i, j, &k) * c->v2(&k) * c->dv(&k);

    answer = fourpi * rowsum;
    return answer;
}
inline double Local_pe(IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Config *c,
                       int *i, int *j)
{
    double answer = Local_Te(v, f0, c, i, j) * Local_ne(v, f0, c, i, j);
    return answer;
}

/*
  Local anisotropic part of pressure tensor
 */
inline double Local_Pi(IMPACT_ParVec *v, IMPACT_Var *f2, IMPACT_Config *c,
                       int *i, int *j, IMPACT_Dim *x1, IMPACT_Dim *x2)
{
    int Nv = c->Nv();
    double rowsum = 0.0;
    double answer = 0.0;
    rowsum = 0.0;
    for (int k = 1; k <= Nv; ++k)
        rowsum += f2->get(v, i, j, &k, x1, x2) * c->v2(&k) * c->v2(&k) * c->dv(&k);
    answer = 2.0 / 5.0 * fourthirdspi * rowsum;
    return answer;
}
inline double Local_je(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c,
                       int *i, int *j, IMPACT_Dim *x1)
{
    int Nv = c->Nv();
    double rowsum = 0.0;
    double answer = 0.0;
    for (int k = 1; k <= Nv; ++k)
        rowsum += f1->get(v, i, j, &k, x1) * c->v2(&k) * c->v(&k) * c->dv(&k);

    answer = e_charge * fourthirdspi * rowsum;
    return answer;
}
inline double Local_qT(IMPACT_ParVec *v, IMPACT_Var *f1, IMPACT_Config *c,
                       int *i, int *j, IMPACT_Dim *x1)
{
    int Nv = c->Nv();
    double rowsum = 0.0;
    double answer = 0.0;
    for (int k = 1; k <= Nv; ++k)
        rowsum += f1->get(v, i, j, &k, x1) * c->v2(&k) * c->v3(&k) * c->dv(&k);

    answer = fourthirdspi * rowsum * 0.5;
    return answer;
}
