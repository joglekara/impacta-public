/*
**********************************************************
Fills sparse matrix from  equations in IMPACTA code
Including boundary conditions and differential info

Version 1.0

AGRT

12/2/07
Notes:
20/2/07
Would be good to improve way matrix is packed - at the moment rudimentary
process of filling in row ParVec, then counting non zeros and
putting in sparse matrix
(see Matrix.h also)

22/2
"count" version counts number of non-zero elements

5/3
The count version should be eliminated at some point to make it
simpler to add in other equations - the outer equation form
can do the job easily enough, if modified.

8/3 Count made more efficient - so unlikely to b e able to eliminate
So much more efficient!! - in comparison, 0.9 s with new algorith,
I'm still waiting for the old one to finish under the same conditions!

9/4 - added switches to switch on and off different equations
these have required the addition of functions to ParSpa

12/4/07 - made slightly more efficient

15/4/07 - added IMPACT_Switch_Elements, which switches
off the diagonal elements in the solution matrix (for
switching off d/dt terms)

23/8/07 - Messages updated so that only the rank 0 processor counts.

25/9/07 - f1 equation RHS has grad f0 term added from gradn and T z

1/10/07 f0 equation RHS has grad f1 term added for self consistancy

6/10/07 - changed so that istart and iend cells are not at boundaries
        - trying to solve boundary problems!
Note both these terms are from formula with no inertia terms
and ignoring f2
11/2/08 - minor change so that f2 depends on f0 through IB heating

13/3/08 - curl b term corrected on n side of equation
2/7/08 - Fixed boundaries implimented so that if chosen, don't
evolve relevent boundaries

jan 2011 - fixed bounds adjusted so that a bit more correct -
now boundary taken in 1 so effectively have ghost cells
**********************************************************
*/

// Check fixed boundaries,
inline int check_fixed(IMPACT_Config *c, int *boundx, int *boundy)
{
    int answer = 0;

    if (*boundx > -1)
    {
        answer += IMPACT_Boundaries::fixed_any[*boundx];
    }
    if (*boundy > -1)
    {
        answer += IMPACT_Boundaries::fixed_any[*boundy];
    }
    return answer;
}
inline void equations_outer(IMPACT_Config *config1, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_ParSpa *S, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Var *f2, IMPACT_Var *f3, IMPACT_Var *E, IMPACT_Var *B, int *i, int *j, int *kmin, int *kmax);

// The code for equations inside boundaries
inline void equations_inner(IMPACT_Config *config1, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_ParSpa *S, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Var *f2, IMPACT_Var *f3, IMPACT_Var *E, IMPACT_Var *B, int *i, int *j)
{
    int bx = -1, by = -1;
    if (*i == 2)
        bx = 0;
    if (*i == (config1->Nx() - 1))
        bx = 1;
    if (*j == 2)
        by = 2;
    if (*j == (config1->Ny() - 1))
        by = 3;

    if (!check_fixed(config1, &bx, &by))
    {
        //___________________________________________________________

        if (equation_switches::f0_equation_on)
            for (int k = 2; k < config1->Nv(); ++k)
            {
                f0equation_inner(config1, vlagged, vtemp, O, f0, E, B, f1, i, j, &k);
                f0equation_Pack_inner(vtemp, S, config1, f0, f1, E, B, i, j, &k);
                f0equation_clean_inner(vtemp, config1, O, f0, E, B, f1, i, j, &k);
            }
        else
            for (int k = 2; k < config1->Nv(); ++k)
                S->nullequation(config1, 0, f0->getrow(i, j, &k));

        //____________________________________________________________

        // f1  equation now
        if (config1->Nf1() > 0)
            for (IMPACT_Dim x1 = 1; x1 <= config1->Nf1(); ++x1)
            {
                if (equation_switches::f1_equation_on)
                    for (int k = 2; k < config1->Nv(); ++k)
                    {
                        f1equation_inner(config1, vlagged, vtemp, O, f0, E, B, f1, f2, i, j, &k, &x1);
                        f1equation_Pack_inner(vtemp, S, config1, f0, f1, f2, E, B, i, j, &k, &x1);
                        // This packs the completed ParVec
                        f1equation_clean_inner(vtemp, config1, O, f0, E, B, f1, f2, i, j, &k, &x1);
                    }
                else
                    for (int k = 2; k < config1->Nv(); ++k)
                        S->nullequation(config1, 1, f1->getrow(i, j, &k, &x1));
            }
        // f2 equation.....
        if (config1->Nf2() > 0)
            for (IMPACT_Dim x1 = 1; x1 <= config1->N3f2(); ++x1)
                for (IMPACT_Dim x2 = x1.get(); x2 <= config1->N3f2(); ++x2)
                    if (x1.get() + x2.get() < 6)
                    {
                        if (equation_switches::f2_equation_on)
                            for (int k = 2; k < config1->Nv(); ++k)
                            {
                                f2equation_inner(config1, vlagged, vtemp, O, f0, f1, f2, f3, E, B,
                                                 i, j, &k, &x1, &x2);
                                f2equation_Pack_inner(vtemp, S, config1, f0, f1, f2, f3, E, B,
                                                      i, j, &k, &x1, &x2);
                                // This packs the completed ParVec
                                f2equation_clean_inner(vtemp, config1, O, E, B, f0, f1, f2, f3,
                                                       i, j, &k, &x1, &x2);
                            }
                        else
                            for (int k = 2; k < config1->Nv(); ++k)
                                S->nullequation(config1, 2, f2->getrow(i, j, &k, &x1, &x2));
                    }
        // E equation now
        if (config1->NE() > 0)
            for (IMPACT_Dim x1 = 1; x1 <= config1->NE(); ++x1)
            {
                if (equation_switches::E_equation_on)
                {
                    Eequation_inner(config1, vlagged, vtemp, O, E, B, f1, i, j, &x1);
                    Eequation_Pack_inner(vtemp, S, config1, f1, E, B, i, j, &x1);
                    // This packs the completed ParVec
                    vtemp->reset();
                    Eequation_clean_inner(vtemp, config1, O, E, B, f1, i, j, &x1);
                }
                else
                    S->nullequation(config1, 4, E->getrow(i, j, &x1));
            }
        //___________________________________________________________

        // finally B equation.
        if (config1->NB() > 0)
            for (IMPACT_Dim x1 = 3; x1 > 3 - config1->NB(); --x1)
            {
                if (equation_switches::B_equation_on)
                {
                    Bequation_inner(config1, vlagged, vtemp, O, E, B, i, j, &x1);
                    Bequation_Pack_inner(vtemp, S, config1, E, B, i, j, &x1);
                    // S->Pack(vtemp,rownumber); //This packs the completed ParVec
                    Bequation_clean_inner(vtemp, config1, O, E, B, i, j, &x1);
                }
                else
                    S->nullequation(config1, 5, B->getrow(i, j, &x1));
            }
        //____________________________________________________________
    }
    else
    {
        int kmin = 1;
        int kmax = config1->Nv();
        equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3, E, B, i, j, &kmin, &kmax);
    }
}

inline void equations_outer(IMPACT_Config *config1, IMPACT_ParVec *vlagged, IMPACT_ParVec *vtemp, IMPACT_ParSpa *S, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Var *f2, IMPACT_Var *f3, IMPACT_Var *E, IMPACT_Var *B, int *i, int *j, int *kmin, int *kmax)
{
    int bx = -1, by = -1;
    if (*i == 1)
        bx = 0;
    if (*i == config1->Nx())
        bx = 1;
    if (*j == 1)
        by = 2;
    if (*j == config1->Ny())
        by = 3;

    //___________________________________________________________
    if (equation_switches::f0_equation_on && !check_fixed(config1, &bx, &by))
        for (int k = *kmin; k <= *kmax; ++k)
        {
            f0equation_outer(config1, vlagged, vtemp, O, f0, E, B, f1, i, j, &k);
            // S->Pack(vtemp,rownumber); //This packs the completed ParVec
            f0equation_Pack_outer(vtemp, S, config1, f0, f1, E, B, i, j, &k);
            f0equation_clean_outer(vtemp, config1, O, f0, E, B, f1, i, j, &k);
        }
    else
        for (int k = *kmin; k <= *kmax; ++k)
            S->nullequation(config1, 0, f0->getrow(i, j, &k));
    //____________________________________________________________

    // f1 and E equations now
    if (config1->Nf1() > 0)
        for (IMPACT_Dim x1 = 1; x1 <= config1->Nf1(); ++x1)
        {
            if (equation_switches::f1_equation_on && !check_fixed(config1, &bx, &by))
                for (int k = *kmin; k <= *kmax; ++k)
                {
                    f1equation_outer(config1, vlagged, vtemp, O, f0, E, B, f1, f2, i, j, &k, &x1);
                    f1equation_Pack_outer(vtemp, S, config1, f0, f1, f2, E, B, i, j, &k, &x1);
                    // This packs the completed ParVec
                    f1equation_clean_outer(vtemp, config1, O, f0, E, B, f1, f2, i, j, &k, &x1);
                }
            else
                for (int k = *kmin; k <= *kmax; ++k)
                    S->nullequation(config1, 1, f1->getrow(i, j, &k, &x1));
        }
    // f2 equation.....
    if (config1->Nf2() > 0)
        for (IMPACT_Dim x1 = 1; x1 <= config1->N3f2(); ++x1)
            for (IMPACT_Dim x2 = x1.get(); x2 <= config1->N3f2(); ++x2)
                if (x1.get() + x2.get() < 6)
                {
                    if (equation_switches::f2_equation_on && !check_fixed(config1, &bx, &by))
                        for (int k = *kmin; k <= *kmax; ++k)
                        {
                            f2equation_outer(config1, vlagged, vtemp, O, f0, f1, f2, f3, E, B,
                                             i, j, &k, &x1, &x2);
                            f2equation_Pack_outer(vtemp, S, config1, f0, f1, f2, f3, E, B,
                                                  i, j, &k, &x1, &x2);
                            // This packs the completed ParVec
                            f2equation_clean_outer(vtemp, config1, O, E, B, f0, f1, f2, f3,
                                                   i, j, &k, &x1, &x2);
                        }
                    else
                        for (int k = *kmin; k <= *kmax; ++k)
                            S->nullequation(config1, 2, f2->getrow(i, j, &k, &x1, &x2));
                }

    if (config1->NE() > 0)
        for (IMPACT_Dim x1 = 1; x1 <= config1->NE(); ++x1)
        {
            if (equation_switches::E_equation_on && !check_fixed(config1, &bx, &by))
            {
                Eequation_outer(config1, vlagged, vtemp, O, E, B, f1, i, j, &x1);
                Eequation_Pack_outer(vtemp, S, config1, f1, E, B, i, j, &x1);
                // This packs the completed ParVec
                vtemp->reset();
                Eequation_clean_outer(vtemp, config1, O, E, B, f1, i, j, &x1);
            }
            else
                S->nullequation(config1, 4, E->getrow(i, j, &x1));
        }
    //___________________________________________________________

    // finally B equation.
    if (config1->NB() > 0)
        for (IMPACT_Dim x1 = 3; x1 > 3 - config1->NB(); --x1)
        {
            if (equation_switches::B_equation_on && !check_fixed(config1, &bx, &by))
            {
                Bequation_outer(config1, vlagged, vtemp, O, E, B, i, j, &x1);
                Bequation_Pack_outer(vtemp, S, config1, E, B, i, j, &x1);
                // This packs the completed ParVec
                Bequation_clean_outer(vtemp, config1, O, E, B, i, j, &x1);
            }
            else
                S->nullequation(config1, 5, B->getrow(i, j, &x1));
        }
    //____________________________________________________________
}

inline void IMPACT_Form_Sparse(IMPACT_Config *config1, IMPACT_MPI_Config *MPIc, IMPACT_ParVec *vlagged, IMPACT_ParSpa *S, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Var *f2, IMPACT_Var *f3, IMPACT_Var *E, IMPACT_Var *B)
{
    // vlagged - ParVec containing lagged data, S- sparse matrix
    /*
      This function packs the parts of the sparse
      matrix S which are f0 equations
    */
    int rank = MPIc->rank(), numprocs = MPIc->size();
    double itemp;
    int checker = 1;
    int istart = MPIc->istart();
    int iend = MPIc->iend();
    std::ostringstream Imessage("");
    IMPACT_ParVec *vtemp;                                        // temporary store of row elements
    vtemp = new IMPACT_ParVec(vlagged->start(), vlagged->end()); // i.e. as long as the side of the matrix

    vtemp->reset(); // This cleans up the temporary vectro for next use
    /* _________________________________________________________
       INNER CELLS
       _________________________________________________________
       simple - just doesn't test whether i-1 etc is outside area
     */
    clock_t timestart_inner, timeend_inner; // for timing matrix count.
    timestart_inner = clock();
    MPI_Barrier(MPI_COMM_WORLD);
    if (iend - istart > 3 && !rank)
        Imessage << ENDFORMAT << "\nIMPACT: Building Sparse Matrix";
    if (iend - istart > 3 && numprocs == 1)
        Imessage << " (Inner loop)\n|<---------------->|";
    if (!rank)
        Imessage << '\n';
    std::cout << Imessage.str();
    int mystart = istart, myend = iend;
    if (rank == 0)
        ++mystart;
    if (rank == numprocs - 1)
        --myend;
    if (config1->Ny() > 2 && config1->Nx() > 2 && config1->Nv() > 2)
        for (int i = mystart; i <= myend; ++i)
        {
            // No communication needed for these cells
            for (int j = 2; j < config1->Ny(); ++j)
            {

                equations_inner(config1, vlagged, vtemp, S, O, f0, f1, f2, f3, E, B, &i, &j);
            }
            itemp = (double)(i - istart) / (iend - istart - 1) * 10.1;
            if ((int)itemp == checker && numprocs == 1)
            {
                std::cout << BRED << "[]" << ENDFORMAT;
                ++checker;
            }
        }

    timeend_inner = clock();

    /* _________________________________________________________
       BOUNDARY CELLS
       _________________________________________________________
       Now we do boundary cells: These form a 3D box shell around
       the inner cells. We will do first each plane A (6 planes), then the
       wireframe line B (12 lines). (i.e. initially missing out the edge of
       each plane) Then finally do the corner cells C (6 of them)

               ______
          |\	    \
          | \ _____\ <- C
          |  |      |
           \ |  A   |<- B
            \|______|

            update 19/2/07
            now the planes at i=1 and Nx are large (go to edges)
            then j planes goes stops at i=2
            FINALLY k planes stop at i=2 and j=2
    */
    /*
        first the i boundary planes
        5/10/07 - Now not used!
       */
    clock_t timestart_outer, timeend_outer; // for timing matrix count.
    int kmin, kmax;
    if (iend - istart > 3 && numprocs == 1)
        std::cout << ENDFORMAT << "\nIMPACT: Building Sparse Matrix (Outer loop)\n|<-->|\n";
    else if (numprocs == 1)
        std::cout << ENDFORMAT << "\nIMPACT: Building Sparse Matrix\n|<-->|\n";
    timestart_outer = clock();
    if (config1->Ny() > 2 && config1->Nx() > 2 && config1->Nv() > 2)
    {
        if (rank == 0 || rank == numprocs - 1)
            for (int j = 2; j < config1->Ny(); ++j)
            {
                kmin = 1;
                kmax = config1->Nv();
                if (rank == 0)
                {
                    int i = istart;
                    equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                                    E, B, &i, &j, &kmin, &kmax);
                }
                if (rank == numprocs - 1)
                {
                    int i = iend;
                    equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                                    E, B, &i, &j, &kmin, &kmax);
                }
            }
        if (numprocs == 1)
            std::cout << BRED << "[]" << ENDFORMAT;
        /*
      now the j boundary planes
        */

        for (int i = istart; i <= iend; ++i)
        {
            kmin = 1;
            kmax = config1->Nv();
            int j = 1;
            equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                            E, B, &i, &j, &kmin, &kmax);
            j = config1->Ny();
            equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                            E, B, &i, &j, &kmin, &kmax);
        }
        if (numprocs == 1)
            std::cout << BRED << "[]" << ENDFORMAT;
        /*
      now the k boundary planes
        */

        for (int i = istart; i <= iend; ++i)
            for (int j = 2; j < config1->Ny(); ++j)
            {
                kmax = kmin = 1;
                equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                                E, B, &i, &j, &kmin, &kmax);
                kmin = kmax = config1->Nv();
                equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                                E, B, &i, &j, &kmin, &kmax);
            }
        if (numprocs == 1)
            std::cout << BRED << "[]" << ENDFORMAT;
    }

    if (config1->Ny() < 3 || config1->Nx() < 3 || config1->Nv() < 3)
    {
        for (int i = istart; i <= iend; ++i)
            for (int j = 1; j <= config1->Ny(); ++j)
            {
                kmin = 1;
                kmax = config1->Nv();
                equations_outer(config1, vlagged, vtemp, S, O, f0, f1, f2, f3,
                                E, B, &i, &j, &kmin, &kmax);
            }
    }
    delete vtemp;
    timeend_outer = clock();
    // MESSAGES
    std::ostringstream Imessage2;
    if (if_time_messages)
    {
        if (numprocs == 1)
        {
            std::cout << std::endl
                      << "IMPACT: Matrix Build - Inner Loop took " << IMPACT_GetTime(double(timeend_inner - timestart_inner) / CLOCKS_PER_SEC) << std::endl;
            std::cout << "                       Outer Loop took " << IMPACT_GetTime(double(timeend_outer - timestart_outer) / CLOCKS_PER_SEC) << std::endl;
            std::cout << ULINE;
        }
        
        Imessage2 << BYELLOW << "\nIMPACT: Processor " << rank << " of " << numprocs << " Matrix Build - Total time = " << IMPACT_GetTime(double(timeend_outer - timestart_inner) / CLOCKS_PER_SEC) << ENDFORMAT;
        if (!rank)
            Imessage2 << '\n'
                      << ULINE;
    }
    if (equation_switches::Cee0_on && !rank)
    {
        Imessage2 << "\n(Average Chang-Cooper delta iterations = " << IMPACT_Diagnostics::total_delta_its / IMPACT_Diagnostics::delta_times << ")\n";
    }
    std::cout << Imessage2.str();
    IMPACT_Diagnostics::total_delta_its = IMPACT_Diagnostics::delta_times = 0.0;
}

/*
This version of the function counts the non-zero elements to initialize the
sparse matrix
*/
inline void IMPACT_Form_Sparse_count(IMPACT_Config *config1, IMPACT_MPI_Config *MPIc, IMPACT_ParVec *vlagged, IMPACT_ParSpa *S, IMPACT_StenOps *O, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Var *f2, IMPACT_Var *f3, IMPACT_Var *E, IMPACT_Var *B)
{
    clock_t timestart, timeend; // for timing matrix count.
    timestart = clock();
    IMPACT_ParVec *vtemp;                                        // temporary store of row elements
    vtemp = new IMPACT_ParVec(vlagged->start(), vlagged->end()); // i.e. as long as the side of the matrix
    vtemp->reset();
    double itemp;
    int checker = 1;
    int istart = MPIc->istart();
    int iend = MPIc->iend();
    int rank = MPIc->rank();
    int numprocs = MPIc->size();
    std::ostringstream Imessage("");
    // As this only needs to
    // int rownumber;
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank)
        Imessage << ENDFORMAT << "\nIMPACT: Counting maximum elements\n";
    if (numprocs == 1)
        Imessage << "|<---------------->|\n";
    std::cout << Imessage.str();
    for (int i = istart; i <= iend; ++i)
    {
        for (int j = 1; j <= config1->Ny(); ++j)
        {
            if (equation_switches::f0_equation_on)
            {
                // IMPACT_Update_Cee0_BC(config1,vlagged,f0,&i,&j);
                for (int k = 1; k <= config1->Nv(); ++k)
                {
                    //___________________________________________________________
                    // f0 equation
                    f0equation_outer(config1, vlagged, vtemp, O, f0, E, B, f1, &i, &j, &k);
                    f0equation_count(vtemp, S, config1, f0, E, B, f1, &i, &j, &k);
                    f0equation_clean_outer(vtemp, config1, O, f0, E, B, f1, &i, &j, &k);
                    //____________________________________________________________
                }
            }
            else
                for (int k = 1; k <= config1->Nv(); ++k)
                    S->count_nullequation(config1, f0->getrow(&i, &j, &k));
            // f1 and E equations now
            if (config1->Nf1() > 0)
                for (IMPACT_Dim x1 = 1; x1 <= config1->Nf1(); ++x1)
                {
                    if (equation_switches::f1_equation_on)
                        for (int k = 1; k <= config1->Nv(); ++k)
                        {
                            f1equation_outer(config1, vlagged, vtemp, O, f0, E, B, f1, f2, &i, &j, &k, &x1);
                            f1equation_count(vtemp, S, config1, f0, E, B, f1, f2, &i, &j, &k, &x1);
                            f1equation_clean_outer(vtemp, config1, O, f0, E, B, f1, f2, &i, &j, &k, &x1);
                        }
                    else
                        for (int k = 1; k <= config1->Nv(); ++k)
                            S->count_nullequation(config1, f1->getrow(&i, &j, &k, &x1));
                }

            // f2 equation
            if (config1->Nf2() > 0)
                for (IMPACT_Dim x1 = 1; x1 <= config1->N3f2(); ++x1)
                    for (IMPACT_Dim x2 = x1.get(); x2 <= config1->N3f2(); ++x2)
                        if (x1.get() + x2.get() < 6)
                        {
                            if (equation_switches::f2_equation_on)
                                for (int k = 1; k <= config1->Nv(); ++k)
                                {
                                    f2equation_outer(config1, vlagged, vtemp, O, f0, f1, f2, f3, E, B,
                                                     &i, &j, &k, &x1, &x2);
                                    f2equation_count(vtemp, S, config1, f0, f1, E, B, f2, f3, &i, &j, &k,
                                                     &x1, &x2);
                                    f2equation_clean_outer(vtemp, config1, O, E, B, f0, f1, f2, f3,
                                                           &i, &j, &k, &x1, &x2);
                                }
                            else
                                for (int k = 1; k <= config1->Nv(); ++k)
                                    S->count_nullequation(config1, f2->getrow(&i, &j, &k, &x1, &x2));
                        }

            if (config1->NE() > 0)
                for (IMPACT_Dim x1 = 1; x1 <= config1->NE(); ++x1)
                {
                    if (equation_switches::E_equation_on)
                    {
                        Eequation_outer(config1, vlagged, vtemp, O, E, B, f1, &i, &j, &x1);
                        Eequation_count(vtemp, S, config1, E, B, f1, &i, &j, &x1);
                        Eequation_clean_outer(vtemp, config1, O, E, B, f1, &i, &j, &x1);
                        vtemp->reset();
                    }
                    else
                        S->count_nullequation(config1, E->getrow(&i, &j, &x1));
                }

            //___________________________________________________________

            // finally B equation.
            if (config1->NB() > 0)
                for (IMPACT_Dim x1 = 3; x1 > 3 - config1->NB(); --x1)
                {
                    if (equation_switches::B_equation_on)
                    {
                        Bequation_outer(config1, vlagged, vtemp, O, E, B, &i, &j, &x1);
                        Bequation_count(vtemp, S, config1, E, B, &i, &j, &x1);
                        Bequation_clean_outer(vtemp, config1, O, E, B, &i, &j, &x1);
                    }
                    else
                        S->count_nullequation(config1, B->getrow(&i, &j, &x1));
                }
            //____________________________________________________________
        }
        itemp = (double)(i - istart) / (iend - istart) * 10.1;
        if ((int)itemp == checker && numprocs == 1)
        {
            std::cout << BRED << "[]" << ENDFORMAT;
            ++checker;
        }
    }
    delete vtemp;
    S->Resize(); // THIS IS  IMPORTANT. Resizes the S vector and sorts out rowvector
    timeend = clock();

    // MESSAGES
    std::ostringstream Imessage2("");
    Imessage2 << ENDFORMAT << "\nIMPACT: Processor " << rank << " of " << numprocs << " - Element Count took " << IMPACT_GetTime(double(timeend - timestart) / CLOCKS_PER_SEC) << "  ";
    Imessage2 << "-- Sparse Matrix memory reserved for " << S->getSl() << " elements";
    if (!rank)
        Imessage2 << '\n';
    std::cout << Imessage2.str();
}
