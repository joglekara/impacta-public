/*
**********************************************************
The Matrix solver - takes an IMPACT_Sparse and converts
to a PETSc AIJ sparse. Solves locally and returns vector.

Version 2.9

AGRT

15/3/07

15/4/07 - NB NB NB!!! vin is RHS (n) vector, vout is LHS (n+1)

25/10/07 - Update! Uses explicit soln as possible start solution

29/10/07 - Explicit soln needs rethinking but importantly the
rtol value for the solver previously = rtol(specified)/N
This was obviously unnecessary and also probably doesn't
helpt with convergence very much

1/2/08  Major Development!!!!!
----------------------------------------------------------

New iterative method of solving matrix. For scalar n,
equation

Ax=b
can be solved in the following way for x

(A+I)x=c

c = Ax+x = b+x

A+I = D

if we solve
D^-1 c = x

We can get x, but c=b+x. Can converge on solution with

x_(i+1) = D^-1 (b+x_i)

This is now implemented.
----------------------------------------------------------
**********************************************************
*/
static char help[] = "Matrix solver which takes IMPACT_Sparse and IMPACT_Vector forms and converts them to Petsc objects before solving the matrix equation";


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Input arguments for matrix solve:
    argc - number of commandline flags
    args - array of flags (strings)
    Ain - input Sparse matrix, A[0][Nz], A[1][N+1], A[2][Nz]
    bin - solution vector
    uin - vector to solve via u=(A)^(-1)b
    Nin - number of rows of Ain
         Nzin - number of non-zero entries in Ain - no! agrt291106

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*
Yale Sparse Matrix Format:

Matrix Aij is represented by 3 vectors:
A - The nonzero values in a long column
I - The position of the first entry of row i in vector A
J - The column values for the entries in vector A
Then to extract Aij from the three vectors, for value Aij
cycle over rows i, then cycle from j=I(i) to I(i+1)
Then extract A(j) and J(j)
*/

inline int MatrixSolve(IMPACT_Config *c, IMPACT_MPI_Config *M, IMPACT_ParSpa *IMPACT_S, IMPACT_ParVec *IMPACT_vin, IMPACT_ParVec *IMPACT_vout, int argc, char **args, int *lagged)
{

    Vec PETSC_vin, PETSC_vout;                 // Petsc Vectors
    /* approx solution, RHS, exact solution */ // no x now agrt251106
    Mat PETSC_S;                               /* PETSC Sparse matrix */
    KSP ksp;                                   /* linear solver context */
    PetscInt Istart, Iend, N, n, its;          // N is size of Matrix n is local length
    // its - iteration numebr
    PetscErrorCode ierr;
    PetscLogDouble time1, time2; // for recording the time!
    PC pc;                       // preconditioner
    // PetscScalar    v; //reusable double for inserting values
    // int *d_nnz,*o_nnz; //These will be arrays for preallocation of memory
    PetscScalar iteratematrixscalar = zerotolerance::iterate_matrix;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    int rank, size;
    int badreturn = 0;
    // Step 1 - set relevent values from MPI_Config:

    rank = M->rank();
    size = M->size();
    Istart = M->start() - 1;
    Iend = M->end() - 1;
    N = c->totalpoints();
    n = M->N();
    if (rank == 0)
        std::cout << "\n";
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Step 2 - initialize Matrix solver:

    PetscInitialize(&argc, &args, (char *)0, help);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           Compute the matrix and right-hand-side vector that define
           the linear system, Ax = b.
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
       (Petsc notes...)
       Create parallel matrix, specifying only its global dimensions.
       When using MatCreate(), the matrix format can be specified at
       runtime. Also, the parallel partitioning of the matrix is
       determined by PETSc at runtime.

       Performance tuning note:  For problems of substantial size,
       preallocation of matrix memory is crucial for attaining good
       performance. See the matrix chapter of the users manual for details.
    */
    // Step 3 - preallocate memory
    MPI_Barrier(MPI_COMM_WORLD);
    // We will first use the  numinrows only -> can be improved later:
    PetscInt *preall;
    preall = new PetscInt[n];

    PetscInt columntemp; // needed because petsc goes from 0 to n-1
    double IfZeroRow;    // needed if whole row is zero - set diag to 1
    PetscScalar one = 1.0;

    for (int i = 0; i < n; ++i)
        preall[i] = IMPACT_S->getrow(i + Istart + 2) - IMPACT_S->getrow(i + Istart + 1);
startofmatrixsolve:
    MatCreateAIJ(PETSC_COMM_WORLD, n, n, N, N, 0, preall, 0, preall, &PETSC_S);
    // MatCreateMPIAIJ(PETSC_COMM_WORLD,n,n,N,N,0,preall,0,preall,&PETSC_S);
    ierr = MatSetOption(PETSC_S, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);
    CHKERRQ(ierr); // PETSC TRUE added for PETSC 3.0.0

    //____________________________________________________________________________
    /*
   Step 4 -Insert Matrix elements from Yale Sparse Matrix Format
    Remember IMPACT is from 1 to N whereas PETSC is from 0 to N-1

   //_____________________________________________________________________________
   */
    ierr = PetscTime(&time1);
    CHKERRQ(ierr); // get time
    PetscScalar my_zero = 0.0;
    for (int rowcount = Istart; rowcount <= Iend; ++rowcount)
    {
        // AGRT 2022 - first set diagonal entry to zero explicitly:
        ierr = MatSetValues(PETSC_S, 1, &rowcount, 1, &rowcount,
                                &my_zero, INSERT_VALUES);
        IfZeroRow = 0.0;
        for (int subrowcount = IMPACT_S->getrow(rowcount + 1);
             subrowcount < IMPACT_S->getrow(rowcount + 2); ++subrowcount)
        {
            IfZeroRow += IMPACT_S->getval(subrowcount);
            columntemp = IMPACT_S->getcol(subrowcount) - 1;
            ierr = MatSetValues(PETSC_S, 1, &rowcount, 1, &columntemp,
                                IMPACT_S->getval_add(subrowcount),
                                INSERT_VALUES);
            CHKERRQ(ierr);
        }
        if (IfZeroRow == 0.0)
            ierr = MatSetValues(PETSC_S, 1, &rowcount, 1, &rowcount,
                                &one, INSERT_VALUES);
    }

    /*
       Assemble matrix, using the 2-step process:
         MatAssemblyBegin(), MatAssemblyEnd()
       Computations can be done while messages are in transition
       by placing code between these two statements.
    */

    ierr = MatAssemblyBegin(PETSC_S, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(PETSC_S, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = PetscTime(&time2);
    CHKERRQ(ierr);
    std::ostringstream Imess1("");
    Imess1 << "\nRank " << rank << ": Sparse matrix transfer IMPACT->PETSC took - " << IMPACT_GetTime(time2 - time1);

    if (if_dump_switches::view_matrix)
        MatView(PETSC_S, PETSC_VIEWER_DRAW_WORLD);

    /*
       Create parallel vectors.
        - We form 1 vector from scratch and then duplicate as needed.
        - When using VecCreate(), VecSetSizes and VecSetFromOptions()
          in this example, we specify only the
          vector's global dimension; the parallel partitioning is determined
          at runtime.
        - When solving a linear system, the vectors and matrices MUST
          be partitioned accordingly.  PETSc automatically generates
          appropriately partitioned matrices and vectors when MatCreate()
          and VecCreate() are used with the same communicator.
        - The user can alternatively specify the local vector and matrix
          dimensions when more sophisticated partitioning is needed
          (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
          below).
    */
    // Step 5 - create vector from impact vector

    ierr = PetscTime(&time1);
    CHKERRQ(ierr); // get time
    VecCreateMPI(MPI_COMM_WORLD, n, N, &PETSC_vin);

    ierr = VecDuplicate(PETSC_vin, &PETSC_vout);
    CHKERRQ(ierr);

    // Insert values for vector PETSC_vin
    for (int i = Istart; i <= Iend; i++)
    {
        ierr = VecSetValues(PETSC_vin, 1, &i, IMPACT_vin->Get_add(i + 1), INSERT_VALUES);
        CHKERRQ(ierr);
        ierr = VecSetValues(PETSC_vout, 1, &i, IMPACT_vout->Get_add(i + 1), INSERT_VALUES);
        CHKERRQ(ierr);
    }
    VecAssemblyBegin(PETSC_vin);
    VecAssemblyBegin(PETSC_vout);
    VecAssemblyEnd(PETSC_vin);
    VecAssemblyEnd(PETSC_vout);
    ierr = PetscTime(&time2);
    CHKERRQ(ierr);
    Imess1 << "\nRank " << rank << ": Vector transfer IMPACT->PETSC took - " << IMPACT_GetTime(time2 - time1);
    if (!rank)
        Imess1 << '\n';
    if (if_time_messages) 
        std::cout << Imess1.str();
    MPI_Barrier(MPI_COMM_WORLD);

    /*  _-----------------------------------------------------------------
        step 5b - If user requests, find explicit solution of equation
        to act as rhs vector as a preconditioner.

        We do two half time steps to get soln
        -> Didn't work very well

        Try linear soln as start point:
        Single step explicit:
    */
    if (zerotolerance::explicit_PC)
    {
        PetscScalar idt = c->idt(), dt = c->dt(), Nsteps = 1.0;

        // These two operations are equivalent to M = 2/dtI - M;
        ierr = MatScale(PETSC_S, -1.0);
        CHKERRQ(ierr);
        ierr = MatShift(PETSC_S, idt);
        CHKERRQ(ierr); // gets rid of diagonal elem
        ierr = MatScale(PETSC_S, dt / Nsteps * 0.5);
        CHKERRQ(ierr);
        ierr = MatShift(PETSC_S, 1.0 / Nsteps);
        CHKERRQ(ierr); // adds diagonal

        // Now get solution M v_in = v_out
        ierr = MatMult(PETSC_S, PETSC_vin, PETSC_vout);
        CHKERRQ(ierr);
        /*
        Vec PETSC_vout_temp;
        ierr = VecDuplicate(PETSC_vout,&PETSC_vout_temp);CHKERRQ(ierr);

        for (int step =1;step<Nsteps;++step)
          {
            VecCopy(PETSC_vout,PETSC_vout_temp);
            ierr = MatMult(PETSC_S,PETSC_vout_temp,PETSC_vout);CHKERRQ(ierr);
          }
        VecDestroy(PETSC_vout_temp);
        */
        // now revert Matrix to original form:
        ierr = MatShift(PETSC_S, -1.0 / Nsteps);
        CHKERRQ(ierr); // adds diagonal
        ierr = MatScale(PETSC_S, idt * Nsteps * 2.0);
        CHKERRQ(ierr);
        ierr = MatShift(PETSC_S, -idt);
        CHKERRQ(ierr);
        ierr = MatScale(PETSC_S, -1.0);
        CHKERRQ(ierr);

        // Fiinally multiply by dt
        // ierr = VecScale(PETSC_vout,dt);CHKERRQ(ierr);
    }

    /*  _-----------------------------------------------------------------
       step 5b - If user requests, iterate towards solution using previous
       close solution
   */

    if (iteratematrixscalar > 0.0)
    {
        iteratematrixscalar = zerotolerance::iterate_matrix;
        if (!rank)
            std::cout << "\nIterating matrix:\nShifting matrix elements...";
        ierr = MatShift(PETSC_S, iteratematrixscalar);
        CHKERRQ(ierr);
        if (!rank)
            std::cout << "Finished\nMultiplying through...";
        ierr = VecAXPY(PETSC_vin, iteratematrixscalar, PETSC_vout);
        if (!rank)
            std::cout << "Finished\n";
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        step 6 - Create the linear solver and set various options
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
       Create linear solver context
    */
    ierr = PetscTime(&time1);
    CHKERRQ(ierr); // get time
    ierr = KSPCreate(MPI_COMM_WORLD, &ksp);
    CHKERRQ(ierr);

    KSPSetType(ksp, KSPBCGS); // Set method of ksp -biconjugate gradient stabilized
    /*
       Set operators. Here the matrix that defines the linear system
       also serves as the preconditioning matrix.
    */
    ierr = KSPSetOperators(ksp, PETSC_S, PETSC_S); //,DIFFERENT_NONZERO_PATTERN);
    CHKERRQ(ierr);
    /*
       Set linear solver defaults for this problem (optional).
       - By extracting the KSP and PC contexts from the KSP context,
         we can then directly call any KSP and PC routines to set
         various options.
       - The following two statements are optional; all of these
         parameters could alternatively be specified at runtime via
         KSPSetFromOptions().  All of these defaults can be
         overridden at runtime, as indicated below.
    */

    // Also on first loop, rtol can be 1e-2
    if (!*lagged && zerotolerance::on_first_lag_low_rtol)
        ierr = KSPSetTolerances(ksp, zerotolerance::low_rtol,
                                zerotolerance::KSP_atol,
                                zerotolerance::KSP_dtol,
                                zerotolerance::MatrixSolver_ItMax);
    else if (*lagged == 1 && zerotolerance::on_first_lag_low_rtol)
        ierr = KSPSetTolerances(ksp, zerotolerance::low_rtol / 10,
                                zerotolerance::KSP_atol,
                                zerotolerance::KSP_dtol,
                                zerotolerance::MatrixSolver_ItMax);
    else
        ierr = KSPSetTolerances(ksp, zerotolerance::KSP_rtol,
                                zerotolerance::KSP_atol,
                                zerotolerance::KSP_dtol,
                                zerotolerance::MatrixSolver_ItMax);
    CHKERRQ(ierr);

    // ierr = PCSetType(pc,PCASM); CHKERRQ(ierr);
    // ierr = PCSetType(pc,PCILU); CHKERRQ(ierr);

    /*
      Set runtime options, e.g.,
          -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
      These options will override those specified above as long as
      KSPSetFromOptions() is called _after_ any other customization
      routines.
    */
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc);
    CHKERRQ(ierr);
    ierr = PCFactorSetZeroPivot(pc, 1e-50);
    CHKERRQ(ierr);
    //  ierr = PCFactorSetShiftType(pc,MAT_SHIFT_POSITIVE_DEFINITE); CHKERRQ(ierr);
    //  ierr = PCFactorSetShiftAmount(pc,PETSC_DECIDE); CHKERRQ(ierr);
    //   ierr = PCFactorSetShiftPd(pc,PETSC_TRUE); CHKERRQ(ierr);
    //   ierr = PCFactorSetShiftNonzero(pc,PETSC_DECIDE); CHKERRQ(ierr);

    ierr = PetscTime(&time2);
    CHKERRQ(ierr);
    std::ostringstream Imess2;
    Imess2 << "\nPETSC: Processor " << rank << " of " << size << " - KSP solver setup took - " << time2 - time1 << " s";
    // if(!rank) Imess2<<'\n';
    if (if_time_messages)
        std::cout << Imess2.str();
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Step 7 -  Solve the linear system
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank)
        std::cout << "\n\n";
    ierr = PetscTime(&time1);
    CHKERRQ(ierr); // get time
    //  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, PETSC_vin, PETSC_vout);
    if (ierr > 0)
    {
        // this changed to help converge matrix
        if (zerotolerance::iterate_matrix == 0.0)
            zerotolerance::iterate_matrix = zerotolerance::init_matrix_it;
        else
            zerotolerance::iterate_matrix *= zerotolerance::matrix_it_multiplier;
        if (zerotolerance::iterate_matrix > 0.01)
        {
            if (!rank)
                std::cout << BRED << "IMPACTA: Matrix error - can't fix\n Exiting...\n"
                          << ENDFORMAT;
            exit(0);
        }
        else
        {
            if (rank == 0)
            {
                std::cout << BRED << "Matrix error (probably zero pivot or zero on the diagonal)\n"
                          << ENDFORMAT << "\nChanging matrix iteration tolerance to "
                          << zerotolerance::iterate_matrix
                          << " and retrying matrix solve\n";
            }
            MPI_Barrier(MPI_COMM_WORLD);
            ierr = VecDestroy(&PETSC_vout);
            CHKERRQ(ierr);
            ierr = KSPDestroy(&ksp);
            CHKERRQ(ierr);
            ierr = VecDestroy(&PETSC_vin);
            CHKERRQ(ierr);
            ierr = MatDestroy(&PETSC_S);
            CHKERRQ(ierr);
            goto startofmatrixsolve;
            //	IMPACT_Soft_Exit(M,1);
        }
    }
    CHKERRQ(ierr);
    // KSPTrueMonitor(ksp,its,converg,PETSC_NULL);
    KSPGetIterationNumber(ksp, &its);
    CHKERRQ(ierr);
    // for diagnostics....
    IMPACT_Diagnostics::total_nl_its += its;
    IMPACT_Diagnostics::nl_times += 1.0;
    //
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);

    // PetscScalar converg; //convergence info
    ierr = PetscTime(&time2);
    CHKERRQ(ierr);
    if (rank == 0)
    {
        std::ostringstream Imess3;
        Imess3 << ULINE << BYELLOW << "PETSC: Matrix Solver - " << its << " iteration(s)";
        Imess3 << " took - " << IMPACT_GetTime(time2 - time1) << std::endl
               << ENDFORMAT << ULINE;
        std::cout << Imess3.str();
    }
    if (reason > 0)
    {
        if (rank == 0)
            PetscPrintf(MPI_COMM_WORLD, "Linear solve converged due to %s\n",
                        KSPConvergedReasons[reason]);
    }
    else
    {
        // this changed to help converge matrix
        if (zerotolerance::iterate_matrix == 0.0)
            zerotolerance::iterate_matrix = zerotolerance::init_matrix_it;
        else
            zerotolerance::iterate_matrix *= zerotolerance::matrix_it_multiplier;
        if (zerotolerance::iterate_matrix > 1.0)
        {
            zerotolerance::adaptive_timesteps *= zerotolerance::adaptive_multiplier;
            badreturn = 1;
            ++zerotolerance::adaptivetimes;
        }
        else
        {
            if (rank == 0)
            {
                std::cout << BRED;
                PetscPrintf(MPI_COMM_WORLD, "Linear solve did not converge due to %s\n",
                            KSPConvergedReasons[reason]);
                std::cout << ENDFORMAT << "\nChanging matrix iteration tolerance to "
                          << zerotolerance::iterate_matrix
                          << " and retrying matrix solve\n";
            }
            MPI_Barrier(MPI_COMM_WORLD);
            ierr = VecDestroy(&PETSC_vout);
            CHKERRQ(ierr);
            ierr = KSPDestroy(&ksp);
            CHKERRQ(ierr);
            ierr = VecDestroy(&PETSC_vin);
            CHKERRQ(ierr);
            ierr = MatDestroy(&PETSC_S);
            CHKERRQ(ierr);
            goto startofmatrixsolve;
            //	IMPACT_Soft_Exit(M,1);
        }
    }
    //__________________________________________________
    if (!badreturn)
    {
        PetscScalar *vec_arr; // pointer to array
        // Step 8 - Return vector
        ierr = VecGetArray(PETSC_vout, &vec_arr);
        CHKERRQ(ierr); // This gets the elements
        for (int i = 0; i < n; ++i)
            if (fabs(vec_arr[i]) > zerotolerance::zerothreshold)
                IMPACT_vout->Set(i + 1 + Istart, vec_arr[i]); // This locally updates vout

        ierr = VecRestoreArray(PETSC_vout, &vec_arr);
        CHKERRQ(ierr);
    }
    //__________________________________________________
    //
    // Step 9 -Clean up the PETSC objects
    //
    ierr = VecDestroy(&PETSC_vout);
    CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&PETSC_vin);
    CHKERRQ(ierr);
    ierr = MatDestroy(&PETSC_S);
    CHKERRQ(ierr);
    delete[] preall;

    // This line leads to readjusting the timestep if things get easier
    if (zerotolerance::iterate_matrix < 1.0 / zerotolerance::matrix_it_multiplier && zerotolerance::adaptivetimes > 0)
        --zerotolerance::adaptivetimes;

    zerotolerance::iterate_matrix = zerotolerance::iterate_matrix_orig_val;
    return badreturn;
}
