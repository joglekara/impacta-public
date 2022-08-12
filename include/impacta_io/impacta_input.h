/*
**********************************************************
IMPACTA
Function for reading input deck and setting up new
config object from it

Version 1.2
AGRT

30/3/07
11/2/08

Unbelievably badly written code - very sorry about that
**********************************************************
*/
std::string StringToUpper(std::string myString);

// Input and initial messages
IMPACT_Config IMPACT_Input(int *argnum, char **argstr)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // First open inputdeck
    /*int dircheck=0;
    for (int i=0;i<*argnum;++i)
      if (!strcmp(argstr[i], "-rd")) {dircheck=i;argstr[i]=NULL;}
    if (dircheck>0)
    {argstr[dircheck+1]; argstr[dircheck+1]=NULL;}*/

    // First sort out file management
    //*************************************************************
    IMPACT_Messages::Root_Dir = IMPACT_Get_CmdLine(argnum, argstr, "-rd");
    if (!strcmp(IMPACT_Messages::Root_Dir.c_str(), "*nofield*"))
        IMPACT_Messages::Root_Dir = "";
    else
        IMPACT_Messages::Root_Dir = IMPACT_Messages::Root_Dir + "/";
    IMPACT_Messages::Data_Directory = IMPACT_Messages::Root_Dir + "IMPACTA_Data/";
    IMPACT_Messages::Input_Directory = IMPACT_Messages::Root_Dir + "IMPACTA_Data_In/";

    if_dump_switches::ifdumpall = 1 - IMPACT_Cmp_CmdLine(*argnum, argstr, "-no_dump");
    // no xnbody wanted:
    IMPACT_Diagnostics::noxnbody = 1 - IMPACT_Cmp_CmdLine(*argnum, argstr, "-using_xnbody");

    IMPACT_Messages::InDeck = IMPACT_Get_CmdLine(argnum, argstr, "-read");
    if (!strcmp(IMPACT_Messages::InDeck.c_str(), "*nofield*"))
        IMPACT_Messages::InDeck = "imstdin";

    IMPACT_Messages::InDeck = IMPACT_Messages::Root_Dir + IMPACT_Messages::InDeck;

    std::string filename = IMPACT_Messages::InDeck;

    std::ifstream infile(filename.c_str());
    if (!infile)
        IMPACT_ferr(filename.c_str());
    //*************************************************************

    std::string input_para;    // input string from infile
    std::string old_para = ""; // for storing the previous string of the loop
    char comment_char[200];    // for reading comments
    int stringlength = 30;
    char *input_char;
    input_char = new char[stringlength + 1]; // for storing data
    int currentvar = 0;                      // for indexing what we are looking for
    int grid_err = 1;
    //*************************************************************
    // Get some command line instructions!
    int probecheck = IMPACT_Cmp_CmdLine(*argnum, argstr, "-probe_input");
    IMPACT_Diagnostics::divcheck = !IMPACT_Cmp_CmdLine(*argnum, argstr, "-no_div_check");
    if (!IMPACT_Diagnostics::divcheck && !rank)
        std::cout << "Ignoring divergence in Picard iteration\n";
    if (IMPACT_Cmp_CmdLine(*argnum, argstr, "-echo_infuncs"))
        IMPACT_Messages::if_show_function_input = 1;
    std::string truth; // for testing truth
    std::istringstream Bimptemp;
    Bimptemp.str(IMPACT_Get_CmdLine(argnum, argstr, "-B_implicit_frac"));
    if (strcmp(Bimptemp.str().c_str(), "*nofield*"))
    {
        Bimptemp >> globalconsts::Bimp;
        if (globalconsts::Bimp < 0.0 || globalconsts::Bimp > 1.0)
        {
            std::cout << "IMPACT: ERROR - Not a valid value for B implicit fraction\n";
            exit(0);
        }
    }
    // For initial DLM Distribution
    std::istringstream DLMtemp;
    DLMtemp.str(IMPACT_Get_CmdLine(argnum, argstr, "-init_DLM_order"));
    if (strcmp(DLMtemp.str().c_str(), "*nofield*"))
    {
        DLMtemp >> equation_switches::Gamma_A;
        if (equation_switches::Gamma_A < 2.0)
        {
            std::cout << "IMPACT: ERROR - Not a valid value for DLM order\n";
            std::cout << equation_switches::Gamma_A << '\n';
            exit(0);
        }
    }
    // dEbydtfix fraction
    std::istringstream Efixtemp;
    Efixtemp.str(IMPACT_Get_CmdLine(argnum, argstr, "-Fix_dE/dt"));
    if (strcmp(Efixtemp.str().c_str(), "*nofield*"))
    {
        Efixtemp >> equation_switches::dEbydtfix;
        if (equation_switches::dEbydtfix < 0.0 || equation_switches::dEbydtfix > 1.0)
        {
            std::cout << "IMPACT: ERROR - Not a valid value for dE/dt fix\n";
            exit(0);
        }
    }
    // Iterate the picard iteration
    std::istringstream PCITMATtemp;
    PCITMATtemp.str(IMPACT_Get_CmdLine(argnum, argstr, "-max_picard_its"));
    if (strcmp(PCITMATtemp.str().c_str(), "*nofield*"))
    {
        PCITMATtemp >> zerotolerance::max_picard_its;
        if (zerotolerance::iterate_matrix < 0)
        {
            std::cout << "IMPACT: ERROR - Not a valid value for max picard its\n";
            exit(0);
        }
    }
    // Change temperature normalization
    std::istringstream Tnormtemp;
    Tnormtemp.str(IMPACT_Get_CmdLine(argnum, argstr, "-renorm_temp"));
    if (strcmp(Tnormtemp.str().c_str(), "*nofield*"))
    {
        Tnormtemp >> equation_switches::NEW_T_NORM;
        if (equation_switches::NEW_T_NORM < 0)
        {
            std::cout << "IMPACT: ERROR - Not a valid value for temperature normalization\n";
            exit(0);
        }
    }
    // Iterate the matrix
    std::istringstream ITMATtemp;
    ITMATtemp.str(IMPACT_Get_CmdLine(argnum, argstr, "-iterate_matrix"));
    if (strcmp(ITMATtemp.str().c_str(), "*nofield*"))
    {
        ITMATtemp >> zerotolerance::iterate_matrix;
        if (zerotolerance::iterate_matrix < 0.0)
        {
            std::cout << "IMPACT: ERROR - Not a valid value for matrix iteration coefficient\n";
            exit(0);
        }
    }
    zerotolerance::iterate_matrix_orig_val = zerotolerance::iterate_matrix;
    if_dump_switches::view_matrix = IMPACT_Cmp_CmdLine(*argnum, argstr, "-view_matrix");
    zerotolerance::linear_solution_PC = IMPACT_Cmp_CmdLine(*argnum, argstr, "-LPC");
    zerotolerance::explicit_PC = IMPACT_Cmp_CmdLine(*argnum, argstr, "-EPC");

    zerotolerance::on_first_lag_low_rtol = IMPACT_Cmp_CmdLine(*argnum, argstr, "-KSP_lower_initial_rtol");
    equation_switches::Bimp_in_E_equation = 1.0 - IMPACT_Cmp_CmdLine(*argnum, argstr, "-B_explicit");
    if (equation_switches::Bimp_in_E_equation == 0.0)
        equation_switches::dEbydtfix = 0.0;
    if (IMPACT_Cmp_CmdLine(*argnum, argstr, "-time_centering"))
    {
        equation_switches::Bimp_in_E_equation = 0.5;
        equation_switches::jimp_in_E_equation = 0.5; // remember - can't be zero
        equation_switches::Eimp_in_B_equation = 0.5;
    }

    //*************************************************************

    // Critical information....
    int nprocs = 0; // number of processors
    int nmax = 0;
    double dt = -1;
    int Nx = 0, Ny = 0, Nv = 0;
    int NB = -1, EOn = -1, f1On = -1, f2On = -1, f3On = -1;
    for (int i = 0; i < IMPACT_Input_Deck::MAXVAR; ++i)
        IMPACT_Input_Deck::Var_check[i] = 0;
    int loopcheck = 0;
    //-1 so we can check it has been specified
    char coords[3] = {0, 0, 0};
    double *xgrid, *ygrid, *vgrid, *xtemp, *ytemp, *ttemp;
    ttemp = new double[nmax + 2];
    for (int i = 0; i < nmax + 2; ++i)
        ttemp[i] = 1.0 * i;

    // non-critical information
    double *Z_gridx, *Z_gridy, *ne_gridx, *ne_gridy, *ni_gridx, *ni_gridy;
    double *Te_gridx, *Te_gridy, *B_gridx, *B_gridy;
    double *heating_x, *heating_y, *heating_t;
    double *Dnz_x, *Dnz_y, *dnz_t, *DTz_x, *DTz_y, *dTz_t;
    double Bhat[3] = {0, 0, 0};
    double xmin = 0.0, ymin = 0.0; // where the simulation lower bound is
    double T0 = 0.0;
    int Ndump = 1;
    double wpeovernuei, vteoverc;
    int grd_val = 0;
    double gridval = 0.0;
    char temptruth;

    while (infile >> input_para && currentvar < IMPACT_Input_Deck::MAXVAR - 1)
    {
        std::ostringstream numtemp("");  // for storing data
        std::istringstream inumtemp(""); // for outputting data

        for (int k = 0; k < stringlength; ++k)
            input_char[k] = 0;            // reset input_char
        if (input_para.c_str()[0] == '%') // to deal with comments
        {
            infile.getline(comment_char, 200, '\n');
            input_para = "";
        }
        // now check if there is an = sign...
        int length = strlen(input_para.c_str());

        for (int i = 0; i < length; ++i)
            if (input_para.c_str()[i] == '=') // to deal with either spaces or not
            {                                 // between = sign and parameter
                // first make it uppercase to get rid of ambiguity

                if (!i) //(if = is the first character)
                {
                    for (int j = 0; j < (int)strlen(old_para.c_str()); ++j)
                    {
                        input_char[j] = (char)toupper(old_para.c_str()[j]);
                    }
                }
                else
                {
                    for (int j = 0; j < i; ++j)
                        input_char[j] = (char)toupper(input_para.c_str()[j]);
                }

                if (i < length - 1) //(if = is not the last character)
                {                   // the procedure is to get the string after the =
                    char tempstore; // and turn it into a number
                    for (int m = 0; m < length - i - 1; ++m)
                    {
                        tempstore = input_para.c_str()[m + i + 1];
                        if (tempstore == ',')
                            break;
                        numtemp << tempstore;
                        inumtemp.str(numtemp.str());
                    }
                }
                else
                {
                    infile >> input_para;
                    inumtemp.str(input_para);
                }
            }

        old_para = input_para;

        /*
    next comes the easy part - the istringstream object is
    emptied of its information depending on the previous
    string name - e.g. if the string before the = is
    xmin - then the data is put into xmin.
         */
        std::string grid_temp = "";
        grid_err = 1; // for extracting grid info;
        std::string dirtemp = "";
        loopcheck = 0;
        int endcheck = IMPACT_Input_Deck::MAXVAR - currentvar;
        if (endcheck > 4)
            endcheck = 4;
        for (int i = 0; i < endcheck; ++i)
            if (!strcmp(input_char, IMPACT_Input_Deck::Var_names[currentvar + i]))
            {
                loopcheck = 1;
                currentvar += i;
                break;
            }
        if (loopcheck)
        {
            switch (currentvar)
            {
            case 0:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> dirtemp;
                IMPACT_Messages::Data_Directory =
                    IMPACT_Messages::Root_Dir + dirtemp + "/";
                break;
            case 1:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> dirtemp;
                IMPACT_Messages::Input_Directory =
                    IMPACT_Messages::Root_Dir + dirtemp + "/";
                break;
            case 2:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> nprocs;
                break;
            case 3:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> nmax;
                break;
            case 4:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> dt;

                ttemp = new double[nmax + 2];
                for (int i = 0; i < nmax + 2; ++i)
                    ttemp[i] = dt * i;

                break;
            case 5:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> Ndump;
                break;
            case 6:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> T0;
                break;
            case 7:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> Nx;
                break;
            case 8:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> Ny;
                break;
            case 9:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> Nv;
                break;
            case 10:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> coords;
                break;
            case 11:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> xmin;
                break;
            case 12:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> ymin;
                break;
            case 13:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                xgrid = new double[Nx + 2];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dx:\t\t";
                grid_err = IMPACT_Make_Axis(infile, grid_temp, xgrid, Nx, xmin, 1);
                // now make xvalues...
                xtemp = new double[Nx + 2];
                xmin += 0.5 * xgrid[1];
                xtemp[1] = xmin;
                // std::cout<<"xtemp[0]: " << xtemp[0] <<" xgrid[0]: " << xgrid[0] << "\n";
                // std::cout<<"xtemp[1]: " << xtemp[1] <<" xgrid[1]: " << xgrid[1] << "\n";
                for (int i = 2; i < Nx + 2; ++i)
                    xtemp[i] = xtemp[i - 1] + 0.5 * (xgrid[i - 1] + xgrid[i]);
                // std::cout<<"xtemp[i]: " << xtemp[i] <<" xgrid[i]: " << xgrid[i] << "\n";}
                break;
            case 14:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                ygrid = new double[Ny + 2];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dy:\t\t";
                grid_err = IMPACT_Make_Axis(infile, grid_temp, ygrid, Ny, ymin, 2);
                // now make y values
                ytemp = new double[Ny + 2];
                ymin += 0.5 * ygrid[1];
                ytemp[1] = ymin;
                for (int i = 2; i < Ny + 2; ++i)
                    ytemp[i] = ytemp[i - 1] + 0.5 * (ygrid[i - 1] + ygrid[i]);
                break;
            case 15:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                vgrid = new double[Nv + 2];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dv:\t\t";
                grid_err = IMPACT_Make_Axis(infile, grid_temp, vgrid, Nv, 0.0, 3);
                vgrid[0] = vgrid[1];
                vgrid[Nv + 1] = vgrid[Nv];
                break;
            case 16:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> NB;
                break;
            case 17:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                EOn = IMPACT_truth(truth);
                break;
            case 18:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                f1On = IMPACT_truth(truth);
                break;
            case 19:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                f2On = IMPACT_truth(truth);
                break;
            case 20:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                f3On = IMPACT_truth(truth);
                break;
            case 21:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::zerothreshold;
                break;
            case 22:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::laggedtolerance;
                break;
            case 23:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::KSP_atol;
                break;
            case 24:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::KSP_rtol;
                break;
            case 25:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::KSP_dtol;
                break;
            case 26:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::MatrixSolver_ItMax;
                break;
            case 27:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                IMPACT_strcpy(argstr[*argnum - IMPACT_Input_Deck::Extra_cmdln], "-pc_type");
                IMPACT_strcpy(argstr[*argnum - IMPACT_Input_Deck::Extra_cmdln + 1], truth);
                break;
            case 28:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                IMPACT_strcpy(argstr[*argnum - IMPACT_Input_Deck::Extra_cmdln + 2], "-ksp_type");
                IMPACT_strcpy(argstr[*argnum - IMPACT_Input_Deck::Extra_cmdln + 3], truth);
                break;
            case 29:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if (IMPACT_truth(truth))
                    IMPACT_strcpy(argstr[*argnum - IMPACT_Input_Deck::Extra_cmdln + 4], "-ksp_monitor");
                break;
            case 30:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> IMPACTA_ions::atomic_mass_number;
                break;
            case 31:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> globalconsts::oneover_atomic_Z;
                globalconsts::oneover_atomic_Z = 1.0 / globalconsts::oneover_atomic_Z;

                break;
            case 32:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> IMPACTA_ions::ion_temperature;
                break;
            case 33:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp, &Initial_Conditions::Z,
                                                  xtemp, ytemp, Nx, Ny);
                Z_gridx = new double[1];
                Z_gridy = new double[1];
                break;
            case 34:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                Z_gridx = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "Z(x):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, Z_gridx, xtemp, Nx);
                break;
            case 35:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                Z_gridy = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "Z(y):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, Z_gridy, ytemp, Ny);
                // this must be parallelized at some point
                IMPACT_Make_2DGrid(&Initial_Conditions::Z, Z_gridx, Z_gridy, Nx, Ny);
                break;
            case 36:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp, &Initial_Conditions::ne,
                                                  xtemp, ytemp, Nx, Ny);
                ne_gridx = new double[1];
                ne_gridy = new double[1];
                break;
            case 37:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                ne_gridx = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "ne(x):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, ne_gridx, xtemp, Nx);
                break;
            case 38:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                ne_gridy = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "ne(y):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, ne_gridy, ytemp, Ny);
                IMPACT_Make_2DGrid(&Initial_Conditions::ne, ne_gridx, ne_gridy, Nx, Ny);
                break;
            case 39:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp, &Initial_Conditions::ni,
                                                  xtemp, ytemp, Nx, Ny);
                ni_gridx = new double[1];
                ni_gridy = new double[1];
                break;
            case 40:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                ni_gridx = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "ni(x):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, ni_gridx, xtemp, Nx);
                break;
            case 41:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                ni_gridy = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "ni(y):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, ni_gridy, ytemp, Ny);

                IMPACT_Make_2DGrid(&Initial_Conditions::ni, ni_gridx, ni_gridy, Nx, Ny);
                break;
            case 42:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                  &Initial_Conditions::Te,
                                                  xtemp, ytemp, Nx, Ny);
                Te_gridx = new double[1];
                Te_gridy = new double[1];
                break;
            case 43:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                Te_gridx = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "Te(x):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, Te_gridx, xtemp, Nx);
                break;
            case 44:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                Te_gridy = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "Te(y):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, Te_gridy, ytemp, Ny);
                IMPACT_Make_2DGrid(&Initial_Conditions::Te, Te_gridx, Te_gridy, Nx, Ny);
                break;
            case 45:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                  &Initial_Conditions::B[0],
                                                  xtemp, ytemp, Nx, Ny);
                break;
            case 46:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                  &Initial_Conditions::B[1],
                                                  xtemp, ytemp, Nx, Ny);
                break;
            case 47:
                for (int i = 0; i < 4; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                inumtemp >> grid_temp;
                currentvar += 3;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                  &Initial_Conditions::B[2],
                                                  xtemp, ytemp, Nx, Ny);
                B_gridx = new double[1];
                B_gridy = new double[1];
                break;
            case 48:
                for (int i = 0; i > -4; --i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Get_Bhat(infile, grid_temp, Bhat);
                break;

            case 49:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                B_gridx = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "B_0(x):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, B_gridx, xtemp, Nx);
                break;
            case 50:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                B_gridy = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "B_0(y):\t\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, B_gridy, ytemp, Ny);
                for (int i = 0; i < 3; ++i)
                {
                    IMPACT_Make_2DGrid(&Initial_Conditions::B[i],
                                       B_gridx, B_gridy, Nx, Ny);
                    Initial_Conditions::B[i].multiplyall(Bhat[i]);
                    // normalize to the unit vector Bhat
                }
                break;

            case 51:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> wpeovernuei;
                break;
            case 52:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> vteoverc;
                break;
            case 53:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::f0_equation_on = IMPACT_truth(truth);
                if (!strcmp(truth.c_str(), "EXMX"))
                {
                    equation_switches::evolvef0 = 1;
                    equation_switches::f0_equation_on = 0;
                }
                break;
            case 54:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::f1_equation_on = IMPACT_truth(truth);
                break;
            case 55:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::f2_equation_on = IMPACT_truth(truth);
                break;
            case 56:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::f3_equation_on = IMPACT_truth(truth);
                break;
            case 57:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::E_equation_on = IMPACT_truth(truth);
                break;
            case 58:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::B_equation_on = IMPACT_truth(truth);
                break;
            case 59:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::df0_by_dt_on = (double)IMPACT_truth(truth);
                break;
            case 60:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::inf0_grad_f1_on = (double)IMPACT_truth(truth);
                break;
            case 61:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::inf0_Edf1dv_on = (double)IMPACT_truth(truth);
                break;
            case 62:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::Cee0_on = (double)IMPACT_truth(truth);
                break;
            case 63:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::e_inert_on = (double)IMPACT_truth(truth);
                break;
            case 64:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::inf1_vgradf0_on = (double)IMPACT_truth(truth);
                break;
            case 65:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::inf1_Edf0dv_on = (double)IMPACT_truth(truth);
                break;
            case 66:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::inf1_f1xB_on = (double)IMPACT_truth(truth);
                break;
            case 67:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::Cei_on = (double)IMPACT_truth(truth);
                break;
            case 68:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::Cee1_on = IMPACT_truth(truth);
                break;
            case 69:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::e_viscosity_on = (double)IMPACT_truth(truth);
                break;
            case 70:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::df3_by_dt_on = (double)IMPACT_truth(truth);
                break;
            case 71:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::disp_j_on = (double)IMPACT_truth(truth);
                break;
            case 72:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::relaxtoeq = IMPACT_truth(truth);
                if (!strcmp(truth.c_str(), "f1only"))
                    equation_switches::relaxtoeq = 2;
                break;
            case 73:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::equil_percent;
                break;
            case 74:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::RB_D_tolerance;
                break;
            case 75:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> zerotolerance::RB_D_itmax;
                break;
            case 76:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                zerotolerance::RB_iterate = IMPACT_truth(truth);
                break;
            case 77:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if (!strcmp(truth.c_str(), "IB"))
                    IMPACT_Heating::IB_type_on = 1.0;
                if (!strcmp(truth.c_str(), "MX"))
                    IMPACT_Heating::MX_type_on = 1.0;
                break;
            case 78:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                if (IMPACT_Heating::IB_type_on == 1.0 || IMPACT_Heating::MX_type_on == 1.0)
                {
                    inumtemp >> grid_temp;
                    grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                      &IMPACT_Heating::Heating_xy, xtemp, ytemp, Nx, Ny);
                    heating_x = new double[1];
                    heating_y = new double[1];
                }
                break;
            case 79:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                if (IMPACT_Heating::IB_type_on == 1.0 || IMPACT_Heating::MX_type_on == 1.0)
                {
                    heating_x = new double[Nx];
                    inumtemp >> grid_temp;
                    if (IMPACT_Messages::if_show_function_input)
                        std::cout << "Heating(x):\t";
                    grid_err = IMPACT_Make_Grid(infile, grid_temp, heating_x, xtemp, Nx);
                }
                break;
            case 80:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                if (IMPACT_Heating::IB_type_on == 1.0 || IMPACT_Heating::MX_type_on == 1.0)
                {
                    heating_y = new double[Ny];
                    inumtemp >> grid_temp;
                    if (IMPACT_Messages::if_show_function_input)
                        std::cout << "Heating(y):\t";
                    grid_err = IMPACT_Make_Grid(infile, grid_temp, heating_y, ytemp, Ny);
                    IMPACT_Make_2DGrid(&IMPACT_Heating::Heating_xy, heating_x, heating_y, Nx, Ny);
                }
                break;
            case 81:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                heating_t = new double[nmax];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "Heating(t):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, heating_t, ttemp, nmax);
                break;
            case 82:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> IMPACT_Heating::vosc_squared;
                IMPACT_Heating::vosc_squared *= IMPACT_Heating::vosc_squared;
                break;
            case 83:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> IMPACT_Heating::polarization;
                break;
            case 84:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                  &IMPACT_Heating::Dnz_xy, xtemp, ytemp, Nx, Ny);
                Dnz_x = new double[1];
                Dnz_y = new double[1];
                break;
            case 85:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                Dnz_x = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dn/dz/n(x):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, Dnz_x, xtemp, Nx);
                break;
            case 86:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                Dnz_y = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dn/dz/n(y):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, Dnz_y, ytemp, Ny);
                IMPACT_Make_2DGrid(&IMPACT_Heating::Dnz_xy, Dnz_x, Dnz_y, Nx, Ny);
                break;
            case 87:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                dnz_t = new double[nmax];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dn/dz/n(t):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, dnz_t, ttemp, nmax);
                break;

            case 88:
                for (int i = 0; i < 3; ++i)
                    IMPACT_Input_Deck::Var_check[currentvar + i] = 1;
                currentvar += 2;
                inumtemp >> grid_temp;
                grid_err = IMPACT_Make_Function2D(infile, grid_temp,
                                                  &IMPACT_Heating::DTz_xy, xtemp, ytemp, Nx, Ny);
                DTz_x = new double[1];
                DTz_y = new double[1];
                break;
            case 89:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                IMPACT_Input_Deck::Var_check[currentvar - 1] = 1;
                DTz_x = new double[Nx];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dT/dz/T(x):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, DTz_x, xtemp, Nx);
                break;
            case 90:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                DTz_y = new double[Ny];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dT/dz/T(y):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, DTz_y, ytemp, Ny);
                IMPACT_Make_2DGrid(&IMPACT_Heating::DTz_xy, DTz_x, DTz_y, Nx, Ny);
                break;
            case 91:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                dTz_t = new double[nmax];
                inumtemp >> grid_temp;
                if (IMPACT_Messages::if_show_function_input)
                    std::cout << "dT/dz/T(t):\t";
                grid_err = IMPACT_Make_Grid(infile, grid_temp, dTz_t, ttemp, nmax);
                break;
            case 92:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_ne = IMPACT_truth(truth);
                break;
            case 93:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_ni = IMPACT_truth(truth);
                break;
            case 94:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_Ci = IMPACT_truth(truth);
                break;
            case 95:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_Z = IMPACT_truth(truth);
                break;
            case 96:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_Te = IMPACT_truth(truth);
                break;
            case 97:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_Ue = IMPACT_truth(truth);
                break;
            case 98:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_je = IMPACT_truth(truth);
                break;
            case 99:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_q = IMPACT_truth(truth);
                break;
            case 100:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_P = IMPACT_truth(truth);
                break;
            case 101:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_E = IMPACT_truth(truth);
                break;
            case 102:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_B = IMPACT_truth(truth);
                break;
            case 103:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_wt = IMPACT_truth(truth);
                break;
            case 104:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_VN = IMPACT_truth(truth);
                break;
            case 105:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_eta = IMPACT_truth(truth);
                break;
            case 106:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_f0 = IMPACT_truth(truth);
                break;
            case 107:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_f1 = IMPACT_truth(truth);
                break;
            case 108:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                if_dump_switches::dump_f2 = IMPACT_truth(truth);
                break;
            case 109:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> if_dump_switches::f_ndump;
                break;
            case 110:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> temptruth;
                switch (temptruth)
                {
                case 'p':
                    IMPACT_Boundaries::boundary_type = 0;
                    break;
                case 'r':
                    IMPACT_Boundaries::boundary_type = 1;
                    break;
                case 'b':
                    IMPACT_Boundaries::boundary_type = 2;
                    break;
                case 'R':
                    IMPACT_Boundaries::boundary_type = 4;
                    break;
                case 'o':
                    IMPACT_Boundaries::boundary_type = 3;
                    break;
                case 'P':
                    IMPACT_Boundaries::boundary_type = 5;
                    break;
                default:
                    std::cout << "\nIMPACTA: ERROR - In input deck, unknown boundaries\n";
                    exit(0);
                }
                break;
            case 111:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> temptruth;
                switch (temptruth)
                {
                case 'p':
                    IMPACT_Boundaries::boundary_type += 0;
                    break;
                case 'r':
                    IMPACT_Boundaries::boundary_type += 10;
                    break;
                case 'b':
                    IMPACT_Boundaries::boundary_type += 20;
                    break;
                case 'o':
                    IMPACT_Boundaries::boundary_type += 30;
                    break;
                default:
                    std::cout << "\nIMPACTA: ERROR - In input deck, unknown boundaries\n";
                    exit(0);
                }

                break;
            case 112:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grd_val;
                IMPACT_Boundaries::fix_f0.set(0, grd_val);
                for (int i = 1; i < 4; ++i)
                {
                    infile >> grd_val;
                    if (grd_val == 1 || grd_val == 0)
                        IMPACT_Boundaries::fix_f0.set(i, grd_val);
                    else
                    {
                        std::cout << "\nIMPACTA: ERROR - In input deck, fixed boundaries\n";
                        exit(0);
                    }
                }
                break;
            case 113:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grd_val;
                IMPACT_Boundaries::fix_f1_E.set(0, grd_val);
                for (int i = 1; i < 4; ++i)
                {
                    infile >> grd_val;
                    if (grd_val == 1 || grd_val == 0)
                        IMPACT_Boundaries::fix_f1_E.set(i, grd_val);
                    else
                    {
                        std::cout << "\nIMPACTA: ERROR - In input deck, fixed boundaries\n";
                        exit(0);
                    }
                }
                break;
            case 114:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grd_val;
                IMPACT_Boundaries::fix_f2.set(0, grd_val);
                for (int i = 1; i < 4; ++i)
                {
                    infile >> grd_val;
                    if (grd_val == 1 || grd_val == 0)
                        IMPACT_Boundaries::fix_f2.set(i, grd_val);
                    else
                    {
                        std::cout << "\nIMPACTA: ERROR - In input deck, fixed boundaries\n";
                        exit(0);
                    }
                }
                break;
            case 115:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grd_val;
                IMPACT_Boundaries::fix_B.set(0, grd_val);
                for (int i = 1; i < 4; ++i)
                {
                    infile >> grd_val;
                    if (grd_val == 1 || grd_val == 0)
                        IMPACT_Boundaries::fix_B.set(i, grd_val);
                    else
                    {
                        std::cout << "\nIMPACTA: ERROR - In input deck, fixed boundaries\n";
                        exit(0);
                    }
                }
                break;
            case 116:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> grd_val;
                IMPACT_Boundaries::fix_ni.set(0, grd_val);
                for (int i = 1; i < 4; ++i)
                {
                    infile >> grd_val;
                    if (grd_val == 1 || grd_val == 0)
                        IMPACT_Boundaries::fix_ni.set(i, grd_val);
                    else
                    {
                        std::cout << "\nIMPACTA: ERROR - In input deck, fixed boundaries\n";
                        exit(0);
                    }
                }
                break;
            case 117:
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Boundaries::fix_Ci.set(0, grd_val);
                for (int i = 1; i < 4; ++i)
                {
                    infile >> gridval;
                    IMPACT_Boundaries::fix_Ci.set(i, grd_val);
                }
                break;
            case 118:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACTA_smoothing::smooth_kernel[0] = gridval;
                double smoothtemp = IMPACTA_smoothing::smooth_kernel[0];
                // std::cout << IMPACTA_smoothing::smooth_kernel[0] << "\n";
                for (int i = 1; i < 5; ++i)
                {
                    infile >> gridval;
                    IMPACTA_smoothing::smooth_kernel[i] = gridval;
                    smoothtemp += IMPACTA_smoothing::smooth_kernel[i];
                    // std::cout << IMPACTA_smoothing::smooth_kernel[i] << "\n";  std::cout << smoothtemp << "\n";
                }
                if (!((smoothtemp - 1.0) < 1e-6))
                {
                    std::cout << "\nIMPACTA: ERROR - In input deck, smoothing kernel must sum to 1.0!\n";
                    exit(0);
                }
                break;
            }
            case 119:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                std::cout << "ionization = " << truth << '\n';

                truth = StringToUpper(truth); // make uppercase
                if (!strcmp(truth.c_str(), "NONE"))
                    IMPACTA_ions::ionization_on = 0;
                if (!strcmp(truth.c_str(), "SAHA"))
                    IMPACTA_ions::ionization_on = 1;
                if (!strcmp(truth.c_str(), "TF"))
                    IMPACTA_ions::ionization_on = 2;
                break;
            }

            case 120:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                IMPACTA_ions::quasin_on = IMPACT_truth(truth);
                break;
            }
            case 121:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                equation_switches::tr_bool = IMPACT_truth(truth);
                break;
            }
            case 122:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::beam_width = gridval;
                break;
            }
            case 123:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::beam_res = gridval;
                break;
            }

            case 124:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> truth;
                truth = StringToUpper(truth); // make uppercase
                if (!strcmp(truth.c_str(), "X"))
                    IMPACT_Heating::direction = 0;
                if (!strcmp(truth.c_str(), "Y"))
                    IMPACT_Heating::direction = 1;
                break;
            }

            case 125:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;

                inumtemp >> truth;
                truth = StringToUpper(truth); // make uppercase
                if (!strcmp(truth.c_str(), "U"))
                    IMPACT_Heating::shape = 0;
                if (!strcmp(truth.c_str(), "S"))
                    IMPACT_Heating::shape = 1;
                if (!strcmp(truth.c_str(), "G"))
                    IMPACT_Heating::shape = 2;
                break;
            }

            case 126:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::n_c = gridval;
                break;
            }

            case 127:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::ray_x0 = gridval;
                break;
            }
            case 128:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::ray_y0 = gridval;
                break;
            }
            case 129:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::theta = gridval * globalconsts::pi / 180.0;

                break;
            }
            case 130:
            {

                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;

                IMPACT_Heating::dtheta = gridval * globalconsts::pi / 180.0;
                break;
            }

            case 131:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::ds = gridval;
                break;
            }

            case 132:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                IMPACT_Heating::diffsteps = gridval;
                break;
            }

            case 133:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                Initial_Conditions::T_mult = gridval;
                break;
            }

            case 134:
            {
                IMPACT_Input_Deck::Var_check[currentvar] = 1;
                inumtemp >> gridval;
                Initial_Conditions::f0delta = gridval;
                break;
            }
            }

            inumtemp.clear();
            if (probecheck)
                std::cout << "->" << IMPACT_Input_Deck::Var_names[currentvar] << '\n';
            ++currentvar;
        }
        // Now check we haven't missed one
        // if(!strcmp(input_char,IMPACT_Input_Deck::Var_names[currentvar+1]))
        if (currentvar > 1)
            if (!IMPACT_Input_Deck::Var_check[currentvar - 1])
            {
                std::cout << "\nIMPACT: ERROR - Input deck is missing parameter or is out of order\n\n";
                std::cout << "Error at " << currentvar << " " << IMPACT_Input_Deck::Var_names[currentvar - 1] << '\n';
                exit(0);
            }
    }
    infile.close();

    /*
      finally we have to check the values are sensible
      note that don't need to check the switches EOn etc
      as an error would lead them to be off.
    */
    if (!grid_err)
        IMPACT_InputDeckError("grid error");
    if (Nx < 1 || Ny < 1 || Nv < 1)
        IMPACT_InputDeckError("Nx or Ny or Nv less than 1");
    if (NB < 0 || NB > 3 || NB == 2)
        IMPACT_InputDeckError("only 0, 1 or 3 B components allowed");
    if (dt < 0.0)
        IMPACT_InputDeckError("dt must be positive"); // can't have time running backwards!
    if (nmax < 1)
        IMPACT_InputDeckError("more than 1 timestep needed!");
    if (EOn < 0 || f1On < 0 || f2On < 0 || f3On < 0)
        IMPACT_InputDeckError("E, f1, f2 components must be 0 or 1");
    if (coords[0] != 'x' || coords[1] != 'y')
        IMPACT_InputDeckError("Coordinate error");
    // if (currentvar< IMPACT_Input_Deck::MAXVAR-1)  IMPACT_InputDeckError("");

    // check parameters all present (input deck not missing field)
    int ifanymissing = 0;
    for (int i = 0; i < IMPACT_Input_Deck::MAXVAR - 1; ++i)
    {
        if (!IMPACT_Input_Deck::Var_check[i])
        {
            std::cout << "IMPACTA: ERROR - input deck missing parameter/statement '" << IMPACT_Input_Deck::Var_names[i] << "'\n";
            ifanymissing++;
        }
    }
    if (ifanymissing)
    {
        std::cout << "input deck out of date or incorrect, please check. Exiting...\n";
        exit(0);
    }

    // check for negative temperature or density or Z
    int ne_error = 0;
    std::string ne_errormessage = "";

    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            if (Initial_Conditions::Z.get(&i, &j) < 0.0)
            {
                ne_errormessage = "Z";
                ne_error = 1;
            }
            if (Initial_Conditions::ne.get(&i, &j) < 0.0)
            {
                ne_errormessage = "ne";
                ne_error = 1;
            }
            if (Initial_Conditions::ni.get(&i, &j) < 0.0)
            {
                ne_errormessage = "ni";
                ne_error = 1;
            }
            if (Initial_Conditions::Te.get(&i, &j) < 0.0)
            {
                ne_errormessage = "Te";
                ne_error = 1;
            }
        }
    if (ne_error)
    {
        std::cout << "IMPACT: ERROR - in input deck, " << ne_errormessage << " is less than zero on at least one gridcell\n";
        exit(0);
    }
    // Now output all this information to the screen so it can be seen
    //
    double xmax = xmin + xgrid[1] * 0.5, ymax = ymin + ygrid[1] * 0.5, vmax = 0.0;
    for (int i = 1; i < Nx - 1; ++i)
        xmax += 0.5 * (xgrid[i + 1] + xgrid[i]);
    for (int i = 1; i < Ny - 1; ++i)
        ymax += 0.5 * (ygrid[i + 1] + ygrid[i]);
    for (int i = 1; i <= Nv; ++i)
        vmax += 0.5 * (vgrid[i + 1] + vgrid[i]);
    xmax += 0.5 * xgrid[Nx - 1];
    ymax += 0.5 * ygrid[Ny - 1];
    if (!rank)
        if (IMPACT_Cmp_CmdLine(*argnum, argstr, "-echo_input"))
        {
            std::cout << '\n'
                      << "From Input Deck:\nNt = " << nmax << "\nNx = " << Nx << "\nNy = " << Ny << "\nNv = " << Nv << '\n';
            std::cout << "t = (" << T0 << ":" << dt * nmax << ")\nx = (" << xmin - .5 * xgrid[1] << ":" << xmax + .5 * xgrid[Nx] << ")\ny = (" << ymin - .5 * ygrid[1] << ":" << ymax + .5 * ygrid[Ny] << ")\nv = (0:" << vmax << ")";
            std::cout << BWHITE << "\nf0: equation is "
                      << BGREEN << onoroff(equation_switches::f0_equation_on);
            std::cout << BWHITE << "\nf1: elements are " << BGREEN
                      << onoroff(f1On) << BWHITE ", equation is " << BGREEN
                      << onoroff(equation_switches::f1_equation_on);
            std::cout << BWHITE << "\nf2: elements are " << BGREEN << onoroff(f2On)
                      << BWHITE << ", equation is "
                      << BGREEN << onoroff(equation_switches::f2_equation_on);
            std::cout << BWHITE << "\nf3: elements are " << BGREEN << onoroff(f3On)
                      << BWHITE << ", equation is " << BGREEN
                      << onoroff(equation_switches::f3_equation_on);
            std::cout << BWHITE << "\nE: elements are " << BGREEN << onoroff(EOn)
                      << BWHITE << ", equation is " << BGREEN
                      << onoroff(equation_switches::E_equation_on);
            std::cout << '\n'
                      << BGREEN << NB << BWHITE
                      << " B-field components, B equation is " << BGREEN
                      << onoroff(equation_switches::B_equation_on) << '\n'
                      << ENDFORMAT;

            std::cout << "Ndump = " << Ndump << '\n';
            std::cout << "zero tolerance = " << zerotolerance::zerothreshold << '\n';
            std::cout << "lagged tolerance = " << zerotolerance::laggedtolerance << '\n';
            std::cout << "KSP atol = " << zerotolerance::KSP_atol << '\n';
            std::cout << "KSP rtol = " << zerotolerance::KSP_rtol << '\n';
            std::cout << "KSP dtol = " << zerotolerance::KSP_dtol << '\n';
            std::cout << "KSP max iterations = " << zerotolerance::MatrixSolver_ItMax << '\n';
            std::cout << "w_pe over nu_ei = " << wpeovernuei << '\n';
            std::cout << "v_te over c = " << vteoverc << '\n';
            std::cout << "delta^2 = " << globalconsts::c_L * globalconsts::c_L / globalconsts::omega_p_nuei_n / globalconsts::omega_p_nuei_n << '\n';
            std::cout << ULINE << '\n';
        }

    IMPACT_Diagnostics::output_precision = nprocs;

    // Make grids by multipling the one d functions

    // MAKE HEATING GRIDS
    IMPACT_Heating::Heating_t.ChangeSize(nmax);
    IMPACT_Heating::DTz_t.ChangeSize(nmax);
    IMPACT_Heating::Dnz_t.ChangeSize(nmax);
    for (int i = 0; i < nmax; ++i)
    {
        IMPACT_Heating::Heating_t.Set(i + 1, heating_t[i]);
        IMPACT_Heating::Dnz_t.Set(i + 1, dnz_t[i]);
        IMPACT_Heating::DTz_t.Set(i + 1, dTz_t[i]);
    }
    if (IMPACT_Heating::IB_type_on == 1.0 || IMPACT_Heating::MX_type_on == 1.0)
    {
        /* for (int i=0;i<Nx;++i)
     heating_x[i]*=Z_gridx[i]*Z_gridx[i]*ni_gridx[i]*oneover6;
     for (int i=0;i<Ny;++i)
     heating_y[i]*=Z_gridy[i]*Z_gridy[i]*ni_gridy[i]*oneover6;
     IMPACT_Make_2DGrid(&IMPACT_Heating::Heating_xy,heating_x,heating_y,Nx,Ny);


     double heating_val=0.0;
     for (int i=1;i<=Nx;++i)
     for (int j=1;j<=Ny;++j)
     {
     heating_val=Initial_Conditions::Z.get(&i,&j)*
     Initial_Conditions::Z.get(&i,&j)*
     Initial_Conditions::ni.get(&i,&j)*oneover6;
     ********************************************
     This now moved into the collision update
     ********************************************
     //IMPACT_Heating::Heating_xy.Iset(IMPACT_Heating::Heating_xy.get(&i,&j)
     //  *heating_val,&i,&j);

     }*/

        /*
    for polarization dependence, the matrix elements are the combo EE-1/3I
    E = cos theta_p xhat + sin theta_p yhat
         */
        /*if (!getcoord(IMPACT_Heating::polarization))
    for (int i=1;i<=3;++i)
      IMPACT_Heating::vosc_hat.set(i,i,0.0);
        else
    {
      if (getcoord(IMPACT_Heating::polarization)==4)
        {
          // Circularly polarized case
          IMPACT_Heating::vosc_hat.set(1,1,2.0*oneover3);
          IMPACT_Heating::vosc_hat.set(2,2,2.0*oneover3);
          IMPACT_Heating::vosc_hat.set(3,3,-oneover3);
        }
      else
        {
          if (getcoord(IMPACT_Heating::polarization)==5)
      {
        // Polarized at 45 degrees to x and y axis
        IMPACT_Heating::vosc_hat.set(1,1,0.5*oneover3);
        IMPACT_Heating::vosc_hat.set(2,2,0.5*oneover3);
        IMPACT_Heating::vosc_hat.set(1,2,0.5);
        IMPACT_Heating::vosc_hat.set(2,1,0.5);
        IMPACT_Heating::vosc_hat.set(3,3,-oneover3);
      }
          else
      {
        for (int i=1;i<=3;++i)
          IMPACT_Heating::vosc_hat.set(i,i,-oneover3);
        IMPACT_Heating::vosc_hat.set(getcoord(IMPACT_Heating::polarization),getcoord(IMPACT_Heating::polarization),2.0*oneover3);
      }
        }
        }*/
        double Elaser[3] = {0.0, 0.0, 0.0}; // unit vectors for laser
        double p = 0.0;                     // 0.0 for circular, 1.0 for linear polarization
        switch (IMPACT_Heating::polarization)
        {
        case 'x':
            Elaser[0] = 1.0;
            p = 1.0;
            break;
        case 'y':
            Elaser[1] = 1.0;
            p = 1.0;
            break;
        case 'z':
            Elaser[2] = 1.0;
            p = 1.0;
            break;
        case 'c':
            p = 0.0;
            break;
        case 'd':
            Elaser[0] = 1.0 / sqrt(2.0);
            Elaser[1] = 1.0 / sqrt(2.0);
            p = 1.0;
            break;
        default:
            break;
        }
        // FORM POLARIZATION MATRICES
        for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= 3; ++j)
                IMPACT_Heating::vosc_hat.set(i, j, Elaser[i - 1] * Elaser[j - 1]);
        IMPACT_Heating::vosc_hat.inc(1, 1, (1.0 - p));
        IMPACT_Heating::vosc_hat.inc(2, 2, (1.0 - p));
        // subtract 1/3 I
        for (int i = 1; i <= 3; ++i)
            IMPACT_Heating::vosc_hat.inc(i, i, -oneover3);
        IMPACT_Heating::Total_heating_contribution.copy(&IMPACT_Heating::Heating_xy);
        // NOW FORM B-FEEDBACK MATRICES
        //-----------------------------------------------------
        if (NB > 0)
        {
            if (NB > 1)
            {
                IMPACT_Heating::vosc_hat_Bx.set(1, 2, -0.5 * Elaser[0] * Elaser[2]);
                IMPACT_Heating::vosc_hat_Bx.set(2, 1, -0.5 * Elaser[0] * Elaser[2]);
                IMPACT_Heating::vosc_hat_Bx.set(1, 3, 0.5 * Elaser[0] * Elaser[1]);
                IMPACT_Heating::vosc_hat_Bx.set(3, 1, 0.5 * Elaser[0] * Elaser[1]);
                IMPACT_Heating::vosc_hat_Bx.set(2, 2, -Elaser[1] * Elaser[2]);
                IMPACT_Heating::vosc_hat_Bx.set(2, 3, 0.5 * ((1.0 - p) + Elaser[1] * Elaser[1] - Elaser[2] * Elaser[2]));
                IMPACT_Heating::vosc_hat_Bx.set(3, 2, 0.5 * ((1.0 - p) + Elaser[1] * Elaser[1] - Elaser[2] * Elaser[2]));
                IMPACT_Heating::vosc_hat_Bx.set(3, 3, Elaser[1] * Elaser[2]);
                //*******************************************************
                IMPACT_Heating::vosc_hat_By.set(1, 1, Elaser[0] * Elaser[2]);
                IMPACT_Heating::vosc_hat_By.set(1, 2, 0.5 * Elaser[1] * Elaser[2]);
                IMPACT_Heating::vosc_hat_By.set(2, 1, 0.5 * Elaser[1] * Elaser[2]);
                IMPACT_Heating::vosc_hat_By.set(1, 3, 0.5 * (-(1.0 - p) + Elaser[2] * Elaser[2] - Elaser[0] * Elaser[0]));
                IMPACT_Heating::vosc_hat_By.set(3, 1, 0.5 * (-(1.0 - p) + Elaser[2] * Elaser[2] - Elaser[0] * Elaser[0]));
                IMPACT_Heating::vosc_hat_By.set(2, 3, -0.5 * Elaser[0] * Elaser[1]);
                IMPACT_Heating::vosc_hat_By.set(3, 2, -0.5 * Elaser[0] * Elaser[1]);
                IMPACT_Heating::vosc_hat_By.set(3, 3, -Elaser[0] * Elaser[2]);
            }
            //*******************************************************
            IMPACT_Heating::vosc_hat_Bz.set(1, 1, -Elaser[0] * Elaser[1]);

            IMPACT_Heating::vosc_hat_Bz.set(1, 2, 0.5 * (Elaser[0] * Elaser[0] - Elaser[1] * Elaser[1]));
            IMPACT_Heating::vosc_hat_Bz.set(2, 1, 0.5 * (Elaser[0] * Elaser[0] - Elaser[1] * Elaser[1]));
            IMPACT_Heating::vosc_hat_Bz.set(1, 3, -0.5 * Elaser[1] * Elaser[2]);
            IMPACT_Heating::vosc_hat_Bz.set(3, 1, -0.5 * Elaser[1] * Elaser[2]);
            IMPACT_Heating::vosc_hat_Bz.set(2, 2, Elaser[0] * Elaser[1]);
            IMPACT_Heating::vosc_hat_Bz.set(2, 3, 0.5 * Elaser[0] * Elaser[2]);
            IMPACT_Heating::vosc_hat_Bz.set(3, 2, 0.5 * Elaser[0] * Elaser[2]);
        }
        //-----------------------------------------------------
    }
    /* IMPACT_Heating::vosc_hat.Print();
  IMPACT_Heating::vosc_hat_Bx.Print();
  IMPACT_Heating::vosc_hat_By.Print();
  IMPACT_Heating::vosc_hat_Bz.Print();
  exit(0);*/
    //-----------------------------------------------------
    if (omega_p_nuei_n < 1.0 && !equation_switches::disp_j_on)
    {
        std::cout << BRED << "IMPACT: Warning - omega_p over nu_ei is less than 1 and displacement current is switched off\n (Matrix solve error may occur)\n";
        std::cout << ENDFORMAT << "Switch on dE/dt? ";
        std::string answer;
        std::cin >> answer;
        if (!strcmp(answer.c_str(), "yes"))
            equation_switches::disp_j_on = 1.0;
    }
    /*
        Setting of global consts
        These have been altered to include temperature normalization
    */

    globalconsts::c_L = fabs(1.0 / vteoverc) / sqrt(equation_switches::NEW_T_NORM);
    globalconsts::omega_p_nuei_n = fabs(wpeovernuei) * pow(equation_switches::NEW_T_NORM, 1.5);
    IMPACT_Heating::vosc_squared /= equation_switches::NEW_T_NORM;

    //-----------------------------------------------------
    globalconsts::deltasquared = globalconsts::c_L * globalconsts::c_L / globalconsts::omega_p_nuei_n / globalconsts::omega_p_nuei_n;

    /*
  IMPACTA_ions::a1bar=me_over_mp/(2.0*oneover_atomic_Z*IMPACTA_ions::atomic_mass_number);
   IMPACTA_ions::a2bar=IMPACTA_ions::a1bar*c_L*c_L/omega_p_lambda_n/omega_p_lambda_n;
    */
    IMPACTA_ions::alpha_ion = (me_over_mp / IMPACTA_ions::atomic_mass_number / (globalconsts::oneover_atomic_Z));

    /* if (!equation_switches::disp_j_on) {equation_switches::disp_j_on+=
                 globalconsts::omega_p_lambda_n*
                 globalconsts::omega_p_lambda_n*
                 zerotolerance::dispj_off_frac;

      // if (equation_switches::disp_j_on<1e-12) equation_switches::disp_j_on=1e-12;
      }   */

    IMPACT_Config config(Nx, Ny, Nv, xgrid, ygrid, vgrid, dt, 1, f1On, f2On, f3On, EOn, NB);
    /* Cee0_flux_store::flux = new IMPACT_Vel_Sten*[Nv+1];
    for (int kk=0;kk<=Nv;++kk)
    Cee0_flux_store::flux[kk]= new IMPACT_Vel_Sten(0.0,0.0,0.0);*/

    config.set_xmin(xmin);
    config.set_ymin(ymin);
    config.set_nmax(nmax);
    config.set_ndump(Ndump);

    if (!config.Nf1())
        equation_switches::f1_equation_on = 0;
    if (!config.Nf2())
        equation_switches::f2_equation_on = 0;
    if (!config.Nf3())
        equation_switches::f3_equation_on = 0;
    if (!config.NE())
        equation_switches::E_equation_on = 0;
    if (!config.NB())
        equation_switches::B_equation_on = 0;

    /* if (!IMPACT_Cmp_CmdLine(*argnum, argstr, "-no_calc_ni"))
      {
        Initial_Conditions::ni.SwitchOffConst();
        Initial_Conditions::ni.ChangeSize(Nx,Ny);
        for (int i=1;i<=Nx;++i)
          for (int j=1;j<=Ny;++j)
      Initial_Conditions::ni.set(i,j,Initial_Conditions::ne.get(&i,&j)/
               Initial_Conditions::Z.get(&i,&j));
               }*/
    // FOR ION MOTION
    //__________________________________________________________________________

    Initial_Conditions::ni.SwitchOffConst(Nx, Ny);

    IMPACTA_ions::ion_motion = !IMPACT_Cmp_CmdLine(*argnum, argstr, "-static_ions");

    if (IMPACTA_ions::ion_motion)
    {

        for (int dir = 0; dir < 3; ++dir)
        {
            Initial_Conditions::C_i[dir].SwitchOffConst();
            Initial_Conditions::C_i[dir].ChangeSize(Nx, Ny);
            for (int i = 1; i <= Nx; ++i)
                for (int j = 1; j <= Ny; ++j)
                    Initial_Conditions::C_i[dir].set(i, j, 0.0);
            Initial_Conditions::C_istar[dir].SwitchOffConst();
            Initial_Conditions::C_istar[dir].ChangeSize(Nx, Ny);
            for (int i = 1; i <= Nx; ++i)
                for (int j = 1; j <= Ny; ++j)
                    Initial_Conditions::C_istar[dir].set(i, j, 0.0);
        }
        Initial_Conditions::DivC_i.SwitchOffConst();
        Initial_Conditions::DivC_i.ChangeSize(Nx, Ny);
        for (int i = 1; i <= Nx; ++i)
            for (int j = 1; j <= Ny; ++j)
                Initial_Conditions::DivC_i.set(i, j, 0.0);
    }
    else
    {
        for (int dir = 0; dir < 3; ++dir)
        {
            Initial_Conditions::C_i[dir].setc(0.0);
            Initial_Conditions::C_istar[dir].setc(0.0);
        }
        Initial_Conditions::DivC_i.setc(0.0);
    }
    //__________________________________________________________________________

    // For ionization
    Initial_Conditions::Z.SwitchOffConst(Nx, Ny);

    // C_istar is used to smooth Z
    if (IMPACTA_ions::ionization_on && !IMPACTA_ions::ion_motion)
    {
        for (int dir = 0; dir < 3; ++dir)
        {
            Initial_Conditions::C_istar[dir].SwitchOffConst();
            Initial_Conditions::C_istar[dir].ChangeSize(Nx, Ny);
            for (int i = 1; i <= Nx; ++i)
                for (int j = 1; j <= Ny; ++j)
                    Initial_Conditions::C_istar[dir].set(i, j, 0.0);
        }
    }

    //__________________________________________________________________________

    if (IMPACT_Cmp_CmdLine(*argnum, argstr, "-show_all_initial"))
    {
        std::cout << BPURPLE << "Z\n"
                  << YELLOW;
        Initial_Conditions::Z.Print();
        std::cout << BPURPLE << "ni\n"
                  << YELLOW;
        Initial_Conditions::ni.Print();
        std::cout << BPURPLE << "ne\n"
                  << YELLOW;
        Initial_Conditions::ne.Print();
        std::cout << BPURPLE << "Te\n"
                  << YELLOW;
        Initial_Conditions::Te.Print();
        for (int dim = 0; dim < 3; ++dim)
        {
            std::cout << BPURPLE << "B_x" << dim + 1 << '\n'
                      << YELLOW;
            Initial_Conditions::B[dim].Print();
        }
        std::cout << ENDFORMAT;
    }
    // For dbydz grids

    IMPACT_Heating::DTbydzgrid.SwitchOffConst(Nx, Ny);
    IMPACT_Heating::Dnbydzgrid.SwitchOffConst(Nx, Ny);

    Initial_Conditions::ne.SwitchOffConst(Nx, Ny);
    Initial_Conditions::Te.SwitchOffConst(Nx, Ny);

    double gridtemp = 0.0;
    for (int i = 1; i <= Nx; ++i)
        for (int j = 1; j <= Ny; ++j)
        {
            gridtemp = IMPACT_Heating::Dnz_xy.get(&i, &j) * IMPACT_Heating::Dnz_t.Get(1) / Initial_Conditions::ne.get(&i, &j);

            IMPACT_Heating::Dnbydzgrid.Iset(gridtemp, &i, &j);

            gridtemp = IMPACT_Heating::DTz_xy.get(&i, &j) * IMPACT_Heating::DTz_t.Get(1) / Initial_Conditions::Te.get(&i, &j);

            IMPACT_Heating::DTbydzgrid.Iset(gridtemp, &i, &j);
        }

    for (int i = 0; i < 4; ++i)
    {
        IMPACT_Boundaries::fixed_any.inc(i, IMPACT_Boundaries::fix_f0[i]);
        IMPACT_Boundaries::fixed_any.inc(i, IMPACT_Boundaries::fix_f1_E[i]);
        IMPACT_Boundaries::fixed_any.inc(i, IMPACT_Boundaries::fix_f2[i]);
        IMPACT_Boundaries::fixed_any.inc(i, IMPACT_Boundaries::fix_B[i]);
        IMPACT_Boundaries::fixed_any.inc(i, IMPACT_Boundaries::fix_ni[i]);
        IMPACT_Boundaries::fixed_any.inc(i, (IMPACT_Boundaries::fix_Ci[i] != 0.0));

        if (IMPACT_Boundaries::fixed_any[i])
        {
            IMPACT_Boundaries::fixed_any.set(i, 1);
        }
        else
        {
            IMPACT_Boundaries::fixed_any.set(i, 0);
        }
    }

    // set smoothing kernel
    for (int i = 0; i < 5; ++i)
    {
        LAX.set(i, IMPACTA_smoothing::smooth_kernel[i]);
    }

    // for diagnosis
    IMPACT_Diagnostics::showmatrix =
        IMPACT_Cmp_CmdLine(*argnum, argstr, "-show_matrix");
    IMPACT_Diagnostics::showvector =
        IMPACT_Cmp_CmdLine(*argnum, argstr, "-show_vector");
    IMPACT_strcpy(argstr[*argnum - IMPACT_Input_Deck::Extra_cmdln + 5], "-mat_view_info");

    delete[] xgrid;
    delete[] ygrid;
    delete[] vgrid;
    delete[] xtemp;
    delete[] ytemp;
    delete[] ttemp;
    delete[] Z_gridx;
    delete[] Z_gridy;
    delete[] ne_gridx;
    delete[] ne_gridy;
    delete[] ni_gridx;
    delete[] ni_gridy;
    delete[] Te_gridx;
    delete[] Te_gridy;
    delete[] B_gridx;
    delete[] B_gridy;
    if (IMPACT_Heating::IB_type_on == 1.0 || IMPACT_Heating::MX_type_on == 1.0)
    {
        delete[] heating_x;
        delete[] heating_y;
        delete[] heating_t;
    }
    delete[] input_char;
    return config;
}

// Routine for uppercasing string
std::string StringToUpper(std::string myString)
{
    const int length = myString.length();
    for (int i = 0; i != length; ++i)
    {
        myString[i] = std::toupper(myString[i]);
    }
    return myString;
}
