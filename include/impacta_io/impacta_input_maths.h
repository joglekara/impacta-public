
void IMPACTA_Get_Func(std::string function, double **grid, int Nx, int Ny, double *gridx, double *gridy);

int IMPACTA_Check_Symbol(std::string symbol)
{
    int answer = 0;
    for (int i = 0; i < nsymbols; ++i)
        if (!strcmp(symbol.c_str(), IMPACTA_symbols[i].c_str()))
            answer = i + 1;
    return answer;
}
int IMPACTA_Check_Number(std::string number)
{
    int answer = 0;
    for (int i = 0; i < nnumbers; ++i)
        if (!strcmp(number.c_str(), IMPACTA_numbers[i].c_str()))
            answer = i + 1;
    return answer;
}
int IMPACTA_Check_Function(std::string function)
{
    int answer = 0;
    for (int i = 0; i < nfuncs; ++i)
        if (!strcmp(function.c_str(), IMPACTA_funcs[i].c_str()))
            answer = i + 1;
    return answer;
}
double IMPACTA_str2dbl(std::string a)
{
    std::istringstream inum("");
    std::ostringstream onum("");
    onum << a;
    inum.str(onum.str());
    double ans;
    inum >> ans;
    return ans;
}
void IMPACTA_Do_Func(double **grid, double **gridtemp, const int Nx, const int Ny, int functiontype)
{
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        {
            switch (functiontype)
            {
            case 1:
                break;
            case 2:
                grid[i][j] = cos(gridtemp[i][j]);
                break;
            case 3:
                grid[i][j] = sin(gridtemp[i][j]);
                break;
            case 4:
                grid[i][j] = tan(gridtemp[i][j]);
                break;
            case 5:
                grid[i][j] = exp(gridtemp[i][j]);
                break;
            case 6:
                grid[i][j] = tanh(gridtemp[i][j]);
                break;
            case 7:
                grid[i][j] = sinh(gridtemp[i][j]);
                break;
            case 8:
                grid[i][j] = cosh(gridtemp[i][j]);
                break;
            case 9:
                grid[i][j] = log(gridtemp[i][j]);
                break;
            }
        }
}

void IMPACTA_Get_Term(double **gridtemp, int operators, std::string v1,
                      double *gridx, double *gridy, int Nx, int Ny)
{
    int behaviour = 0;
    if (!strcmp(v1.c_str(), "X"))
        behaviour = 1;
    if (!strcmp(v1.c_str(), "Y"))
        behaviour = 2;
    if (!strcmp(v1.c_str(), "LX"))
        behaviour = 3;
    if (!strcmp(v1.c_str(), "LY"))
        behaviour = 4;
    if (!strcmp(v1.c_str(), "PI"))
        behaviour = 5;
    if (!strcmp(v1.c_str(), "FUNCTION"))
        behaviour = 6; // do nothing
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            switch (behaviour)
            {
            case 0:
                gridtemp[i][j] = IMPACTA_str2dbl(v1);
                break;
            case 1:
                gridtemp[i][j] = gridx[i];
                break;
            case 2:
                gridtemp[i][j] = gridy[j];
                break;
            case 3:
                gridtemp[i][j] = (gridx[Nx - 1] - gridx[0]) * (double)Nx / (double)(Nx - 1);
                break;
            case 4:
                gridtemp[i][j] = (gridy[Ny - 1] - gridy[0]) * (double)Ny / (double)(Ny - 1);
                break;
            case 5:
                gridtemp[i][j] = pi;
                break;
            }
}
void IMPACTA_Operate(double *grid1, double *grid2, int operators)
{
    switch (operators)
    {
    case 1:
        *grid1 += *grid2;
        break;
    case 2:
        *grid1 -= *grid2;
        break;
    case 3:
        *grid1 *= *grid2;
        break;
    case 4:
        *grid1 /= *grid2;
        break;
    case 5:
        *grid1 = pow(*grid1, *grid2);
        break;
    }
}
void IMPACTA_Do_Algebra(double **grid1, double **grid2, int operators,
                        int Nx, int Ny)
{
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            IMPACTA_Operate(&grid1[i][j], &grid2[i][j], operators);
}
void IMPACTA_Function(double **grid, const int Nx, const int Ny, int functiontype, std::string func_pow)
{
    double **gridtemp;
    gridtemp = new double *[Nx];
    for (int i = 0; i < Nx; ++i)
        gridtemp[i] = new double[Ny];
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            gridtemp[i][j] = grid[i][j];
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            IMPACTA_Do_Func(grid, gridtemp, Nx, Ny, functiontype);
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            grid[i][j] = pow(grid[i][j], IMPACTA_str2dbl(func_pow));
    delete[] gridtemp;
}
void IMPACTA_Count_Maths(std::string function, int *num)
{

    // Split functions into operators/maths
    int length = strlen(function.c_str());
    std::string check_function = "";
    int num_terms = 0;
    for (int i = 0; i < length; ++i)
    {
        check_function += function[i];
        if (IMPACTA_Check_Number(check_function))
            check_function = "";
        if (IMPACTA_Check_Symbol(check_function))
        {
            check_function = "";
            ++num_terms;
        }
    }
    *num = num_terms;
}
void IMPACTA_Do_Maths(std::string function, int *operators, std::string *values, int *num)
{

    // Split functions into operators/maths
    int length = strlen(function.c_str());
    std::string check_function = "";
    int num_terms = *num;
    for (int i = 0; i < num_terms; ++i)
    {
        operators[i] = 0;
        values[i] = "";
    }
    int stepper = 0;
    for (int i = 0; i < length; ++i)
    {
        check_function += function[i];
        if (IMPACTA_Check_Number(check_function))
        {
            while (!IMPACTA_Check_Symbol(check_function) && i < length)
            {
                values[stepper] += check_function;
                ++i;
                check_function = function[i];
            }
            --i;
            check_function = "";
            ++stepper;
        }
        if (IMPACTA_Check_Symbol(check_function))
        {
            operators[stepper] = IMPACTA_Check_Symbol(check_function);
            check_function = "";
        }
    }
}
int IMPACTA_Get_Func_n(std::string function)
{
    int length = strlen(function.c_str());
    std::string check_function = "";

    int n_sep_funcs = 1;
    std::string tempstring;
    // First check function for errors
    for (int i = 0; i < length; ++i)
    {
        check_function += function[i];
        if (IMPACTA_Check_Number(check_function))
            check_function = "";
        if (IMPACTA_Check_Symbol(check_function))
            check_function = "";
        if (IMPACTA_Check_Function(check_function))
            check_function = "";
        tempstring = "";
        tempstring += check_function;
        if (!strcmp(tempstring.c_str(), ")"))
            check_function = "";
    }

    if (strcmp(check_function.c_str(), ""))
    {
        std::cout << "IMPACTA: ERROR - In mathematical expression in input deck\n";
        std::cout << "Error: " << check_function << '\n';
        exit(0);
    }

    // then get number of functions of first order (i.e. first iteration
    // of brackets) - doesn't matter if there are too many

    int templen = 1;
    for (int i = 0; i < length; ++i)
    {

        check_function += function[i];
        if (IMPACTA_Check_Number(check_function))
            check_function = "";
        if (IMPACTA_Check_Symbol(check_function))
            check_function = "";
        if (IMPACTA_Check_Function(check_function))
        {
            templen = strlen(IMPACTA_funcs[IMPACTA_Check_Function(check_function) - 1].c_str());
            tempstring = "";
            if (i > templen + 1)
            {
                ++n_sep_funcs;
                tempstring += function[i - templen - 1];
                if (!strcmp(tempstring.c_str(), ")"))
                    --n_sep_funcs;
            }
            /**/

            check_function = "";
            int bracketcheck = 1; // one open bracket;
            while (bracketcheck > 0 && i < length)
            {
                ++i;
                check_function += function[i];
                if (!strcmp(check_function.c_str(), ")"))
                {
                    --bracketcheck;
                    check_function = "";
                }
                if (IMPACTA_Check_Number(check_function))
                    check_function = "";
                if (IMPACTA_Check_Symbol(check_function))
                    check_function = "";
                if (IMPACTA_Check_Function(check_function))
                {
                    ++bracketcheck;
                    check_function = "";
                }
            }
            check_function = "";
            tempstring = "";
            if (i < length - 1)
            {
                tempstring += function[i + 1];
                if (IMPACTA_Check_Number(tempstring) ||
                    IMPACTA_Check_Symbol(tempstring))
                    ++n_sep_funcs;
            }
        }
    }
    return n_sep_funcs;
}
void IMPACTA_Get_Func_i(std::string function, double **grid, int Nx, int Ny, double *gridx, double *gridy)
{
    int length = strlen(function.c_str());
    std::string check_function = "";

    int n_sep_funcs = IMPACTA_Get_Func_n(function);

    std::string tempstring;
    int templen = 1;

    // Make parts for building up function
    std::string *thefunctions;
    thefunctions = new std::string[n_sep_funcs];

    std::string *func_pow;
    int *func_ops;
    func_ops = new int[n_sep_funcs];         // operator on function
    func_pow = new std::string[n_sep_funcs]; // power function raised to
    int *functiontype;
    functiontype = new int[n_sep_funcs];
    // initialize values
    for (int i = 0; i < n_sep_funcs; ++i)
    {
        thefunctions[i] = "";
        func_ops[i] = 0;
        functiontype[i] = 0;
        func_pow[i] = "1.0";
    }
    int stepper = 0; // to step through the function parts
    // Now extract functions of first order

    for (int i = 0; i < length; ++i)
    {

        check_function += function[i];
        if (i == 0 && IMPACTA_Check_Symbol(check_function) != 1 &&
            IMPACTA_Check_Symbol(check_function) != 2)
            func_ops[stepper] = 1; // This to give every one an operation
        if (i == 0)
            func_ops[0] = IMPACTA_Check_Symbol(check_function);
        if (IMPACTA_Check_Number(check_function))
        {
            thefunctions[stepper] +=
                IMPACTA_numbers[IMPACTA_Check_Number(check_function) - 1];
            check_function = "";
        }
        if (IMPACTA_Check_Symbol(check_function))
        {
            if (stepper < n_sep_funcs - 1)
                func_ops[stepper + 1] = IMPACTA_Check_Symbol(check_function);
            check_function = "";
            if (i < length - 1)
            {
                tempstring = "";
                tempstring += function[i + 1];
                int test = 0;
                test += IMPACTA_Check_Number(tempstring);
                test += IMPACTA_Check_Symbol(tempstring);
                if (i < length - 2)
                    tempstring += function[i + 2];
                test += IMPACTA_Check_Number(tempstring);
                if (i < length - 3)
                    tempstring += function[i + 3];
                test += IMPACTA_Check_Number(tempstring);
                if (test)
                    thefunctions[stepper] += function[i];
            }
            else
                thefunctions[stepper] += function[i];
        }
        if (IMPACTA_Check_Function(check_function))
        {

            tempstring = "";
            templen = strlen(IMPACTA_funcs[IMPACTA_Check_Function(check_function) - 1].c_str());
            tempstring = "";
            if (i > templen + 1)
            {
                ++stepper;
                tempstring += function[i - templen - 1];
                if (!strcmp(tempstring.c_str(), ")"))
                {
                    --stepper;
                    tempstring = "";
                    tempstring += function[i - templen];
                    func_ops[stepper] = IMPACTA_Check_Symbol(tempstring);
                }
                // if (stepper<n_sep_funcs-1)
                // func_ops[stepper]=func_ops[stepper+1];}
            }

            thefunctions[stepper] = "";
            if (stepper < n_sep_funcs)
                functiontype[stepper] = IMPACTA_Check_Function(check_function);

            check_function = "";
            int bracketcheck = 1; // one open bracket;

            while (bracketcheck > 0 && i < length)
            {
                ++i;

                check_function += function[i];
                if (strcmp(check_function.c_str(), ")") || bracketcheck > 1)
                    thefunctions[stepper] += function[i];
                if (!strcmp(check_function.c_str(), ")"))
                {
                    --bracketcheck;
                    check_function = "";
                }

                if (IMPACTA_Check_Number(check_function))
                    check_function = "";
                if (IMPACTA_Check_Symbol(check_function))
                    check_function = "";
                if (IMPACTA_Check_Function(check_function))
                {
                    ++bracketcheck;
                    check_function = "";
                }
            }
            if (i < length - 2)
            {
                tempstring = function[i + 1];
                if (IMPACTA_Check_Symbol(tempstring) == 5)
                {
                    tempstring = function[i + 2];
                    func_pow[stepper] = tempstring;
                    ++i;
                    ++i;
                }
            }
            tempstring = "";
            if (i < length - 1)
            {
                tempstring += function[i + 1];
                if (IMPACTA_Check_Number(tempstring) ||
                    IMPACTA_Check_Symbol(tempstring))
                    ++stepper;
            }
            check_function = "";
        }
    }

    // Now set operation, function type and power for these functions

    // Stage 2 - now do the maths
    //
    int **operators;
    std::string **values;
    int *num_terms;
    operators = new int *[n_sep_funcs];
    values = new std::string *[n_sep_funcs];
    num_terms = new int[n_sep_funcs];
    for (int i = 0; i < n_sep_funcs; ++i)
    {
        num_terms[i] = 1;
        if (!functiontype[i])
            IMPACTA_Count_Maths(thefunctions[i], &num_terms[i]);
        operators[i] = new int[num_terms[i]];
        values[i] = new std::string[num_terms[i]];
        for (int j = 0; j < num_terms[i]; ++j)
        {
            operators[i][j] = func_ops[i];
            values[i][j] = "function";
        }
    }
    for (int i = 0; i < n_sep_funcs; ++i)
        if (!functiontype[i])
            IMPACTA_Do_Maths(thefunctions[i], operators[i],
                             values[i], &num_terms[i]);

    // Now everything is ordered.
    // now put into one long line
    int n_terms = 0;
    for (int i = 0; i < n_sep_funcs; ++i)
    {
        n_terms += num_terms[i];
    }
    int *opline;
    std::string *valline;
    int *done_job;
    opline = new int[n_terms];
    done_job = new int[n_terms];
    for (int i = 0; i < n_terms; ++i)
        done_job[i] = 1; // if term used already
    valline = new std::string[n_terms];
    stepper = 0;
    int *map_func_i;
    map_func_i = new int[n_sep_funcs];

    for (int i = 0; i < n_sep_funcs; ++i)
    {
        map_func_i[i] = 0;
        if (functiontype[i])
            map_func_i[i] = stepper;
        for (int j = 0; j < num_terms[i]; ++j)
        {
            opline[stepper] = operators[i][j];
            valline[stepper] = values[i][j];
            ++stepper;
        }
    }
    /*std::cout<<"nterms="<<n_terms<<'\n';
    for (int i=0;i<n_terms;++i)
      std::cout<<"op="<<opline[i]<<", val="<<valline[i]<<'\n';
    for (int i=0;i<n_sep_funcs;++i) std::cout<<"mapfunc="<<map_func_i[i]<<',';
    std::cout<<'\n';*/
    double ***gridtemp;
    gridtemp = new double **[n_terms];
    for (int i = 0; i < n_terms; ++i)
    {
        gridtemp[i] = new double *[Nx];
        for (int j = 0; j < Nx; ++j)
        {
            gridtemp[i][j] = new double[Ny];
            for (int k = 0; k < Ny; ++k)
                gridtemp[i][j][k] = 0.0;
        }
    }
    // Work out terms
    // Start with Functions
    //  if (n_sep_funcs>1)
    for (int i = 0; i < n_sep_funcs; ++i)
        if (functiontype[i])
            IMPACTA_Get_Func(thefunctions[i], gridtemp[map_func_i[i]],
                             Nx, Ny, gridx, gridy);

    for (int i = 0; i < n_sep_funcs; ++i)
    {
        if (functiontype[i])
            IMPACTA_Function(gridtemp[map_func_i[i]],
                             Nx, Ny, functiontype[i], func_pow[i]);
    }

    // Now evaluate terms
    for (int i = 0; i < n_terms; ++i)
        if (strcmp(valline[i].c_str(), "function"))
            IMPACTA_Get_Term(gridtemp[i], opline[i], valline[i], gridx, gridy, Nx, Ny);
    //

    /* START OF CALCULATIONS - Using Bodmas ordering*/
    // Now do powers
    for (int i = 0; i < n_terms; ++i)
        if (opline[i] == 5)
        {
            IMPACTA_Do_Algebra(gridtemp[i - 1], gridtemp[i],
                               opline[i], Nx, Ny);
            done_job[i] = 0;
        }

    // now do the divisions
    for (int i = 0; i < n_terms; ++i)
        if (opline[i] == 4 && done_job[i] == 1)
        {
            stepper = i - 1;
            while (!done_job[stepper] && stepper > 0)
                --stepper;
            IMPACTA_Do_Algebra(gridtemp[stepper], gridtemp[i],
                               opline[i], Nx, Ny);
            done_job[i] = 0;
        }

    // now do the multiplication
    for (int i = 1; i < n_terms; ++i)
        if (opline[i] == 3 && done_job[i] == 1)
        {
            stepper = i - 1;
            while (!done_job[stepper] && stepper > 0)
                --stepper;
            IMPACTA_Do_Algebra(gridtemp[stepper], gridtemp[i],
                               opline[i], Nx, Ny);
            done_job[i] = 0;
        }

    // finally sum all terms
    // first term, if - need to change
    if (opline[0] == 2)
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                gridtemp[0][i][j] = -gridtemp[0][i][j];
    // now sum all the rest onto first one
    for (int i = 1; i < n_terms; ++i)
        if ((opline[i] == 1 || opline[i] == 2) && done_job[i] == 1)
        {
            stepper = i - 1;
            while (!done_job[stepper] && stepper > 0)
                --stepper;
            IMPACTA_Do_Algebra(gridtemp[stepper], gridtemp[i],
                               opline[i], Nx, Ny);
            done_job[i] = 0;
        }
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            grid[i][j] = gridtemp[0][i][j];

    // Clean up
    delete[] thefunctions;
    delete[] func_pow;
    delete[] func_ops;
    delete[] functiontype;
    delete[] opline;

    delete[] done_job;
    delete[] map_func_i;

    for (int i = 0; i < n_sep_funcs; ++i)
    {
        delete[] operators[i];
        delete[] values[i];
    }
    delete[] operators;
    delete[] values;
    delete[] num_terms;

    for (int i = 0; i < n_terms; ++i)
    {
        for (int j = 0; j < Nx; ++j)
            delete[] gridtemp[i][j];
        delete[] gridtemp[i];
    }
    delete[] gridtemp;
}
void IMPACTA_Get_Func(std::string function, double **grid, int Nx, int Ny, double *gridx, double *gridy)
{
    // First set to uppercase
    int length = strlen(function.c_str());

    for (int i = 0; i < length; ++i)
        function[i] = (char)toupper(function.c_str()[i]);
    std::string tempfunc = "";
    tempfunc += function[0];
    if (IMPACTA_Check_Symbol(tempfunc) != 1 && IMPACTA_Check_Symbol(tempfunc) != 2)
    {
        tempfunc = "+" + function;
        function = tempfunc;
    }

    IMPACTA_Get_Func_i(function, grid, Nx, Ny, gridx, gridy);
}
