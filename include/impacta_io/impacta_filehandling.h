/*
**********************************************************
File handling - functions for writing and reading moments

Version 2.6
AGRT

27/3/07

25/7/07 - Finally updated so that the actual distribution functions
          are output.

24/9/07 - Now copies the input deck into the data directory as a
record of what the parameters are
7/4/08 - Added a dump for wt i.e. sqrt(B^2)*T^3/2/Z/ne
**********************************************************
*/

// To set up directory structure
int IMPACT_Tree(IMPACT_Config *c)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = 0;
    if (!rank)
    {
        std::string unix_com = "test -d " + IMPACT_Messages::Data_Directory;
        result = system(unix_com.c_str());
        if (!result)
        {
            if (!IMPACTA_Version_HPC)
            {
                std::cout << "IMPACT: ERROR - Directory structure '" << IMPACT_Messages::Data_Directory << "' already exists\n";
                std::cout << "\nOverwrite? (type 'yes' for yes) ";
                std::string answer;
                std::cin >> answer;
                if (!strcmp(answer.c_str(), "yes"))
                {
                    std::string root = "rm -r " + IMPACT_Messages::Data_Directory;
                    result = system(root.c_str());
                    if (result != 0)
                    {
                        std::cout << "\nCould not remove directory";
                        IMPACT_Bad_Exit(0);
                    }
                }
                else
                {
                    IMPACT_Bad_Exit(0);
                }
            }
            else
            {
                std::cout << "IMPACT: ERROR - Directory structure '" << IMPACT_Messages::Data_Directory << "' already exists\n";
                IMPACT_Bad_Exit(0);
            }
        }
        std::string root = "mkdir " + IMPACT_Messages::Data_Directory;
        result = system(root.c_str());
        unix_com = root + IMPACT_Messages::Field_Dir;
        result = system(unix_com.c_str());
        unix_com = root + IMPACT_Messages::DistFunc_Dir;
        result = system(unix_com.c_str());
        unix_com = root + IMPACT_Messages::Moment_Dir;
        result = system(unix_com.c_str());
        unix_com = root + IMPACT_Messages::Constants_Dir;
        result = system(unix_com.c_str());

        // Copy input deck to directory
        time_t rawtime;
        struct tm *timeinfo;

        time(&rawtime);
        timeinfo = localtime(&rawtime);
        std::string dateStr = asctime(timeinfo);

        unix_com = IMPACT_Messages::Data_Directory + "imstdin--" + dateStr;
        for (int i = 0; i < (int)unix_com.length(); ++i)
        {
            if (unix_com[i] == ' ')
                unix_com.replace(i, 1, "_");
        }
        unix_com = "cp " + IMPACT_Messages::InDeck + " " + unix_com;
        result = system(unix_com.c_str());

        // now put in appropriate directories for the Fields
        if (c->NE() > 0 && if_dump_switches::dump_E)
            for (int x1 = 1; x1 <= c->NE(); ++x1)
            {
                unix_com = root + IMPACT_Messages::Field_Dir + "E" + getcoord(x1);
                result = system(unix_com.c_str());
            }
        if (c->NB() > 0 && if_dump_switches::dump_B)
            for (int x1 = 3; x1 > 3 - c->NB(); --x1)
            {
                unix_com = root + IMPACT_Messages::Field_Dir + "B" + getcoord(x1);
                result = system(unix_com.c_str());
            }

        // Now the directories for the ions
        if (if_dump_switches::dump_ni && IMPACTA_ions::ion_motion)
        {
            unix_com = root + IMPACT_Messages::Constants_Dir + "ni";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_Ci && IMPACTA_ions::ion_motion)
            for (int x1 = 1; x1 <= c->Nf1(); ++x1)
            {
                unix_com = root + IMPACT_Messages::Constants_Dir + "C" + getcoord(x1);
                result = system(unix_com.c_str());
            }
        if (if_dump_switches::dump_Z && (IMPACTA_ions::ion_motion || IMPACTA_ions::ionization_on))
        {
            unix_com = root + IMPACT_Messages::Constants_Dir + "Z";
            result = system(unix_com.c_str());
        }
        // Now the dirs for the moments
        if (if_dump_switches::dump_ne)
        {
            unix_com = root + IMPACT_Messages::Moment_Dir + "ne";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_Ue)
        {
            unix_com = root + IMPACT_Messages::Moment_Dir + "Ue";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_Te)
        {
            unix_com = root + IMPACT_Messages::Moment_Dir + "Te";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_eta)
        {
            unix_com = root + IMPACT_Messages::Moment_Dir + "eta";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_wt)
        {
            unix_com = root + IMPACT_Messages::Moment_Dir + "wt";
            result = system(unix_com.c_str());
        }
        if (c->Nf1() > 0 && if_dump_switches::dump_je)
            for (int x1 = 1; x1 <= c->Nf1(); ++x1)
            {
                unix_com = root + IMPACT_Messages::Moment_Dir + "j" + getcoord(x1);
                result = system(unix_com.c_str());
            }
        if (c->Nf1() > 0 && if_dump_switches::dump_q)
            for (int x1 = 1; x1 <= c->Nf1(); ++x1)
            {
                unix_com = root + IMPACT_Messages::Moment_Dir + "qT" + getcoord(x1);
                result = system(unix_com.c_str());
                unix_com = root + IMPACT_Messages::Moment_Dir + "qe" + getcoord(x1);
                result = system(unix_com.c_str());
            }
        if (c->Nf1() > 0 && if_dump_switches::dump_VN)
            for (int x1 = 1; x1 <= c->Nf1(); ++x1)
            {
                unix_com = root + IMPACT_Messages::Moment_Dir + "VN" + getcoord(x1);
                result = system(unix_com.c_str());
            }
        if (c->Nf2() > 0 && if_dump_switches::dump_P)
            for (int x1 = 1; x1 <= c->N3f2(); ++x1)
                for (int x2 = x1; x2 <= c->N3f2(); ++x2)
                    if (x1 + x2 < 6)
                    {
                        unix_com = root + IMPACT_Messages::Moment_Dir + "P" + getcoord(x1) + getcoord(x2);
                        result = system(unix_com.c_str());
                    }
        // Finally the dirs for the distribution functions
        if (if_dump_switches::dump_f0)
        {
            unix_com = root + IMPACT_Messages::DistFunc_Dir + "f0";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_f1)
        {
            unix_com = root + IMPACT_Messages::DistFunc_Dir + "f1";
            result = system(unix_com.c_str());
        }
        if (if_dump_switches::dump_f2)
        {
            unix_com = root + IMPACT_Messages::DistFunc_Dir + "f2";
            result = system(unix_com.c_str());
        }
        /*   if(if_dump_switches::dump_f3)
      {
        unix_com=root+IMPACT_Messages::DistFunc_Dir+"f3";
        result = system(unix_com.c_str());
        }*/
    }
    return result;
}
void IMPACT_ferr(const char *file)
{
    std::cout << "IMPACT: ERROR - file " << file << " not found" << std::endl;
    chkMPI();
}
std::string IMPACT_Out(std::string file)
{
    // Output data to IMPACT_Messages::Data_Directory
    std::string name;
    std::string dir = IMPACT_Messages::Data_Directory;
    name = dir + file;
    /*std::ofstream fileout;
    fileout.open(name.c_str());
    if (!fileout) IMPACT_ferr(name.c_str());*/
    return name.c_str();
}
std::string IMPACT_In(std::string file)
{
    // Open files from IMPACT_Messages::Input_Directory
    std::string name;
    std::string dir = IMPACT_Messages::Input_Directory;
    name = dir + file;
    return name.c_str();
}
std::string IMPACT_Out(std::string dir, std::string file, int index)
{
    // Open files from IMPACT_Messages::Input_Directory
    std::ostringstream numtemp;
    numtemp << index;
    std::string number;
    if (index < 10)
        number = '0';
    if (index < 100)
        number += '0';
    number += numtemp.str();
    std::string name = dir + file + number + ".imd";
    return name.c_str();
}
// Output data to IMPACT_Messages::Data_Directory
void IMPACT_Write(struct IMPACT_Moment *M, std::string str)
{
    std::ofstream outfile;
    std::string name = IMPACT_Out(str);
    outfile.open(name.c_str());
    if (!outfile)
        IMPACT_ferr(name.c_str());
    outfile << M;
    outfile.close();
}
// This format automatically outputs a .imd - so no need to add extension
void IMPACT_Write(struct IMPACT_Moment *M, std::string str, int index)
{
    std::ofstream outfile;
    std::ostringstream numtemp;
    numtemp << index;
    std::string number;
    if (index < 10)
        number = '0';
    if (index < 100)
        number += '0';
    number += numtemp.str();
    std::string name = IMPACT_Out(str) + number + ".imd";
    outfile.open(name.c_str());
    if (!outfile)
        IMPACT_ferr(name.c_str());
    outfile << M;
    outfile.close();
}
// Open files from IMPACT_Messages::Input_Directory
void IMPACT_Read(struct IMPACT_Moment *M, std::string str)
{
    std::ifstream infile;
    std::string name = IMPACT_In(str);
    infile.open(name.c_str());
    if (!infile)
        IMPACT_ferr(name.c_str());
    infile >> M;
    infile.close();
}
void IMPACT_Write(struct IMPACT_Moment *M, std::string str, std::string dir)
{
    std::ofstream outfile;
    std::string name = dir + str;
    outfile.open(name.c_str());
    if (!outfile)
        IMPACT_ferr(name.c_str());
    outfile << M;
    outfile.close();
}
// Open files from IMPACT_Messages::Input_Directory
void IMPACT_Read(struct IMPACT_Moment *M, std::string str, std::string dir)
{
    std::ifstream infile;
    std::string name = dir + str;
    infile.open(name.c_str());
    if (!infile)
        IMPACT_ferr(name.c_str());
    infile >> M;
    infile.close();
}
void IMPACT_Write(struct IMPACT_Moment *M, std::string str, std::string dir, int timestep)
{
    std::ofstream outfile;
    std::string name = IMPACT_Out(dir, str, timestep);
    outfile.open(name.c_str());
    if (!outfile)
        IMPACT_ferr(name.c_str());
    outfile << M;
    outfile.close();
}
void IMPACT_WriteInt(std::string str, std::string dir, int timestep)
{
    std::ofstream outfile;
    std::string name = IMPACT_Out(dir, str, timestep);
    outfile.open(name.c_str());
    if (!outfile)
        IMPACT_ferr(name.c_str());
    // outfile<<M;
    IMPACT_Heating::i_mat.Print(outfile);
    outfile.close();
}
void IMPACT_Write(IMPACT_Dist *M, std::string str, std::string dir, int timestep)
{
    std::ofstream outfile;
    std::string name = IMPACT_Out(dir, str, timestep);
    outfile.open(name.c_str());
    if (!outfile)
        IMPACT_ferr(name.c_str());
    outfile << M;
    outfile.close();
}
// To dump all the data at a timestep we use....
void IMPACT_Dump(IMPACT_Config *c, IMPACT_MPI_Config *MPIc, IMPACT_ParVec *v, IMPACT_Var *f0, IMPACT_Var *f1, IMPACT_Var *f2, IMPACT_Var *f3, IMPACT_Var *E, IMPACT_Var *B, int timestep)
{

    // If new T normalization invoked, need to change all vals
    IMPACTA_Rescale_T(c, MPIc, v, f0, f1, f2, f3, E, B, 1.0 / equation_switches::NEW_T_NORM);

    // easy way of checking whether this timestep should be dumped...

    const int ifdumpfields = (timestep % c->N_Dump());
    const int ifdumpdists = (timestep % (if_dump_switches::f_ndump * c->N_Dump()));
    std::ostringstream Imessage;

    if ((!ifdumpfields || timestep == c->n_max()) && if_dump_switches::ifdumpall == 1)
    {
        clock_t timestart, timeend; // for timing
        timestart = clock();
        if (!MPIc->rank())
        {
            Imessage << "\n"
                     << ENDFORMAT << "IMPACT: Writing data\n|<";
            int total = (2 + 2 * c->Nf1() + c->Nf2() + c->NE() + c->NB()) - 2;

            for (int i = 0; i < total; ++i)
                Imessage << "--";
            Imessage << ">|\n";
            std::cout << Imessage.str();
        }

        std::string dir;
        std::string name;
        struct IMPACT_Moment ZM;

        new_Constant(&Initial_Conditions::Z, c, MPIc, &ZM, "Z");

        if (!timestep && !MPIc->rank())
        {

            struct IMPACT_Moment niM;
            new_Constant(&Initial_Conditions::ni, c, MPIc, &niM, "ni");

            dir = IMPACT_Messages::Data_Directory + IMPACT_Messages::Constants_Dir;
            if (if_dump_switches::dump_Z)
                IMPACT_Write(&ZM, "Z_", dir, timestep);
            if (if_dump_switches::dump_ni && !IMPACTA_ions::ion_motion)
                IMPACT_Write(&niM, "ni_", dir, timestep);
        }
        if (IMPACTA_ions::ion_motion && !MPIc->rank())
        {
            struct IMPACT_Moment CiM;

            for (int direc = 0; direc < c->Nf1(); ++direc)
            {
                name = "C";
                name += getcoord(direc + 1);
                new_Constant(&Initial_Conditions::C_i[direc], c, MPIc, &CiM, name);
                dir = IMPACT_Messages::Data_Directory +
                      IMPACT_Messages::Constants_Dir + name + "/";
                name += "_";
                if (if_dump_switches::dump_Ci)
                    IMPACT_Write(&CiM, name, dir, timestep);
            }
            struct IMPACT_Moment niM;
            new_Constant(&Initial_Conditions::ni, c, MPIc, &niM, "ni");
            dir = IMPACT_Messages::Data_Directory +
                  IMPACT_Messages::Constants_Dir + "ni/";

            if (if_dump_switches::dump_ni)
                IMPACT_Write(&niM, "ni_", dir, timestep);

            struct IMPACT_Moment ZM;
            new_Constant(&Initial_Conditions::Z, c, MPIc, &ZM, "Z");
            dir = IMPACT_Messages::Data_Directory +
                  IMPACT_Messages::Constants_Dir + "Z/";
            if (if_dump_switches::dump_Z)
                IMPACT_Write(&ZM, "Z_", dir, timestep);
        }

        if (IMPACTA_ions::ionization_on && !MPIc->rank() && !IMPACTA_ions::ion_motion)
        {
            struct IMPACT_Moment ZM;
            new_Constant(&Initial_Conditions::Z, c, MPIc, &ZM, "Z");
            dir = IMPACT_Messages::Data_Directory +
                  IMPACT_Messages::Constants_Dir + "Z/";
            if (if_dump_switches::dump_Z)
                IMPACT_Write(&ZM, "Z_", dir, timestep);
        }
        //________________________________________________________________
        // First output Moments
        MPI_Barrier(MPI_COMM_WORLD);

        struct IMPACT_Moment ne;

        MPI_Barrier(MPI_COMM_WORLD);

        new_ne(v, f0, c, MPIc, &ne);

        struct IMPACT_Moment Ue;
        MPI_Barrier(MPI_COMM_WORLD);

        new_Ue(v, f0, c, MPIc, &Ue);

        struct IMPACT_Moment Te = Divide(&Ue, &ne);
        for (int i = 1; i <= c->Nx(); ++i)
            for (int j = 1; j <= c->Ny(); ++j)
                Te.values.set(i, j, Te.values.get(i, j) * 2.0 / 3.0);

        struct IMPACT_Moment eta;
        MPI_Barrier(MPI_COMM_WORLD);
        new_eta(v, f0, c, MPIc, &eta);

        MPI_Barrier(MPI_COMM_WORLD);

        struct IMPACT_Moment m05;

        MPI_Barrier(MPI_COMM_WORLD);

        new_m0five(v, f0, c, MPIc, &m05);
        MPI_Barrier(MPI_COMM_WORLD);

        struct IMPACT_Moment m03;

        MPI_Barrier(MPI_COMM_WORLD);

        new_m0three(v, f0, c, MPIc, &m03);

        if (!MPIc->rank())
        {
            dir = IMPACT_Messages::Data_Directory + IMPACT_Messages::Moment_Dir;
            if (if_dump_switches::dump_ne)
                IMPACT_Write(&ne, "ne_", dir + "ne/", timestep);
            std::cout << BRED << "[]";
            if (if_dump_switches::dump_Ue)
                IMPACT_Write(&Ue, "Ue_", dir + "Ue/", timestep);
            std::cout << BRED << "[]";
            IMPACT_Write(&m05, "m05_", dir + "Ue/", timestep);
            IMPACT_Write(&m03, "m03_", dir + "Ue/", timestep);

            if (if_dump_switches::dump_Te)
            {
                strcpy(Te.name, "T_e");
                IMPACT_Write(&Te, "Te_", dir + "Te/", timestep);

                // Ray Tracing Output
                IMPACT_WriteInt("Int_", dir + "Te/", timestep);
            }

            if (if_dump_switches::dump_eta)
            {
                strcpy(eta.name, "eta");
                IMPACT_Write(&eta, "eta_", dir + "eta/", timestep);
            }
        }

        struct IMPACT_Moment je[3], q[3], VN[3];

        if (c->Nf1() > 0)
            for (IMPACT_Dim x1 = 1; x1 <= c->Nf1(); ++x1)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                new_je(v, f1, c, MPIc, &je[x1.get() - 1], &x1);

                name = "j";
                name += getcoord(x1.get());

                if (!MPIc->rank())
                {
                    if (if_dump_switches::dump_je)
                        IMPACT_Write(&je[x1.get() - 1], name + "_", dir + name + "/", timestep);
                    std::cout << BRED << "[]";
                }
                MPI_Barrier(MPI_COMM_WORLD);

                // Total heat flow
                new_qT(v, f1, c, MPIc, &q[x1.get() - 1], &x1);
                name = "qT";
                name += getcoord(x1.get());

                if (!MPIc->rank())
                {
                    if (if_dump_switches::dump_q)
                        IMPACT_Write(&q[x1.get() - 1], name + "_", dir + name + "/", timestep);
                    std::cout << BRED << "[]";
                }

                // Now dump intrinsic heat flow
                for (int i = 1; i <= c->Nx(); ++i)
                    for (int j = 1; j <= c->Ny(); ++j)
                        q[x1.get() - 1].values.set(i, j, q[x1.get() - 1].values.get(i, j) + 1.25 * je[x1.get() - 1].values.get(i, j));
                // end of loop!
                name = "qe";
                name += getcoord(x1.get());
                if (!MPIc->rank())
                {
                    if (if_dump_switches::dump_q)
                        IMPACT_Write(&q[x1.get() - 1], name + "_", dir + name + "/", timestep);
                    std::cout << BRED << "[]";
                }

                MPI_Barrier(MPI_COMM_WORLD);
                // ve is electron flow velocity
                struct IMPACT_Moment ve = Divide(&je[x1.get() - 1], &ne);

                MPI_Barrier(MPI_COMM_WORLD);

                new_VN(v, f0, f1, c, MPIc, &VN[x1.get() - 1], &x1);

                for (int i = 1; i <= c->Nx(); ++i)
                    for (int j = 1; j <= c->Ny(); ++j)
                    {

                        VN[x1.get() - 1].values.set(i, j, VN[x1.get() - 1].values.get(i, j) + ve.values.get(i, j));
                    }
                // end of loop!
                name = "VN";
                name += getcoord(x1.get());

                if (!MPIc->rank())
                {
                    if (if_dump_switches::dump_VN)
                        IMPACT_Write(&VN[x1.get() - 1], name + "_", dir + name + "/", timestep);
                    std::cout << BRED << "[]";
                }
            }

        struct IMPACT_Moment P[5];
        for (IMPACT_Dim x1 = 1; x1 <= c->N3f2(); ++x1)
            for (IMPACT_Dim x2 = x1.get(); x2 <= c->N3f2(); ++x2)
                if (x1.get() + x2.get() < 6)
                {
                    MPI_Barrier(MPI_COMM_WORLD);
                    new_P(v, f2, c, MPIc, &P[x1.get() + x2.get() - 2], &x1, &x2);
                    name = "P";
                    name += getcoord(x1.get());
                    name += getcoord(x2.get());
                    if (!MPIc->rank())
                    {
                        if (if_dump_switches::dump_P)
                            IMPACT_Write(&P[x1.get() + x2.get() - 2], name + "_", dir + name + "/", timestep);
                        std::cout << BRED << "[]";
                    }
                    //________________________________________________________________
                    // NEW bit AGRT 2012 - eventually include in input deck but for now....
                    // This outputs nonlocal version of divP term in Ohm's law
                    std::string name2;
                    MPI_Barrier(MPI_COMM_WORLD);
                    new_divP(v, f0, f2, c, MPIc, &P[x1.get() + x2.get() - 2], &x1, &x2);
                    name2 = "divP";
                    name2 += getcoord(x1.get());
                    name2 += getcoord(x2.get());
                    if (!MPIc->rank())
                    {
                        if (if_dump_switches::dump_P)
                            IMPACT_Write(&P[x1.get() + x2.get() - 2], name2 + "_", dir + name + "/", timestep);
                        std::cout << BRED << "[]";
                    }
                    //________________________________________________________________
                }
        //________________________________________________________________
        // Now output fields
        dir = IMPACT_Messages::Data_Directory + IMPACT_Messages::Field_Dir;
        struct IMPACT_Moment EM[3], BM[3];
        if (c->NE() > 0)
            for (IMPACT_Dim x1 = 1; x1 <= c->NE(); ++x1)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                new_E(v, E, c, MPIc, &EM[x1.get() - 1], &x1);
                name = "E";
                name += getcoord(x1.get());
                if (!MPIc->rank())
                {
                    if (if_dump_switches::dump_E)
                        IMPACT_Write(&EM[x1.get() - 1], name + "_", dir + name + "/", timestep);
                    if (!MPIc->rank())
                        std::cout << BRED << "[]";
                }
            }
        struct IMPACT_Moment wt;
        double wt_temp = 0.0;
        if (c->NB() > 0)
        {
            for (IMPACT_Dim x1 = 3; x1 > 3 - c->NB(); --x1)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                new_B(v, B, c, MPIc, &BM[x1.get() - 1], &x1);
                name = "B";
                name += getcoord(x1.get());
                if (!MPIc->rank())
                {
                    if (if_dump_switches::dump_B)
                        IMPACT_Write(&BM[x1.get() - 1], name + "_", dir + name + "/", timestep);
                    if (!MPIc->rank())
                        std::cout << BRED << "[]";
                }
            }
            if (if_dump_switches::dump_wt)
            {
                name = "wt";
                IMPACT_Dim x1_temp(1);
                new_B(v, B, c, MPIc, &wt, &x1_temp);
                strcpy(wt.name, name.c_str());

                dir = IMPACT_Messages::Data_Directory + IMPACT_Messages::Moment_Dir;
                for (int i = 1; i <= c->Nx(); ++i)
                    for (int j = 1; j <= c->Ny(); ++j)
                    {
                        wt_temp = 0.0;
                        for (IMPACT_Dim x1 = 3; x1 > 3 - c->NB(); --x1)
                            wt_temp += BM[x1.get() - 1].values.get(i, j) * BM[x1.get() - 1].values.get(i, j);
                        wt_temp = sqrt(wt_temp) / ne.values.get(i, j) / ZM.values.get(i, j);
                        wt_temp *= pow(Te.values.get(i, j), 1.5) * 3.0 / 4.0 * sqrt(globalconsts::pi);
                        wt.values.set(i, j, wt_temp);
                    }

                if (!MPIc->rank())
                {
                    IMPACT_Write(&wt, name + "_", dir + name + "/", timestep);
                }
            }
        }

        // ****************************************
        // New - now output complete distribution functions

        if (!ifdumpdists)
        {
            dir = IMPACT_Messages::Data_Directory + IMPACT_Messages::DistFunc_Dir;
            if (if_dump_switches::dump_f0)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                IMPACT_Dist f0_out(c, 0);
                f0_out.fill(v, c, MPIc);
                name = "f0";
                if (!MPIc->rank())
                {
                    IMPACT_Write(&f0_out, name + "_", dir + name + "/", timestep);
                }
            }
            if (if_dump_switches::dump_f1 && c->Nf1() > 0)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                IMPACT_Dist f1_out(c, 1);
                f1_out.fill(v, c, MPIc);
                name = "f1";
                if (!MPIc->rank())
                {
                    IMPACT_Write(&f1_out, name + "_", dir + name + "/", timestep);
                }
            }
            if (if_dump_switches::dump_f2 && c->Nf2() > 0)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                IMPACT_Dist f2_out(c, 2);
                f2_out.fill(v, c, MPIc);
                name = "f2";
                if (!MPIc->rank())
                {
                    IMPACT_Write(&f2_out, name + "_", dir + name + "/", timestep);
                }
            }
        }
        // ****************************************

        timeend = clock();

        std::ostringstream Imessage2;

        Imessage2 << "\n"
                  << ENDFORMAT << "\nIMPACT: Data dump - Total time = " << double(timeend - timestart) / CLOCKS_PER_SEC << " s" << std::endl
                  << ULINE;

        if (!MPIc->rank())
            std::cout << Imessage2.str();
    }

    // Put back quantities to original values
    IMPACTA_Rescale_T(c, MPIc, v, f0, f1, f2, f3, E, B, equation_switches::NEW_T_NORM);
}
