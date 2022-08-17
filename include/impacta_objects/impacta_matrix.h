/*
**********************************************************
IMPACT Version 1.0
Matrix and Sparse matrix classes

Version 1.3
AGRT

19/2/07

Latest notes:
Added method for inserting vector into sparse matrix at
row rownumber
IMPORTANT - this can be made much more efficient possibly in future.

Also - stencil class added - this allows ddx stencils etc.

5/3 - Had to change the count member function. This is because
the order of the rownumber calling was critical previously

!!!IMPORTANT!!!!
Resize MUST be  called after all rows have been counted as it now
sorts out the num in rows which WILL BE COMPLETELY WRONG OTHERWISE!!!!

7/3 - New count routines for each equation which count only the relevent
elements

18/4/07 - new added velocity stencil - 3 point for inserting
velocity differentials

8/5 - velocity stencil array added for storing RB coefficents
**********************************************************

*/

class IMPACT_Matrix // matrix class
{
private:
    friend class IMPACT_Sparse;
    IMPACT_Vector **Matrix;
    int N, M; // number of rows/cols
public:
    IMPACT_Matrix();
    IMPACT_Matrix(int Nin, int Min);      // Square Matrix_size - -self explanatory			//note - no default constructor!!
    IMPACT_Matrix(IMPACT_Config config1); // Square Matrix_size - -self explanatory			//note - no default constructor!!
    virtual ~IMPACT_Matrix();
    double get(int, int);
    void set(int, int, double);
    void inc(int, int, double);
    void setall(double);
    void Print();
    void Print(std::ofstream &outfile);
    void ChangeSize(int new_N, int new_M);
    void multiply(double value);
    void size(int *Nout, int *Mout);
    int chksize(int *Nout, int *Mout);
    void reset();
};
inline IMPACT_Matrix::IMPACT_Matrix()
{
    //(adapted from Vector)
    Matrix = new IMPACT_Vector *[1];
    Matrix[0] = new IMPACT_Vector();
    N = M = 1;
}
inline IMPACT_Matrix::IMPACT_Matrix(int Nin, int Min)
{ // N by M matrix
    //(adapted from Vector)
    Matrix = new IMPACT_Vector *[Nin];
    for (int i = 0; i < Nin; i++)
        Matrix[i] = new IMPACT_Vector(Min);
    N = Nin, M = Min;
}
inline IMPACT_Matrix::IMPACT_Matrix(IMPACT_Config config1)
{
    N = (config1.conf_N_f0 + config1.conf_N_f1 + config1.conf_N_f2 + config1.conf_N_f3) * config1.conf_Nv + config1.conf_NE + config1.conf_NB;
    M = N = N * config1.conf_Nx * config1.conf_Ny;
    Matrix = new IMPACT_Vector *[N];
    for (int i = 0; i < N; i++)
        Matrix[i] = new IMPACT_Vector(N);
}
inline IMPACT_Matrix::~IMPACT_Matrix()
{ // This deallocates memory resources
    for (int i = 0; i < N; i++)
        delete Matrix[i]; // When the object goes out of scope
    delete[] Matrix;
}
inline double IMPACT_Matrix::get(int i, int j)
{ // i and j are from 1 to N not 0 to N-1
    double answer;
    answer = Matrix[i - 1]->Get(j);
    return answer;
}
inline void IMPACT_Matrix::set(int i, int j, double value)
{
    // i and j are from 1 to N not 0 to N-1
    Matrix[i - 1]->Set(j, value);
}
inline void IMPACT_Matrix::inc(int i, int j, double value)
{
    // i and j are from 1 to N not 0 to N-1
    Matrix[i - 1]->Set(j, value + get(i, j));
}
inline void IMPACT_Matrix::setall(double value)
{
    // i and j are from 1 to N not 0 to N-1
    for (int i = 0; i < N; i++)
        Matrix[i]->SetAll(value);
}
inline void IMPACT_Matrix::ChangeSize(int Nin, int Min)
{
    // This deallocates memory resources
    for (int i = 0; i < N; i++)
        delete Matrix[i]; // When the object goes out of scope
    delete[] Matrix;

    Matrix = new IMPACT_Vector *[Nin];
    for (int i = 0; i < Nin; i++)
        Matrix[i] = new IMPACT_Vector(Min);
    N = Nin, M = Min;
}
inline void IMPACT_Matrix::multiply(double value)
{
    for (int i = 1; i <= N; i++)
        for (int j = 1; j <= M; j++)
            set(i, j, value * get(i, j));
}
inline void IMPACT_Matrix::size(int *Nout, int *Mout)
{
    *Nout = N;
    *Mout = M;
}
inline int IMPACT_Matrix::chksize(int *Nout, int *Mout)
{
    int ans = 0;
    ans += (*Nout == N);
    ans += (*Mout == M);
    ans = ans / 2;
    return ans;
}
inline void IMPACT_Matrix::Print()
{
    // Needs no explanation?
    for (int i = 1; i <= N; i++)
    {
        std::cout << i << ": ";
        for (int j = 1; j <= M; j++)
        {
            if (i == j)
                std::cout << BCYAN;
            else
                std::cout << ENDFORMAT;
            std::cout << get(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl
              << std::endl
              << std::endl;
    std::cout << ENDFORMAT;
}
inline void IMPACT_Matrix::Print(std::ofstream &outfile)
{
    // Needs no explanation?
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= M; j++)
        {
            outfile << get(i, j) << " ";
        }
        outfile << std::endl;
    }
    outfile << std::endl
            << std::endl
            << std::endl;
}
inline void IMPACT_Matrix::reset()
{
    for (int i = 0; i < N; ++i)
        Matrix[i]->reset();
}
class IMPACT_Sparse
{
    // These functions are for efficient counting and packing but need Var
    // inputs so have to be outside the class.
    /*  friend void f0equation_count(IMPACT_ParVec * v, IMPACT_Sparse *S,IMPACT_Config *c, IMPACT_Var * f0,IMPACT_Var * E,IMPACT_Var * f1,int*i,int*j,int*k);
    friend void f1equation_count(IMPACT_ParVec * v, IMPACT_Sparse *S,IMPACT_Config *c, IMPACT_Var * f0,IMPACT_Var * E,IMPACT_Var * B,IMPACT_Var * f1,int*i,int*j,int*k,IMPACT_Dim * x1);
    friend void Eequation_count(IMPACT_ParVec * v, IMPACT_Sparse *S,IMPACT_Config *c,IMPACT_Var * E,IMPACT_Var *B,IMPACT_Var * f1,int*i,int*j,IMPACT_Dim * x1);
    friend void Bequation_count(IMPACT_ParVec * v, IMPACT_Sparse *S,IMPACT_Config *c,IMPACT_Var * E,IMPACT_Var *B,int*i,int*j,IMPACT_Dim * x1);*/
    // friend void PackSparse_int(IMPACT_Sparse *S,IMPACT_ParVec *v,int packstart,int packend,int *runningtotal);
    // friend void PackSparse(IMPACT_Sparse *S,IMPACT_ParVec *v,int i,int *runningtotal);
    // friend void chksparse_chk(IMPACT_Sparse *S, int *rownumber,int *runningtotal);
private:
    IMPACT_Vector *values; // the values to be input
    int *numinrows;        // array of the start of the rowsin the matrix
    int *colnums;          // the column labels
    int N;                 // the dimensions of the NxN matrix
    int sparselength;      // length of the sparse matrix
public:
    IMPACT_Sparse()
    {
        values = new IMPACT_Vector(1);
        numinrows = new int[2];
        colnums = new int[1];
        sparselength = 1;
        N = 1;
    }
    virtual ~IMPACT_Sparse()
    {
        delete values;
        delete[] numinrows;
        delete[] colnums;
    }
    IMPACT_Sparse(IMPACT_Config config1)
    {
        N = (config1.conf_N_f0 + config1.conf_N_f1 + config1.conf_N_f2 + config1.conf_N_f3) * config1.conf_Nv + config1.conf_NE + config1.conf_NB;
        N *= config1.conf_Nx * config1.conf_Ny;
        numinrows = new int[N + 1];
        colnums = new int[1];
        values = new IMPACT_Vector(1);
        sparselength = 1;
    }
    IMPACT_Sparse(IMPACT_Config config1, int slength) // the length of sparse is defined
    {
        N = (config1.conf_N_f0 + config1.conf_N_f1 + config1.conf_N_f2 + config1.conf_N_f3) * config1.conf_Nv + config1.conf_NE + config1.conf_NB;
        N *= config1.conf_Nx * config1.conf_Ny;
        numinrows = new int[N + 1];
        colnums = new int[slength];
        values = new IMPACT_Vector(slength);
        sparselength = slength;
    }
    // This next one is to enable the parallel sparse matrix
    IMPACT_Sparse(int rowlength)
    // the length of the rowvector is defined
    {
        N = rowlength;
        numinrows = new int[rowlength + 1];
        colnums = new int[1];
        values = new IMPACT_Vector(1);
        sparselength = 1;
    }
    inline int getN()
    {
        return N;
    }
    inline int getSl()
    {
        return sparselength;
    }
    void printvals()
    {
        values->PrintArray();
    }
    void printrows()
    {
        for (int i = 0; i < N; i++)
            std::cout << "(" << i + 1 << "," << numinrows[i] << "),";
        std::cout << std::endl;
    }
    void printvalsbyrow()
    {
        for (int i = 0; i < N; i++)
        {
            std::cout << "Row " << i + 1 << ": ";
            for (int j = numinrows[i]; j < numinrows[i + 1]; j++)
                std::cout << values->Get(j) << ",";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    void printelembyrow()
    {
        for (int i = 0; i < N; i++)
        {
            std::cout << "Row " << i + 1 << ": ";
            for (int j = numinrows[i]; j < numinrows[i + 1]; j++)
                std::cout << j << ",";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    void printcols()
    {
        for (int i = 0; i < sparselength; i++)
            std::cout << colnums[i] << ",";
        std::cout << std::endl;
    }
    void printcolsbyrow()
    {
        for (int i = 0; i < N; i++)
        {
            std::cout << "Row " << i + 1 << ": ";
            for (int j = numinrows[i]; j < numinrows[i + 1]; j++)
                std::cout << colnums[j - 1] << ",";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    inline int getrow(int index) // returns value in row index
    {
        // error checking
        //  if(index<1||index>N+1){std::cout<<"in Sparse getrow - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        return numinrows[index - 1];
    }
    inline int getrow_nocheck(int index) // returns value in row index
    {
        // error checking
        //  if(index<0||index>N+1){std::cout<<"in Sparse getrow - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        return numinrows[index - 1];
    }
    inline void setrow(int index, int value) // sends value to row index
    {
        // error checking
        //  if(index<1||index>N+1){std::cout<<"in Sparse setrow - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        numinrows[index - 1] = value;
    }
    inline void incrow(int index, int value) // sends value incrimentally to row index
    {
        // error checking
        //  if(index<1||index>N+1){std::cout<<"in Sparse incrow - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        numinrows[index - 1] += value;
    }
    inline int getcol(int index) // returns value in col index
    {
        // error checking
        //  if(index<1||index>sparselength){std::cout<<"in Sparse getcol - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        return colnums[index - 1];
    }
    inline int *getcol_add(int index) // returns value in col index
    {
        // error checking
        //  if(index<1||index>sparselength){std::cout<<"in Sparse getcol - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        return &colnums[index - 1];
    }
    inline void setcol(int index, int value) // sends value to col index
    {
        // if(index<1||index>sparselength){std::cout<<"in Sparse setcol - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        colnums[index - 1] = value;
    }
    inline double getval(int index) // returns value in val index
    {
        // error checking
        //  if(index<1||index>sparselength){std::cout<<"in Sparse getval - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        return values->Get(index);
    }
    inline double *getval_add(int index)
    {
        // error checking
        // if(index<1||index>sparselength){std::cout<<"in Sparse getval - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        return values->Get_add(index);
    }
    inline void setval(int index, double value) // sends value to val index
    {
        // if(index<1||index>sparselength){std::cout<<"in Sparse setval - index outside max value"<<std::endl<<"element = "<< index<<std::endl;exit(0);}
        values->Set(index, value);
    }
    inline void Pack(IMPACT_Matrix *M) // pack a matrix into the sparse matrix
    {
        if (N != M->N) // If Matrix dimensions do not agree then
        {
            N = M->N; // resize sparse
            delete[] numinrows;
            numinrows = new int[N + 1];
        }
        //
        double zerothresh = zerotolerance::zerothreshold;
        double val;
        int total = 0;
        numinrows[0] = 1;            // first element is at 1!
        for (int i = 1; i <= N; ++i) // This bit counts non-zeros
        {
            numinrows[i] = numinrows[i - 1];
            for (int j = 1; j <= N; ++j)
            {
                val = M->get(i, j);
                if (fabs(val) > zerothresh)
                    ++numinrows[i];
            }

            total += (numinrows[i] - numinrows[i - 1]);
        }
        // Now we can make new vectors for values and colnums
        if (total != sparselength)
        {
            sparselength = total;
            delete[] colnums;
            values->ChangeSize(total);
            colnums = new int[total];
        }

        int runningtotal = 0;        // to keep tabs of the
        for (int i = 1; i <= N; ++i) // Now pack sparse matrix
        {
            for (int j = 1; j <= N; ++j)
            {
                val = M->get(i, j);
                if (val > zerothresh)
                {
                    values->Set(runningtotal + 1, val);
                    colnums[runningtotal] = j;
                    runningtotal++;
                }
            } // j
        }     // i
    }
    inline void reset() // this makes everything zero EXCEPT the rownums
    {
        for (int i = 1; i <= sparselength; ++i)
        {
            colnums[i - 1] = 0;
            values->Set(i, 0.0);
        }
    }
    inline int getfree(int rownumber)
    // finds free spot in sparse matrix at row rownumber
    {
        int freespot = 0;
        int ii = getrow(rownumber);
        while (freespot == 0)
        {
            if (getcol(ii) == 0)
                freespot = ii;
            // if (ii>S->getrow(rownumber+1))
            //   {std::cout<<"Error, no free space in sparse matrix";exit(0);}
            ++ii;
        }
        return freespot;
    }
    //_________________________________________________________________
    // Counts non zeros in a row and packs them:

    inline void Pack(IMPACT_ParVec *v, int rownumber)
    {
        int index = getrow(rownumber);
        double zerothresh = zerotolerance::zerothreshold;
        double val;
        int runningtotal = index - 1;          // to keep tabs of the
        for (int i = 1; i <= v->length(); ++i) // Now pack sparse matrix
        {
            val = v->Get(i);
            if (fabs(val) > zerothresh)
            {
                values->Set(runningtotal + 1, val);
                colnums[runningtotal] = i;
                runningtotal++;
            }
        } // i
    }

    // this counts the non-zero elements in *v and resizes accordingly
    //  NB!!! Probably needs a reset after to set to zero
    inline void Count(IMPACT_ParVec *v, int rownumber)
    {
        double zerothresh = zerotolerance::zerothreshold;
        double val;
        // store old rowcount

        numinrows[rownumber] = 0; // numinrows[rownumber-1];
        /*
          NOTE NOTE NOTE This means numinrows is WRONG until
          Resize is called!!
         */
        for (int j = 1; j <= N; ++j)
        { // cycle through and count

            val = v->Get(j);
            if (fabs(val) > zerothresh)
                ++numinrows[rownumber];
        }
    }
    inline void Resize() // after count, resizes column and row vectors
    {
        numinrows[0] = 1; // first element is at 1!
        for (int j = 0; j < N; ++j)
        {
            numinrows[j + 1] += numinrows[j]; // adds rows correctly
        }
        sparselength = (numinrows[N] - 1); //-1 as next(non-existent) row would be at N
        // updates length of sparse vector
        delete[] colnums;
        values->ChangeSize(sparselength);

        colnums = new int[sparselength];
    }
};

IMPACT_Matrix S2M(IMPACT_Sparse *S) // for checking
{

    int N = S->getN();
    IMPACT_Matrix M(N, N);
    for (int i = 1; i <= N; ++i)
        for (int j = S->getrow(i); j < S->getrow(i + 1); ++j)
        {
            if (S->getcol(j) != 0)
                M.set(i, S->getcol(j), S->getval(j));
        }
    return M;
}
int getfree(int rownumber, IMPACT_Sparse *S)
// finds free spot in sparse matrix at row rownumber
{
    int freespot = 0;
    int ii = S->getrow(rownumber);
    while (freespot == 0)
    {
        if (S->getcol(ii) == 0)
            freespot = ii;
        // if (ii>S->getrow(rownumber+1))
        //   {std::cout<<"Error, no free space in sparse matrix";exit(0);}
        ++ii;
    }
    return freespot;
}
/*
This next class has 5 components (can obviously be extended if required)
These represent a point in x,y space and the neighbouring points.
e.g. (x,y),(x+1,y),(x-1,y), (x,y+1),(x,y-1)
 */
class IMPACT_stencil
{
private:
    double stencil[5]; // (x,y),(x+1,y),(x-1,y), (x,y+1),(x,y-1);
public:
    IMPACT_stencil()
    {
        for (int i = 0; i < 5; ++i)
        {
            stencil[i] = 0;
        }
    }
    ~IMPACT_stencil()
    {
    }
    IMPACT_stencil(double s1, double s2, double s3, double s4, double s5)
    {
        stencil[0] = s1;
        stencil[1] = s2;
        stencil[2] = s3;
        stencil[3] = s4;
        stencil[4] = s5;
    }
    /*  IMPACT_stencil(double s1,double s2,double s3,double s4,double s5,double m)
      {
        stencil[0]=s1;
        stencil[1]=s2;
        stencil[2]=s3;
        stencil[3]=s4;
        stencil[4]=s5;
        multiplier=m;
        }*/
    /* IMPACT_stencil(double * s,double m)
      {
     for (int i=0;i<5;++i)
      {
        stencil[i]=s[i];
      }
        multiplier=m;
        }*/
    inline double operator()(int index)
    {
        return stencil[index];
    }
    inline double get(int index)
    {
        return stencil[index];
    }
    inline void set(int index, double val)
    {
        stencil[index] = val;
    }
    inline void invert(int index)
    {
        stencil[index] = -stencil[index];
    }
    inline IMPACT_stencil operator*(double m) // updates the multiplier
                                              // returns a stencil so that it can be used
    {
        // multiplier=m;
        IMPACT_stencil answer(stencil[0] * m, stencil[1] * m, stencil[2] * m, stencil[3] * m, stencil[4] * m);
        return answer;
    }
    inline IMPACT_stencil operator+(IMPACT_stencil toadd) // adds stencils
    {
        IMPACT_stencil answer(get(0) + toadd.get(0), get(1) + toadd.get(1), get(2) + toadd.get(2), get(3) + toadd.get(3), get(4) + toadd.get(4));
        return answer;
    }
    inline void flip(IMPACT_Dim *x1)
    {
        double tempstore;
        switch (x1->get())
        {
        case 1:
            tempstore = stencil[1];
            stencil[1] = stencil[2];
            stencil[2] = tempstore;
            break;
        case 2:
            tempstore = stencil[3];
            stencil[3] = stencil[4];
            stencil[4] = tempstore;
            break;
        }
    }
    void print()
    {
        for (int i = 0; i < 5; ++i)
            std::cout << stencil[i] << ',';
        std::cout << '\n';
    }
};

class IMPACT_Vel_Sten // velocity space stencil
{
private:
    double stencil[3]; // (v-1,v,v+1)
public:
    IMPACT_Vel_Sten()
    {
        stencil[0] = stencil[1] = stencil[2] = 0.0;
    }
    ~IMPACT_Vel_Sten()
    {
    }
    IMPACT_Vel_Sten(double s1, double s2, double s3)
    {
        stencil[0] = s1;
        stencil[1] = s2;
        stencil[2] = s3;
    }
    inline double get(int index)
    {
        return stencil[index];
    }
    inline void set(int index, double val)
    {
        stencil[index] = val;
    }
    inline void inc(int index, double val)
    {
        stencil[index] += val;
    }
    inline void operator+=(IMPACT_Vel_Sten toadd)
    {
        stencil[0] = toadd.get(0);
        stencil[1] = toadd.get(1);
        stencil[2] = toadd.get(2);
    }
    inline void operator*(double value)
    {
        for (int i = 0; i < 3; ++i)
            stencil[i] *= value;
    }
    inline void printelements()
    {
        std::cout << "(" << stencil[0] << ',' << stencil[1] << ',' << stencil[2] << ")\n";
    }
};

class IMPACT_VS_Array // velocity space stencil ARRAY
{
private:
    IMPACT_Vel_Sten ****VSten;
    int Ny, Nv, xstart, Nx; // as parallel=>x start to end
public:
    IMPACT_VS_Array()
    {
        Ny = Nv = xstart = Nx = 1;
        VSten = new IMPACT_Vel_Sten ***[Nx]; // x
        for (int i = 0; i < Nx; ++i)
        {
            VSten[i] = new IMPACT_Vel_Sten **[Ny + 1];
            for (int j = 0; j <= Ny; ++j)
            {
                VSten[i][j] = new IMPACT_Vel_Sten *[Nv + 1];
                for (int k = 0; k <= Nv; ++k)
                    VSten[i][j][k] = new IMPACT_Vel_Sten(0.0, 0.0, 0.0);
            }
        }
    }
    ~IMPACT_VS_Array()
    {
        delete[] VSten;
    }
    IMPACT_VS_Array(IMPACT_Config *c, IMPACT_MPI_Config *M)
    {
        Ny = c->Ny();
        Nv = c->Nv();
        xstart = M->istart();
        Nx = M->iend() + 1 - xstart;
        VSten = new IMPACT_Vel_Sten ***[Nx]; // x
        for (int i = 0; i < Nx; ++i)
        {
            VSten[i] = new IMPACT_Vel_Sten **[Ny + 1];
            for (int j = 0; j <= Ny; ++j)
            {
                VSten[i][j] = new IMPACT_Vel_Sten *[Nv + 1];
                for (int k = 0; k <= Nv; ++k)
                    VSten[i][j][k] = new IMPACT_Vel_Sten(0.0, 0.0, 0.0);
            }
        }
    }
    void set(IMPACT_Vel_Sten *vs, int *i, int *j, int *k)
    {
        for (int u = 0; u < 3; ++u)
            VSten[*i - xstart][*j][*k]->set(u, vs->get(u));
    }
    void inc(IMPACT_Vel_Sten *vs, int *i, int *j, int *k)
    {
        for (int u = 0; u < 3; ++u)
            VSten[*i - xstart][*j][*k]->inc(u, vs->get(u));
    }
    void inc_kminus(IMPACT_Vel_Sten *vs, int *i, int *j, int *k)
    {
        for (int u = 1; u < 3; ++u)
            VSten[*i - xstart][*j][*k]->inc(u, vs->get(u - 1));
    }
    void multiply(double *value, int *i, int *j, int *k)
    {
        for (int u = 0; u < 3; ++u)
            VSten[*i - xstart][*j][*k]->set(u, VSten[*i - xstart][*j][*k]->get(u) * (*value));
    }
    IMPACT_Vel_Sten *get(int *i, int *j, int *k)
    {
        return VSten[*i - xstart][*j][*k];
    }
};
class IMPACT_IMatrix : public IMPACT_Matrix
{
private:
    double constvalue;
    int ifconst;

public:
    IMPACT_IMatrix() : IMPACT_Matrix()
    {
        constvalue = 1;
        ifconst = 1;
    }
    IMPACT_IMatrix(double value) : IMPACT_Matrix()
    {
        constvalue = value;
        ifconst = 1;
    }
    IMPACT_IMatrix(double *gridx, double *gridy, int Nx, int Ny) : IMPACT_Matrix(Nx, Ny)
    {
        constvalue = 0;
        ifconst = 0;
        double tempinsert;
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
            {
                tempinsert = gridx[i] * gridy[j];
                IMPACT_Matrix::set(i + 1, j + 1, tempinsert);
            }
    }
    void copy(IMPACT_IMatrix *IM)
    {
        int Nin, Min;
        IM->size(&Nin, &Min);
        if (Nin * Min == 1)
            setc(IM->get(&Nin, &Min));
        else
        {
            ifconst = 0;
            if (!chksize(&Nin, &Min))
                IMPACT_Matrix::ChangeSize(Nin, Min);
            for (int i = 1; i <= Nin; ++i)
                for (int j = 1; j <= Min; ++j)
                    IMPACT_Matrix::set(i, j, IM->get(&i, &j));
        }
        /*if (Nin*Min==1) {ifconst=1; constvalue=IM->get(&Nin,&Min);}
        else
          {
        ifconst=0;
        IMPACT_Matrix::ChangeSize(Nin,Min);
        for (int i=1;i<=Nin;++i)
          for (int j=1;j<=Min;++j)
          set(i,j,IM->get(&i,&j));}*/
    }
    void copy(IMPACT_Matrix *IM)
    {
        int Nin, Min;
        IM->size(&Nin, &Min);
        if (Nin * Min == 1)
            setc(IM->get(Nin, Min));
        else
        {
            ifconst = 0;
            if (!chksize(&Nin, &Min))
                IMPACT_Matrix::ChangeSize(Nin, Min);
            for (int i = 1; i <= Nin; ++i)
                for (int j = 1; j <= Min; ++j)
                    IMPACT_Matrix::set(i, j, IM->get(i, j));
        }
        /*if (Nin*Min==1) {ifconst=1; constvalue=IM->get(&Nin,&Min);}
        else
          {
      ifconst=0;
      IMPACT_Matrix::ChangeSize(Nin,Min);
      for (int i=1;i<=Nin;++i)
        for (int j=1;j<=Min;++j)
        set(i,j,IM->get(&i,&j));}*/
    }
    void setc(double value) // set a constant value
    {
        constvalue = value;
        ifconst = 1;
    }
    void setg(double *gridx, double *gridy, int Nx, int Ny)
    {
        constvalue = 0;
        ifconst = 0;
        double tempinsert;
        IMPACT_Matrix::ChangeSize(Nx, Ny);
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
            {
                tempinsert = gridx[i] * gridy[j];
                IMPACT_Matrix::set(i + 1, j + 1, tempinsert);
            }
    }
    double get(int *i, int *j)
    {
        double result = constvalue;
        if (!ifconst)
            result = IMPACT_Matrix::get(*i, *j);
        return result;
    }
    void Iset(double value, int *i, int *j)
    {
        constvalue = value;
        if (!ifconst)
            IMPACT_Matrix::set(*i, *j, value);
    }
    void multiplyall(double value)
    {
        if (ifconst)
            constvalue *= value;
        if (!ifconst)
            IMPACT_Matrix::multiply(value);
    }
    void Print()
    {
        if (!ifconst)
            IMPACT_Matrix::Print();
        if (ifconst)
            std::cout << constvalue << '\n';
    }
    void Printf(std::ofstream &outfile)
    {

        if (!ifconst)
            IMPACT_Matrix::Print(outfile);
        if (ifconst)
            std::cout << constvalue << '\n';
    }
    void SwitchOffConst()
    {
        ifconst = 0;
    }
    void SwitchOffConst(int Nxin, int Nyin)
    {
        if (ifconst)
        {
            ifconst = 0;
            ChangeSize(Nxin, Nyin);
            for (int i = 1; i <= Nxin; ++i)
                for (int j = 1; j <= Nyin; ++j)
                    Iset(constvalue, &i, &j);
        }
    }
    void SwitchOnConst()
    {
        ifconst = 1;
    }
    bool sign(int *i, int *j)
    {
        return (get(i, j) >= 0.0);
    }
};
void IMPACTA_Share_Moment(IMPACT_Config *c, IMPACT_MPI_Config *M,
                          IMPACT_IMatrix *Moment)
{
    if (M->size() > 1)
    {
        double *jstring;
        jstring = new double[c->Ny()];

        for (int i = 1; i <= c->Nx(); ++i)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            for (int j = 1; j <= c->Ny(); ++j)
                jstring[j - 1] = Moment->get(&i, &j);
            // MPI::COMM_WORLD.Bcast(jstring,c->Ny(),MPI::DOUBLE,M->rank(&i));
            MPI_Bcast(jstring, c->Ny(), MPI_DOUBLE, M->rank(&i), MPI_COMM_WORLD);
            for (int j = 1; j <= c->Ny(); ++j)
                Moment->set(i, j, jstring[j - 1]);
        }
        delete[] jstring;
    }
}