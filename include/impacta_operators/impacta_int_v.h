/*
**********************************************************
 Contains stencils for integrating over v - for collision 
 operators etc.

Version 1.0.0
1
AGRT
21/2/07

**********************************************************
*/

class IMPACT_vint //class for inserting integrals over v in matrix
{
 private:
  double *int_terms; // store terms
  int Numv; // number of v points
 public:
  IMPACT_vint() //default constructor
    {
      int_terms=new double[1];
      Numv=1;
    }

  IMPACT_vint(int N)
    {
      int_terms=new double[N];
      Numv=N;
    }
  IMPACT_vint(IMPACT_Config *c) // takes size from Nv() in config
    {
      int_terms=new double[c->Nv()+2]; // +2 to be same as v,dv arrays
      Numv=c->Nv();
    }
  ~IMPACT_vint() //default destructor
    {
      delete[] int_terms;
    }
  //resize vector;
 void resize(IMPACT_Config *c)
    {
      delete[] int_terms;
      int_terms=new double[c->Nv()+2]; // +2 to be same as v,dv arrays
      Numv=c->Nv();
    }
  //access method to get particular element
  double get(int k)
  {
    return int_terms[k]; 
  }
  //access method to set particular element
  void set(int k,double val)
  {
    int_terms[k]=val; 
  }
  int Nv()
  {
    return Numv;
  }
  void print()
  {
    for (int k=0;k<Numv+2;++k)
      std::cout<<int_terms[k]<<",";
    std::cout<<std::endl;
  }
 
};
