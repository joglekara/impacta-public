/*
**********************************************************
Obtain distribution functions of the equations

Version 1.0
AGRT

25/7/07
**********************************************************
*/

class IMPACT_Dist
{
 private:
  struct IMPACT_Moment **dist;
  int num_components; // the number of components of this distribution function
  int num_v; // the number of velocity cells
  double* _dv;
  int num_x;
  int num_y;
  int f_;
 public:
  IMPACT_Dist()
    {
      num_components=0;
      num_v=num_x=num_y=0;
      f_=0;
      _dv=new double[1];
      dist=new struct IMPACT_Moment*[1];
      dist[0]=new struct IMPACT_Moment[num_v];
    }
  virtual ~IMPACT_Dist()
    {
      for (int i=0;i< num_components;++i)
	delete[] dist[i];
      delete[] dist;
      delete[] _dv;
    }
  IMPACT_Dist(IMPACT_Config *c, int fin) //f is the f of dist (i.e. f1,f2 etc.)
    {
      num_components = c->Nf(fin);
      num_v=c->Nv();
      _dv=new double[num_v+1];
      for (int k=0;k<=num_v;++k)
	_dv[k]=c->dv(&k);
      num_x=c->Nx();
      num_y=c->Ny();
      dist=new struct IMPACT_Moment*[num_components];
      for (int tt=0;tt<num_components;++tt)
	dist[tt]=new struct IMPACT_Moment[num_v];
      f_=fin;
      std::string namestr = "f_";
      std::ostringstream numtemp;
      numtemp<<fin;
      namestr+=numtemp.str();
      for (int tt=0;tt<num_components;++tt)
      for (int k=0;k<num_v;++k){
	strcpy(dist[tt][k].name,namestr.c_str());
	strcpy(dist[tt][k].coords,IMPACT_Coords);
	dist[tt][k].Nx=num_x;
	dist[tt][k].Ny=num_y;
	dist[tt][k].values.ChangeSize(num_x,num_y);
	GetGrid(c,&dist[tt][k]);
      }
    }
  // Now fill it with values

  inline void fill(IMPACT_ParVec *v, IMPACT_Config *c, IMPACT_MPI_Config *MPIc)
  {
    double *data;
    data = new double[c->cellpoints()];
    for (int tt=0;tt<num_components;++tt)
      {
	int kstart=0;
	if (f_>0)
	  {
	  for (int temp=0;temp<f_;++temp)
	    kstart+=c->Nv()*c->Nf(temp);
	  kstart+=tt*c->Nv(); //add on the cells from previous components
	  }
	else kstart=1;
	for (int i=1;i<=num_x;++i)
	  {
	    for (int j=1;j<=num_y;++j)
	      {
		Gather_kstring(v,c,MPIc,&i,&j,data);
		for (int k=0;k<num_v;++k)
		  dist[tt][k].values.set(i,j,data[k+kstart-1]);
	      }
	  } //end of i loop
      }//end of tt loop
    delete[] data;
  }
  // Access methods
  struct IMPACT_Moment* getadd(int component, int k)
  {
    return &dist[component][k-1];
  }
  int Nv()
  {
    return num_v;
  }
  double dv(int *k)
  {
    return _dv[*k];
  }
  int N()
  {
    return num_components;
  }
  int f()
  {
    return f_;
  }
  std::string comp_str(int component)
    {
      std::string answer="";
      switch (f_)
	{
	case 0:
	  break;
	case 1:
	  answer = GetDirChar(component+1); break;
	case 2:
	  answer = Get2DirChar(component+1); break;
	}
      return answer;
    }
};
//std::ofstream& operator << (std::ofstream& str_out, IMPACT_Moment *IM)
void operator<< (std::ofstream& str_out, IMPACT_Dist *IM)
{
  // Set precision of output
  str_out.precision(IMPACT_Diagnostics::output_precision);
  str_out.setf(std::ios::scientific,std::ios::floatfield); 

  //Header for data file.
  str_out << IMPACT_Messages::DataHeader[1]<<'\n';
  str_out << IMPACT_Messages::Version<<TAB<<":"<<TAB<< IMPACT_Messages::Date<<'\n';
  
  str_out <<"0\n"; //marks end of header
  
  str_out << IM->getadd(0,1)->name <<'\n';
  str_out <<IM->getadd(0,1)->coords[0]<<IM->getadd(0,1)->coords[1]<<'\n';
  // Now print size of 3 grid
  str_out <<IM->getadd(0,1)->Nx<<TAB<<IM->getadd(0,1)->Ny<<TAB<<IM->Nv()<<'\n';
   // Now print axes 
  for (int i=1;i<=IM->getadd(0,1)->Nx;++i)
    str_out << IM->getadd(0,1)->xaxis.Get(i)<<TAB;
  str_out<<'\n';
  for (int i=1;i<=IM->getadd(0,1)->Ny;++i)
    str_out << IM->getadd(0,1)->yaxis.Get(i)<<TAB;
  str_out<<'\n';
  double vsum=0.0; int iminus;
  for (int i=1;i<=IM->Nv();++i)
    {
      iminus=i-1;
      vsum+=(IM->dv(&i)+IM->dv(&iminus))*0.5;
      str_out << vsum<<TAB;
      
    }
  str_out<<'\n';
    // END OF HEADER INFORMATION
  for (int tt=0;tt<IM->N();++tt)
    {
      str_out <<'\n'<<ULINE<<"f_"<<IM->f()<<(IM->comp_str(tt)).c_str()<<'\n';
    for (int i=1;i<=IM->getadd(0,1)->Nx;++i)
      for (int j=1;j<=IM->getadd(0,1)->Ny;++j)
	{
	  str_out << "\n("<<i<<','<<j<<")\n";
	  for (int k=1;k<=IM->Nv();++k)
	    str_out << IM->getadd(tt,k)->values.get(i,j)<<TAB;
     	}
	str_out << '\n';
    }
}

/*
//std::ofstream& operator >> (std::ifstream& str_in, IMPACT_Moment *IM)
void operator >> (std::ifstream& str_in, IMPACT_Moment *IM)
{
  char temp[20];
  //First skip header
  do{
  str_in >> temp;
  }while (temp[0]!='0');
  str_in >> IM->name;
  str_in >>IM->coords[0];
  str_in>>IM->coords[1];
  str_in >>IM->Nx>>IM->Ny;
  
  //Resize the vectors in the moment appropriately
  IM->xaxis.ChangeSize(IM->Nx);
  IM->yaxis.ChangeSize(IM->Ny);
  IM->values.ChangeSize(IM->Nx,IM->Ny);
 
  //Now put in values
   double tempx;
  //IM->xaxis = new double[IM->Nx];
  for (int i=0;i<IM->Nx;++i)
    { str_in >> tempx;
      IM->xaxis.Set(i+1,tempx);
    }
   double tempy;
  //IM->yaxis = new double[IM->Ny];
  for (int i=0;i<IM->Ny;++i)
    { str_in >> tempy;
      IM->yaxis.Set(i+1,tempy);
    }
  
  double tempvalues;
  for (int i=0;i<IM->Nx;++i)
      for (int j=0;j<IM->Ny;++j)
	{
	str_in >> tempvalues;
	IM->values.set(i+1,j+1,tempvalues);
	}
	};


*/
