/*
**********************************************************
Input and Output of Moments of the equations

Version 1.0
AGRT

23/3/07
**********************************************************
*/

//std::ofstream& operator << (std::ofstream& str_out, IMPACT_Moment *IM)
void operator<< (std::ofstream& str_out, IMPACT_Moment *IM)
{
  
  //Header for data file.
  str_out << IMPACT_Messages::DataHeader[0]<<'\n';
  str_out << IMPACT_Messages::Version<<TAB<<":"<<TAB<< IMPACT_Messages::Date<<'\n';
  
  str_out <<"0\n"; //marks end of header
  
  str_out << IM->name <<'\n';
  str_out <<IM->coords[0]<<IM->coords[1]<<'\n';
  str_out <<IM->Nx<<TAB<<IM->Ny<<'\n';

// Set precision of output
  str_out.precision(IMPACT_Diagnostics::output_precision);
  str_out.setf(std::ios::scientific,std::ios::floatfield); 

  for (int i=1;i<=IM->Nx;++i)
    str_out << IM->xaxis.Get(i)<<TAB;
  str_out<<'\n';
  for (int i=1;i<=IM->Ny;++i)
    str_out << IM->yaxis.Get(i)<<TAB;
  str_out<<'\n';
  for (int i=0;i<IM->Nx;++i)
    {
      for (int j=0;j<IM->Ny;++j)
	str_out << IM->values.get(i+1,j+1)<<TAB;
      str_out << '\n';
    }
}
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
}



