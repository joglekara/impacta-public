/*
**********************************************************
Functions that make IMPACT_Input work

Unfortunately its very messy programming at present which
should be improved at some point - but probably won't be

Version 1.1
AGRT

6/4/07
10/4/08 - Had to change make axis function as my assumption
that only periodic bounds used gridpoints out of the domain 
was incorrect.
**********************************************************
*/
//Turns an input string into a one if input is true or zero otherwise
int IMPACT_truth(std::string truth)
{
  int result = 0;
  char truthchar[5]={0,0,0,0,0};
  for (int i=0;i<4;++i)
    truthchar[i]=toupper(truth.c_str()[i]);
  if(!strcmp(truthchar,"TRUE")) result=1;
  return result;
}
void IMPACTA_ZapNan(double *grid,int N)
{
  for (int i=0;i<N;++i)
    if (!(grid[i]==grid[i])) grid[i]=0.0;
}
void plusormultiplygrid(char type, double* grid, int N, double value)
{
  switch (type)
    {
    case '+':
      for (int i=0;i<N;++i)
	grid[i]+=value;
      break;
    case '-':
      for (int i=0;i<N;++i)
	grid[i]-=value;
      break;
    case '*':
      for (int i=0;i<N;++i)
	grid[i]*=value;
      break;
    case '/':
      for (int i=0;i<N;++i)
	grid[i]/=value;
      break;
    }
}
void plusormultiplygrids(char type, double* grid, int N, double* othergrid)
{
  switch (type)
    {
    case '+':
      for (int i=0;i<N;++i)
	grid[i]+=othergrid[i];
      break;
    case '-':
      for (int i=0;i<N;++i)
	grid[i]-=othergrid[i];
      break;
    case '*':
      for (int i=0;i<N;++i)
	grid[i]*=othergrid[i];
      break;
    case '/':
      for (int i=0;i<N;++i)
	grid[i]/=othergrid[i];
      break;
    }
}
void IMPACT_InputDeckError(std::string errmes)
{
  std::cout << "IMPACT: ERROR - In Input Deck, "<<errmes<<"\n";
  exit(0);
}
std::string onoroff(int on)
{
  std::string result = "Off";
  if (on) result = "On";
  return result;
}

/* The next function turns a string within brackets 
   into a maths function on a grid*/
void IMPACT_brktomth(char* function,double* grid, int N,double* gridref)
{
  int length = strlen(function);
  int count=0,countinner=0;
  char func_para[10]="";
  int ifx=0;
  double multiplier, power; //for the functions
  std::istringstream numtemp;
  char tempsignstore='+';
  char powchar[10];

  double * gridtemp;
  gridtemp = new double[N];
  for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;}//.reset grids
  
  while(count<length)
    {
      ifx=0;
      if (function[count]=='-') {tempsignstore='-';++count;}
	      while (!(function[count]=='*'||function[count]=='/'||
		     function[count]=='+'||function[count]=='-'||
		       function[count]==')')&&count<length)
		{
		  func_para[countinner]=function[count];
		  if (toupper(function[count])=='X'
		      &&toupper(function[count-1])!='A') ifx=1;
		  count++;			
		  countinner++;
		}
	      if (!ifx)
		{
		  multiplier=IMPACT_isConst(func_para,gridref[N-1]-gridref[0],N);
		  plusormultiplygrid(tempsignstore,gridtemp,N,multiplier);
		}
	      if (ifx)
		{
		  power =1;
		  for (int i=0;i<countinner;++i)
		    if (func_para[i]=='^')
		      {
			int j=i+1;
			for (int k=0;k<10;k++)
			  powchar[k]=0;
			while (func_para[j]!='+'&&func_para[j]!='-'&&
			       func_para[j]!='*'&&func_para[j]!='/'
			       &&func_para[j]!=')'&&j<length)
			  {
			    powchar[j-i-1]=func_para[j];
			    ++j;
			  }
			numtemp.str(powchar);
			numtemp>>power;
			break;
			}
		  double *powergrid;
		  powergrid = new double[N];
		  for (int i=0;i<N;++i)
		    {
		      powergrid[i]=pow(gridref[i],power);
		    }
		  plusormultiplygrids(tempsignstore,gridtemp,N, powergrid);
		  delete[] powergrid;
		}
	      tempsignstore=function[count];
	      if (tempsignstore==')') --count;
	      ++count;
	      countinner=0;
	      if (tempsignstore=='+'||tempsignstore=='-'||count>=length)
		{
		  for (int i=0;i<N;++i) //unload gridtemp
		    {
		      grid[i]+=gridtemp[i];
		      gridtemp[i]=0.0;
		    }
		}
	     
	      for (int i=0;i<10;++i)
		func_para[i]=0;
    }
  delete[] gridtemp;
}
void IMPACT_brastomth(char* function,double* grid, int N,double* gridref)
{
  int length = strlen(function);
  int count = 0,countinit=0,countinit2=0;
  double *gridtemp;
  gridtemp = new double[N];
  char funcinbra[200]="";
  char tempsignstore='+',signstore='+';
  std::istringstream numtemp;
  double power;
  char numt[10];
  for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;}//.reset grids
  while (count<length)
    {
      countinit=count;
      while (function[count]!='('&&function[count]!=')'&&count<length)
	{
	  funcinbra[count-countinit] = function[count];
	  count++;
	}
      IMPACT_brktomth(funcinbra,gridtemp,N,gridref);
      signstore=tempsignstore;
      if(function[count]=='+'||function[count]=='-'||function[count]=='*'||
	 function[count]=='/') tempsignstore=function[count];
      if (function[count]=='('&&count>0) tempsignstore=function[count-1];
      if (function[count]==')') 
	{
	  count++;
	   if (function[count]=='^')
	    {
	      ++count;
	      countinit2=count;
	      for(int i=0;i<10;++i)
		numt[i]=0;
	      while((function[count]!='+'&&function[count]!='-'&&
		    function[count]!='*'&&function[count]!='/'
		     &&function[count]!=')')&&count<length)
		{
		  numt[count-countinit2]=function[count];
		  ++count;
		}
	      numtemp.str(numt);
	      numtemp>>power;
	      for (int i=0;i<N;++i)
		gridtemp[i]=pow(gridtemp[i],power);
	    }
	  tempsignstore=function[count];
	  if (function[count+1]=='(') count++;
	}

      plusormultiplygrids(signstore,grid,N, gridtemp);
      count++;
      for (int i=0;i<200;++i)
	funcinbra[i]=0;
    }
  delete[] gridtemp;
}

void IMPACT_fnctomth(char* function,double* grid, int N,double* gridref)
{
  int length = strlen(function);
  int count = 0,countinit=0,countfunc2=0;
  double *gridtemp;
  gridtemp=new double[N];
  double *anothergrid;
  anothergrid=new double[N]; 
  
  char func[5]="";
  char functemp[200]="";
  char funckytemp[100]="";
  char sign='+',tempsign='+',anothersign='*';
  int nbrackets=0,result=0;
  for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;anothergrid[i]=0.0;}//.reset grids
  // first get the number of functions...
  countinit=count;
  while (count<length)
    {
	  functemp[count-countinit]=function[count];
	  func[0]=func[1];
	  func[1]=func[2];
	  func[2]=func[3];
	  func[3]=function[count];
	  func[4]=0;
	  result = IMPACT_isFunc(func);
	  if(result>0)
	    {
	      countfunc2=count;
	      while (function[countfunc2]!='+'&&function[countfunc2]!='-'&&
		     function[countfunc2]!='*'&&function[countfunc2]!='/'
		     &&countfunc2>0)
		{
		  countfunc2--; 
		}
	      if (function[countfunc2]=='*'||function[countfunc2]=='/')
		{
		  anothersign=function[countfunc2];
		  if (function[countfunc2-1]==')')
		    {
		      do
			{
			  countfunc2--;
			  for (int j=99;j>0;--j)
			    funckytemp[j]=funckytemp[j-1];
			  funckytemp[0]=function[countfunc2]; 
			}
			while(function[countfunc2]!='('&&countfunc2>0);
		    }
		  else
		    {
		      do
			{
			  countfunc2--;
			  for (int j=99;j>0;--j)
			    funckytemp[j]=funckytemp[j-1];
			  funckytemp[0]=function[countfunc2]; 
			}
			while(function[countfunc2]!='+'&&function[countfunc2]!='-'
			      &&countfunc2>0);
		    } 
		}
	      for (int i=0;i<200;++i)
		functemp[i]=0;
	      if (countfunc2-countinit>0)
		for (int i=countinit;i<countfunc2;++i)
		  functemp[i]=function[i];

	      IMPACT_brastomth(functemp,gridtemp,N,gridref); 
	      plusormultiplygrids(sign,grid,N, gridtemp);
	      ++count;
	      if (function[count]!='(')
		{
		  std::cout<<"IMPACT: ERROR - in function in input deck\n";
		  std::cout<<function[count];
		  exit(0);
		}
	      tempsign = function[count-1];
	      if (!(tempsign=='+'||tempsign=='-'||tempsign=='/'
		  ||tempsign=='*')) tempsign=sign;
	      nbrackets=1;
	      ++count;
	      for (int i=0;i<200;++i)
		functemp[i]=0;
	      countinit=count;
	      while (nbrackets>0&&nbrackets<3&&count<length)
		{
		  if (function[count]=='(') nbrackets++;
		  if (function[count]==')') nbrackets--;
		  if(nbrackets>0) functemp[count-countinit]=function[count];
		  ++count;
		}

	      IMPACT_brastomth(functemp,gridtemp,N,gridref);
	      IMPACT_dofunction(gridtemp,N,result);
	      if (funckytemp[0]!=0) 
		{
		  IMPACT_brastomth(funckytemp,anothergrid,N,gridref);
		  plusormultiplygrids(anothersign,anothergrid,N, gridtemp);
		  for(int i=0;i<N;++i)
		    gridtemp[i]=anothergrid[i];
		}
	      plusormultiplygrids(tempsign,grid,N, gridtemp);
	      sign = function[count];
	      result=0;
	      for (int i=0;i<200;++i)
		functemp[i]=0;
	      countinit=count+1;
	    }
	  ++count;
	}
  IMPACT_brastomth(functemp,gridtemp,N,gridref);
  plusormultiplygrids(sign,grid,N, gridtemp);
  delete[] gridtemp;
  delete[] anothergrid;
}
/*
  The next function creates a gridpoint array for putting in the config file
  from the information given in the 
 */
void IMPACT_strtomth(char* function,double* grid, int N,double* gridref)
{
  int length = strlen(function);
  int nbrackets = 0,count=0,countinit=0,countinit2=0;
  char functemp[200]=""; char sign='+';
  char numt[10];
  double *gridtemp;
  gridtemp=new double[N];
  double power=1; int anothertag=0;
  int tag=IMPACT_Messages::if_show_function_input;
  std::istringstream numtemp;
  for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;}//.reset grids
  while (count<length)
    {
      if(function[count]=='-') {sign = '-';++count;}
      countinit=count;
      anothertag=0;
      while (anothertag<1&&count<length)
	{
	  if (count>0)
	    if (!(function[count]=='('&&(function[count-1]=='+'||
		  function[count-1]=='-'|| function[count-1]=='/'||
		  function[count-1]=='*'))) { anothertag--; }
	  if (function[count]=='(') anothertag++;
	  functemp[count-countinit]=function[count];
	  count++;
	}
      if (functemp[0]!=0&&tag) std::cout<<sign;
      if (functemp[0]!=0&&tag) std::cout<<functemp<<'\n';
      IMPACT_fnctomth(functemp,gridtemp,N,gridref);
      plusormultiplygrids(sign,grid,N, gridtemp);
      if(function[count]=='+'||function[count]=='-'||function[count]=='*'||
	 function[count]=='/') sign=function[count]; 
      for (int i=0;i<200;++i)
	functemp[i]=0;
      
      nbrackets=0;
      if (function[count]=='(') {nbrackets++;count++;}
      countinit=count;
      while (nbrackets>0&&nbrackets<4&&count<length)
	{
	  if (function[count]=='(') nbrackets++;
	  if (function[count]==')') nbrackets--;
	  if(nbrackets>0) functemp[count-countinit]=function[count];
	  ++count;
	}
      if (functemp[0]!=0&&tag) std::cout<<sign;
      if (functemp[0]!=0&&tag) std::cout<<'('<<functemp<<')';
      IMPACT_fnctomth(functemp,gridtemp,N,gridref);
      if (function[count]=='^')
	{
	  ++count;
	  countinit2=count;
	  for(int i=0;i<10;++i)
	    numt[i]=0;
	  while((function[count]!='+'&&function[count]!='-'&&
		 function[count]!='*'&&function[count]!='/')&&count<length)
	    {
	      numt[count-countinit2]=function[count];
	      ++count;
	    }
	  numtemp.str(numt);
	  numtemp>>power;
	  for (int i=0;i<N;++i)
	    gridtemp[i]=pow(gridtemp[i],power);
	  if (functemp[0]!=0&&tag) std::cout<<"^"<<power;
	}   
      // if (tag) std::cout<<'\n';
      plusormultiplygrids(sign,grid,N, gridtemp);
      if(function[count]=='+'||function[count]=='-'||function[count]=='*'||
	 function[count]=='/') sign=function[count]; 
      for (int i=0;i<200;++i)
	functemp[i]=0;
      ++count;
    }
  delete[] gridtemp;
}

int IMPACT_Make_Function(std::ifstream &infile, std::string grid_temp,double* grid, double* gridref, int N,double min)
{
  int ok_message = 0; //a return to show nothing went wrong.
  char type[4] = {0,0,0,0};
  int tag=IMPACT_Messages::if_show_function_input; //if echo function input
  for (int i=0;i<3;++i)
    type[i] = (char)toupper(grid_temp.c_str()[i]);

  if (!strcmp(type,"MAX")) // a maximum value is given
    {
      ok_message=1;
      double max_val=0.0;
      infile >>max_val;
      double cst_value=fabs(max_val-min)/(N);
      /*if (!(N))
	{
	  std::cout<< "IMPACT: ERROR - in input deck, Max value specified for one gridpoint\n\n";
	  exit(0);
	  }*/
      for (int i=0;i<N;++i)
	grid[i]=cst_value;
      if (tag) std::cout << "Const. val = "<<cst_value<<'\n';
    }
  if (!strcmp(type,"CST")) // a constant value is given
    {
      ok_message=1;
      double cst_value=0.0;
      infile >>cst_value;
      for (int i=0;i<N;++i)
	grid[i]=cst_value;
      IMPACTA_ZapNan(grid,N);
      if (tag) std::cout << "Const. val = "<<cst_value<<'\n';
    }
  if (!strcmp(type,"GRD")) // a grid of values is specified
    {
      ok_message=1;
      double grd_value=0.0;
      for (int i=0;i<N;++i)
	{
	  infile >>grd_value;
	    grid[i]=grd_value;
	}
      IMPACTA_ZapNan(grid,N);
      if (tag) std::cout << "User specified grid\n";
    }
  if (!strcmp(type,"FLE")) // get grid from text file
    {
      ok_message=1;
      std::string filename;
      infile>>filename;
      filename = IMPACT_Messages::Input_Directory + filename;
      //std::cout<<filename<<'\n';
      std::ifstream fin(filename.c_str());
      int poscheck=0;
      double grd_value;
      while (fin>>grd_value)
	{
	  grid[poscheck] = grd_value;
	  ++poscheck;
	  if (poscheck>N-1) break;
	}
      if (poscheck<N-1) 
	{
	  std::cout<<"IMPACT: ERROR - In specified grid file, not enough points or file '"<<filename.c_str()<<"' does not exist.";
	  exit(0);
	}
      fin.close();
      if (tag) std::cout << "User specified file - "<<filename<<"\n";
      IMPACTA_ZapNan(grid,N);
    }
  if (!strcmp(type,"FNC")) // get grid from function
    {
      ok_message=1;
      if (!(N-1))
	{
	  std::cout<< "IMPACT: ERROR - in input deck, function specified for single gridpoint.\n\n";
	  exit(0);
	}
      char function[200];
      infile.getline(function,200,'{');
      infile.getline(function,200,'}');
      //now we need to break the function into a series of strings
      //IMPACT_strtomth(function,grid,N,gridref);
      // NEW NEW - Use IMPACTA_Get_Func
      double *gridy;
      double **gridtemp;
      gridtemp =new double*[N];
      for (int i=0;i<N;++i)
	{
	  gridtemp[i]=new double[1];
	  gridtemp[i][0]=0.0;
	}
      gridy=new double[1];
      gridy[0]=1.0;
      IMPACTA_Get_Func(function,gridtemp,N,1,gridref,gridy);
      for (int i=0;i<N;++i)
	grid[i]=gridtemp[i][0];
      delete[] gridy;
      for (int i=0;i<N;++i)
	delete[] gridtemp[i];
      delete[] gridtemp;
      IMPACTA_ZapNan(grid,N);
    }
  if (!ok_message) 
    {
      std::cout<< "IMPACT: ERROR - in initial condition in input deck.\n\n";
      exit(0);
    }
  IMPACTA_ZapNan(grid,N);
  return ok_message;
}
int IMPACT_Make_Axis(std::ifstream &infile, std::string grid_temp,double* grid, int N,double min,int direction) //direction added to get bounds correct
{
  int result;
  double *gridref;
  gridref=new double[N];
  double *gridsmall;
  gridsmall=new double[N];

  for (int i=0;i<N;++i)
    gridref[i]=1.0+i;
  result = IMPACT_Make_Function(infile, grid_temp,gridsmall, gridref, N,min);
  
 
  for (int i=0;i<N;++i)
    grid[i+1]=gridsmall[i];
  // finally to deal with periodic boundary conditions....
  // noting that non periodic boundaries don't use points outside the boundary
  // the above is incorrect - now adjusted to correct for different bounds
  if (direction==3||!strcmp(boundaries[direction-1],"periodic")) 
    {
      grid[0]=grid[N];
      grid[N+1]=grid[1];
    }
  else
    {
      grid[0]=grid[1];
      grid[N+1]=grid[N];
    }
  IMPACTA_ZapNan(grid,N);
  delete[] gridref;
  delete[] gridsmall;
  return result;
}
int IMPACT_Make_Grid(std::ifstream &infile, std::string grid_temp,double* grid,double* refgrid,int N)
{
  double *refgridsmall;
  refgridsmall = new double[N];
  for (int i=0;i<N;++i)
    refgridsmall[i]=refgrid[i+1];
  int result;
  result = IMPACT_Make_Function(infile, grid_temp,grid, refgridsmall, N,0.0);
  IMPACTA_ZapNan(grid,N);
  delete[] refgridsmall;
  return result;
}
void IMPACT_Make_2DGrid(IMPACT_IMatrix *answer, double* gridx, double* gridy,int Nx,int Ny)
{
  IMPACTA_ZapNan(gridx,Nx);
  IMPACTA_ZapNan(gridy,Ny);
  double concheck=0.0;
  for (int i=0;i<Nx;++i)
    concheck += fabs(gridx[i]-gridx[0]);
  for (int i=0;i<Ny;++i)
    concheck += fabs(gridy[i]-gridy[0]);
  if (fabs(concheck)<zerotolerance::zerothreshold) 
    {
      answer->setc(gridx[0]*gridy[0]);
    }
  else
  {
      answer->setg(gridx,gridy,Nx,Ny);
  }
}
int IMPACT_Get_Bhat(std::ifstream &infile, std::string grid_temp,double *Bhat)
{
  char type[4] = {0,0,0,0};
  int tag=IMPACT_Messages::if_show_function_input; //if echo function input
  for (int i=0;i<3;++i)
    type[i] = (char)toupper(grid_temp.c_str()[i]);
  std::istringstream numtemp1,numtemp2,numtemp3;
   int result=0;
   char function[200];
   if (!strcmp(type,"BX")) // a constant value is given
     {
       result=1;
       Bhat[0]=1;
     }
   if (!strcmp(type,"BY")) // a constant value is given
     {
       result=1;
       Bhat[1]=1;
     }
   if (!strcmp(type,"BZ")) // a constant value is given
     {
       result=1;
       Bhat[2]=1;
     }
   if (!strcmp(type,"VEC")) // a constant value is given
     {
       result=1;
       infile.getline(function,200,'(');
       infile.getline(function,200,',');
       numtemp1.str(function);
       numtemp1>>Bhat[0];
       infile.getline(function,200,',');
       numtemp2.str(function);
       numtemp2>>Bhat[1];
       infile.getline(function,200,')');
       numtemp3.str(function);
       numtemp3>>Bhat[2];
     }
   double Bmag = sqrt(Bhat[0]*Bhat[0]+Bhat[1]*Bhat[1]+Bhat[2]*Bhat[2]);
   if (Bmag!=0.0)
     for (int i=0;i<3;++i)
       Bhat[i]/=Bmag;
   if (tag) std::cout<< "B hat:\t\t("<<Bhat[0]<<", "<<Bhat[1]<<", "<<Bhat[2]<<")\n";
   return result;
}
void IMPACT_strcpy(char *argstr,std::string str)
{
  int length = strlen(str.c_str());
  for (int i=0;i<length;++i)
      argstr[i]=str.c_str()[i];
  // std::cout<<str.c_str()[i]<<','<<argstr[i]<<'\n';}
}

