/*
**********************************************************
Functions that make IMPACT_Input work

Unfortunately its very messy programming at present which
should be improved at some point - but probably won't be

Version 1.1
AGRT

6/4/07
20/5/07 - this version contains code to make proper 2D grids
from a file, complete overhaul
**********************************************************
*/
//Turns an input string into a one if input is true or zero otherwise

void IMPACT_dofunction2D(IMPACT_Matrix *gridtemp,int Nx,int Ny,int result)
{
  for (int i=1;i<=Nx;++i)
    for (int j=1;j<=Ny;++j)
      gridtemp->set(i,j,IMPACT_function(gridtemp->get(i,j),result));
}
void plusormultiplygrid2D(char type, IMPACT_Matrix* grid, int Nx,int Ny, double value)
{
  switch (type)
    {
    case '+':
      for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      grid->set(i,j,grid->get(i,j)+value);
      break;
    case '-':
     for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	grid->set(i,j,grid->get(i,j)-value);
      break;
    case '*':
    for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      grid->set(i,j,grid->get(i,j)*value);
      break;
    case '/':
      for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      grid->set(i,j,grid->get(i,j)/value);
      break;
    }
}
void plusormultiplygrids2D(char type, IMPACT_Matrix* grid, 
			   IMPACT_Matrix* othergrid,int Nx, int Ny)
{
  switch (type)
    {
    case '+':
      for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      grid->set(i,j,grid->get(i,j)+othergrid->get(i,j));
      break;
    case '-':
     for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	grid->set(i,j,grid->get(i,j)-othergrid->get(i,j));
      break;
    case '*':
    for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      grid->set(i,j,grid->get(i,j)*othergrid->get(i,j));
      break;
    case '/':
      for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      grid->set(i,j,grid->get(i,j)/othergrid->get(i,j));
      break;
    }
}


/* The next function turns a string within brackets 
   into a maths function on a grid*/

void IMPACT_brktomth2D(char* function,IMPACT_Matrix* grid, double* gridrefx,
		       double* gridrefy,int Nx,int Ny)
{
  int length = strlen(function);
  int count=0,countinner=0;
  char func_para[10]="";
  int ifx=0,ify=0;
  double multiplier, power; //for the functions
  std::istringstream numtemp;
  char tempsignstore='+';
  char powchar[10];

  IMPACT_Matrix gridtemp(Nx,Ny);
  // for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;}//.reset grids
  grid->reset();
  gridtemp.reset();
  while(count<length)
    {
      ifx=ify=0;
      if (function[count]=='-') {tempsignstore='-';++count;}
	      while (!(function[count]=='*'||function[count]=='/'||
		     function[count]=='+'||function[count]=='-'||
		       function[count]==')')&&count<length)
		{
		  func_para[countinner]=function[count];
		  if (toupper(function[count])=='X'
		      &&toupper(function[count-1])!='A') ifx=1;
		  if (toupper(function[count])=='Y') ify=1;
		  count++;			
		  countinner++;
		}
	      if (!ifx&&!ify)
		{
		  multiplier=IMPACT_isConst2D(func_para,gridrefx[Nx-1]-gridrefx[0]
					      ,gridrefy[Ny-1]-gridrefy[0],Nx,Ny);
		  plusormultiplygrid2D(tempsignstore,&gridtemp,Nx,Ny,multiplier);
		}
	      //******************
	      // if x
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
		  IMPACT_Matrix powergrid(Nx,Ny);
		  for (int i=1;i<=Nx;++i)
		    for (int j=1;j<=Ny;++j)
		    {
		      powergrid.set(i,j,pow(gridrefx[i-1],power));
		    }
		  plusormultiplygrids2D(tempsignstore,&gridtemp, &powergrid,Nx,Ny);
		}
	      // NOW y____________________________
	      if (ify)
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
		  IMPACT_Matrix powergrid(Nx,Ny);
		  for (int i=1;i<=Nx;++i)
		    for (int j=1;j<=Ny;++j)
		    {
		      powergrid.set(i,j,pow(gridrefy[j-1],power));
		    }
		  plusormultiplygrids2D(tempsignstore,&gridtemp, &powergrid,Nx,Ny);
		}
              //******************

	      tempsignstore=function[count];
	      if (tempsignstore==')') --count;
	      ++count;
	      countinner=0;
	      if (tempsignstore=='+'||tempsignstore=='-'||count>=length)
		{
		  for (int i=1;i<=Nx;++i)
		    for (int j=1;j<=Ny;++j)
		    {
		      grid->set(i,j,grid->get(i,j)+gridtemp.get(i,j));
		      gridtemp.set(i,j,0.0);
		    }
		}
	     
	      for (int i=0;i<10;++i)
		func_para[i]=0;
    }
}

void IMPACT_brastomth2D(char* function,IMPACT_Matrix * grid, double* gridrefx,
			double* gridrefy,int Nx,int Ny)
{
  int length = strlen(function);
  int count = 0,countinit=0,countinit2=0;
  //double gridtemp[N];
  IMPACT_Matrix gridtemp(Nx,Ny);
  char funcinbra[200]="";
  char tempsignstore='+',signstore='+';
  std::istringstream numtemp;
  double power;
  char numt[10];
  grid->reset();
  gridtemp.reset();
  // for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;}//.reset grids
  while (count<length)
    {
      countinit=count;
      while (function[count]!='('&&function[count]!=')'&&count<length)
	{
	  funcinbra[count-countinit] = function[count];
	  count++;
	}
      IMPACT_brktomth2D(funcinbra,&gridtemp,gridrefx,gridrefy,Nx,Ny);
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
	      for (int i=1;i<=Nx;++i)
		for (int j=1;j<=Ny;++j)
		  gridtemp.set(i,j,pow(gridtemp.get(i,j),power));
	    }
	  tempsignstore=function[count];
	  if (function[count+1]=='(') count++;
	}

      plusormultiplygrids2D(signstore,grid, &gridtemp,Nx,Ny);
      count++;
      for (int i=0;i<200;++i)
	funcinbra[i]=0;
    }
  
}
void IMPACT_strtomth2D(char* function,IMPACT_Matrix* grid, double* gridrefx,
		       double * gridrefy, int Nx,int Ny);
void IMPACT_fnctomth2D(char* function, IMPACT_Matrix* grid, double* gridrefx,
		       double * gridrefy, int Nx,int Ny)
{
  int length = strlen(function);
  int count = 0,countinit=0,countfunc2=0;
  IMPACT_Matrix gridtemp(Nx,Ny),anothergrid(Nx,Ny); 
  char func[6]="";
  char functemp[200]="";
  char funckytemp[100]="";
  char sign='+',tempsign='+',anothersign='*';
  int nbrackets=0,result=0;
  grid->reset(); 
  gridtemp.reset();
  anothergrid.reset();
//for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;anothergrid[i]=0.0;}//.reset grids
  // first get the number of functions...
  countinit=count;
  while (count<length)
    {
      functemp[count-countinit]=function[count];
      func[0]=func[1];
      func[1]=func[2];
      func[2]=func[3];
      func[3]=func[4];
      func[4]=function[count];
      func[5]=0;
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
	  IMPACT_brastomth2D(functemp,&gridtemp,gridrefx,gridrefy,Nx,Ny); 
	  plusormultiplygrids2D(sign,grid, &gridtemp,Nx,Ny);
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
	  for (int i=0;i<(int)strlen(functemp);++i)
		{
		  func[0]=func[1];
		  func[1]=func[2];
		  func[2]=func[3];
		  func[3]=func[4];
		  func[4]=function[i];
		  func[5]=0;
		  result = IMPACT_isFunc(func);
		  if (result>0) break;
		}

	  if (result>0)IMPACT_strtomth2D(functemp,&gridtemp,gridrefx, 
					 gridrefy,Nx,Ny);
	  else IMPACT_brastomth2D(functemp,&gridtemp,gridrefx,gridrefy,Nx,Ny); 
	  IMPACT_dofunction2D(&gridtemp,Nx,Ny,result);
	  if (funckytemp[0]!=0) 
	    {
	      IMPACT_brastomth2D(funckytemp,&anothergrid,gridrefx,
				 gridrefy,Nx,Ny); 
	      plusormultiplygrids2D(anothersign,&anothergrid, &gridtemp,Nx,Ny);
	      for (int i=1;i<=Nx;++i)
		for (int j=1;j<=Ny;++j)
		  gridtemp.set(i,j,anothergrid.get(i,j));
	    }
	  plusormultiplygrids2D(sign,grid, &gridtemp,Nx,Ny);
	  sign = function[count];
	  result=0;
	  for (int i=0;i<200;++i)
	    functemp[i]=0;
	  countinit=count+1;
	}
      ++count;
    }

  IMPACT_brastomth2D(functemp,&gridtemp,gridrefx,gridrefy,Nx,Ny); 
  plusormultiplygrids2D(sign,grid, &gridtemp,Nx,Ny);
}
/*
  The next function creates a gridpoint array for putting in the config file
  from the information given in the 
 */
void IMPACT_strtomth2D(char* function,IMPACT_Matrix* grid, double* gridrefx,
		       double * gridrefy, int Nx,int Ny)
{
  int length = strlen(function);
  int nbrackets = 0,count=0,countinit=0,countinit2=0;
  char functemp[200]=""; char sign='+';
  char numt[10];
  IMPACT_Matrix gridtemp(Nx,Ny);
  double power=1; int anothertag=0;
  int tag=IMPACT_Messages::if_show_function_input;
  std::istringstream numtemp;
  grid->reset();
  gridtemp.reset();
  // for (int i=0;i<N;++i) {grid[i]=0.0; gridtemp[i]=0.0;}//.reset grids
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
					 function[count-1]=='*'))) anothertag--;
	  if (function[count]=='(') anothertag++;
	  functemp[count-countinit]=function[count];
	  count++;
	}
      if (functemp[0]!=0&&tag) std::cout<<sign;
      if (functemp[0]!=0&&tag) std::cout<<functemp<<'\n';
      IMPACT_fnctomth2D(functemp,&gridtemp,gridrefx,gridrefy,Nx,Ny);
      plusormultiplygrids2D(sign,grid, &gridtemp,Nx,Ny);
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
      IMPACT_fnctomth2D(functemp,&gridtemp,gridrefx,gridrefy,Nx,Ny);
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
	  for (int i=1;i<=Nx;++i)
	    for (int j=1;j<=Ny;++j)
	      gridtemp.set(i,j,pow(gridtemp.get(i,j),power));
	  if (functemp[0]!=0&&tag) std::cout<<"^"<<power;
	}   
      // if (tag) std::cout<<'\n';
      plusormultiplygrids2D(sign,grid, &gridtemp,Nx,Ny);
      if(function[count]=='+'||function[count]=='-'||function[count]=='*'||
	 function[count]=='/') sign=function[count]; 
      for (int i=0;i<200;++i)
	functemp[i]=0;
      ++count;
    }
}

int IMPACT_Make_Function2D(std::ifstream &infile, std::string grid_temp,IMPACT_IMatrix *Initializer, double* gridrefx, double* gridrefy,int Nx,int Ny)
{
  int ok_message = 0; //a return to show nothing went wrong.
  char type[4] = {0,0,0,0};
  int tag=IMPACT_Messages::if_show_function_input; //if echo function input
  
  for (int i=0;i<3;++i)
    type[i] = (char)toupper(grid_temp.c_str()[i]);

  if (!strcmp(type,"CST")) // a constant value is given
    {
      ok_message=1;
      Initializer->SwitchOnConst();
      double cst_value=0.0;
      infile >>cst_value;
      /*for (int i=0;i<N;++i)
	grid[i]=cst_value;*/
      Initializer->setc(cst_value);
      if (tag) std::cout << "Const. val = "<<cst_value<<'\n';
    }
  /*if (!strcmp(type,"GRD")) // a grid of values is specified
    {
      ok_message=1;
      double grd_value=0.0;
      for (int i=0;i<N;++i)
	{
	  infile >>grd_value;
	    grid[i]=grd_value;
	}
      if (tag) std::cout << "User specified grid\n";
      }*/
  if (!strcmp(type,"FLE")) // get grid from text file
    {
      ok_message=1;
      Initializer->SwitchOffConst();
      Initializer->ChangeSize(Nx,Ny);
      std::string filename;
      infile>>filename;
      filename = IMPACT_Messages::Input_Directory + filename;
      
      std::ifstream fin(filename.c_str());
      int ipos=1,jpos=1;
      double grd_value;
      while (fin>>grd_value)
	{
	  Initializer->set(ipos,jpos, grd_value);
	  ++ipos;
	  if (ipos>Nx) {ipos=1;++jpos;}
	  if (jpos>Ny) break;
	}
      if (jpos<Ny) 
	{
	  std::cout<<"IMPACT: ERROR - In specified grid file, "<<filename<<": Not enough points or file not found.\n";
	  exit(0);
	}
      fin.close();
      if (tag) std::cout << "User specified file - "<<filename<<"\n";
    }
  
  if (!strcmp(type,"FNC")) // get grid from function
    {
      ok_message=1;
      Initializer->SwitchOffConst();
      Initializer->ChangeSize(Nx,Ny);
      /*if (!(Nx-1))
	{
	  std::cout<< "IMPACT: ERROR - in input deck, function specified for single gridpoint\n\n";
	  exit(0);
	  }*/
      char function[200];
      infile.getline(function,200,'{');
      infile.getline(function,200,'}');
      double *refgridsmallx;
      refgridsmallx = new double[Nx];
      for (int i=0;i<Nx;++i)
	refgridsmallx[i]=gridrefx[i+1];
      double *refgridsmally;
      refgridsmally = new double[Nx];
      for (int i=0;i<Ny;++i)
	refgridsmally[i]=gridrefy[i+1];
      //now we need to break the function into a series of strings
      //IMPACT_strtomth2D(function,Initializer,refgridsmallx,refgridsmally,Nx,Ny);
      double **gridtemp;
      gridtemp =new double*[Nx];
      for (int i=0;i<Nx;++i)
	{
	  gridtemp[i]=new double[Ny];
	  //gridtemp[i][0]=0.0;
	}
      IMPACTA_Get_Func(function,gridtemp,Nx,Ny,refgridsmallx,refgridsmally);
      for (int i=0;i<Nx;++i)
	for (int j=0;j<Ny;++j)
	  Initializer->set(i+1,j+1,gridtemp[i][j]);
      for (int i=0;i<Nx;++i) {
	delete[] gridtemp[i];
      }
      delete[] gridtemp;
      delete[] refgridsmallx;
      delete[] refgridsmally;
    }
  return ok_message;
}
