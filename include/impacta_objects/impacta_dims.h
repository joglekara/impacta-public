/*
**********************************************************
"Dimension" class - basically an integer but allows
unambiguous definitions of Variable class function calls

Version 2.9
AGRT

21/3/07
update 11/2/08
**********************************************************
*/

class IMPACT_Dim //err, this is to make sense of the variable class
{
 private:
  int x1;
 public:
  IMPACT_Dim()
    {x1=1;}
  ~IMPACT_Dim()
    {}
IMPACT_Dim(int dimension)
    {x1=dimension;}

  void operator=(int dimension)
    {
      x1=dimension;
    }
  int get()
  {return x1;}
 void set(int value)
  {x1=value;}
  void operator++()
    {
      ++x1;
    }
  void operator--()
    {
      --x1;
    }
  void operator+=(int dimension)
    {
      x1+=dimension;
    }
  IMPACT_Dim operator+(int dimension)
    {
      return IMPACT_Dim(x1+dimension);
    }
  IMPACT_Dim operator*(int dimension)
    {
      return IMPACT_Dim(x1*dimension);
    }

  IMPACT_Dim operator-(int dimension)
    {
      return IMPACT_Dim(x1-dimension);
    }
 void operator+=(IMPACT_Dim dimension)
    {
      x1+=dimension.get();
    }
  IMPACT_Dim operator+(IMPACT_Dim dimension)
    {
      return IMPACT_Dim(x1+dimension.get());
    }
  
  IMPACT_Dim operator-(IMPACT_Dim dimension)
    {
      return IMPACT_Dim(x1-dimension.get());
    }
IMPACT_Dim operator*(IMPACT_Dim dimension)
    {
      return IMPACT_Dim(x1*dimension.get());
    }
  int operator==(int dimension)
    {
      return x1==dimension;
    }
int operator<(int dimension)
    {
      return x1<dimension;
    }
int operator>(int dimension)
    {
      return x1>dimension;
    }
int operator<=(int dimension)
    {
      return x1<=dimension;
    }
int operator>=(int dimension)
    {
      return x1>=dimension;
    }
 char coord()
 {
   char ans;
   switch (x1)
     {
     case 1:
       ans='x';
       break;
     case 2:
       ans='y';
       break;
     case 3:
       ans='z';
       break;
     default:
       std::cout<<"Not a coordinate x,y,z\n";
	 exit(0);
     }
   return ans;
 }

};
// gets components orthogonal to x1
inline void GetOrthogonal(IMPACT_Dim *x1,IMPACT_Dim *x2,IMPACT_Dim *x3)
{
  *x2=*x1+1;
  if (*x2>3) *x2=IMPACT_Dim(1);
  *x3=*x2+1;
  if (*x3>3) *x3=IMPACT_Dim(1); 
}
void GetOrthogonal(IMPACT_Dim *x1,IMPACT_Dim *x2)
{
  *x2=*x1+1;
  if (*x2>2) *x2=IMPACT_Dim(1);
}
void Set_Diff_Dir(IMPACT_Dim *x1,int *iplus,int *iminus,
			int *jplus,int *jminus)
{
  if (*x1==1) {*iplus=*iplus+1; *iminus=*iminus-1;}
  if (*x1==2) {*jplus=*jplus+1; *jminus=*jminus-1;}
}
	
std::string GetDirChar(int x1)
{
  std::string answer;
  switch (x1)
    {
    case 1:
      answer="x"; break;
    case 2:
      answer="y"; break;
    case 3:
      answer="z"; break;
    }
  return answer;
}
std::string Get2DirChar(int x1)
{
  std::string answer="";
  switch (x1)
    {
    case 1:
      answer="xx"; break; //0: 0 0
    case 2:
      answer="xy"; break; //1: 0 1
    case 3:
      answer="yy"; break; //2: 1 1
    case 4:
      answer="xz"; break;  //3: 0 3
    case 5:
      answer="yz"; break; //4: 1 3
    }
  return answer;
}
char getcoord(int i)
{
  char x1='x';
  switch (i)
    {
    case 1:
      x1='x';break;
    case 2:
      x1='y';break;
    case 3:
      x1='z';break;
    default:
      std::cout<<"IMPACT: ERROR - in getcoord";
    }
  return x1;
}
int getcoord(char x)
{
  int x1=0;
  switch (x)
    {
    case 'x':
      x1 = 1;break;
    case 'y':
      x1 = 2;break;
    case 'z':
      x1 = 3;break;
    case 'c':
      x1 = 4;break; // for circular polarization
    case 'd':
      x1 = 5;break; // for diagonal polarization
    default:
      x1 = 0;break;
      //std::cout<<"IMPACT: ERROR - in getcoord";
    }
  return x1;
}
