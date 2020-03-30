/*
**********************************************************
IMPACT Vector class

Version 1.4
1
AGRT

22/2/07 moved here
27/2 - Need to write parallel version of class....
7/3 - chkzero member added
3/5/07 - altered so that error checking can easily be 
switched off
**********************************************************

*/


class IMPACT_Vector{		//Vector class (arbitrary length)
 public:		      
  IMPACT_Vector(int array_size); //array_size - -self explanatory	       
  IMPACT_Vector();			
  virtual ~IMPACT_Vector();		 
  void ChangeSize(int array_size);
  int length();
  void Set(int element, double value);
  void SetAll(double value);
  void MultiplyAll(double value);
  double Get(int element);
  void checkerr(int *element);
  void PrintArray();
  void PrintArraytoFile(std::ofstream &outfile);
  void PrintNonZero();
  void reset();
  void unitary();
  void Duplicate(IMPACT_Vector *vec);
  int chkzero(int element);
  double VecSum();
  double **GetVec();
  double * Get_add(int element);
 private:
  int size;
  double *IVec;			
  		
  
};
inline IMPACT_Vector::IMPACT_Vector(int array_size){
  //(adapted from a program off the internet)
  IVec = new double[array_size];
  size = array_size;
  reset();		
}
inline IMPACT_Vector::IMPACT_Vector(){
  //(adapted from a program off the internet)
  IVec = new double[1];				
  size = 1;
  IVec[0]=0.0;		
}
inline IMPACT_Vector::~IMPACT_Vector(){  //This deallocates memory resources		
  delete[] IVec;            //When the object goes out of scope
}
//__________________________________________________
// Below here are access methods
/*The next part changes the size of the vecotr
it maintains elements which overlap with the 
new vector

25/2 - now it doesn't maintain elements as there is no need
*/
inline void IMPACT_Vector::ChangeSize(int array_size){		
  // double* array_temp;
  // array_temp = new double[size];
  // for(int i=0; i<size; ++i) 
  // array_temp[i]=IVec[i];
  delete[] IVec;		
  IVec = new double[array_size];		
  // for(int i=0; i<size; ++i)			
  //     IVec[i] = array_temp[i];	
  size = array_size;
  // delete[] array_temp;		
}
//This function sets the value of an element:
//element is the element number (1 to N);
// value is the value of that elemnt.
inline void IMPACT_Vector::Set(int element, double value)
{
  checkerr(&element);  //check within range
  IVec[element-1]=value;    //set element
}
inline void IMPACT_Vector::SetAll(double value)
{
  for (int element=0;element<size;++element)
    IVec[element]=value;    //set element
}
inline void IMPACT_Vector::MultiplyAll(double value)
{
  for (int element=0;element<size;++element)
    IVec[element]*=value;    //multiply element
}
//This function gets the value of an element:
//element is the element number (1 to N);
//returns value of that elemnt
inline double IMPACT_Vector::Get(int element)
{
  checkerr(&element);
  return IVec[element-1];
}
inline void IMPACT_Vector::PrintArray(){		
  for(int i=0; i<size; ++i)		  
    std::cout<<"("<<i+1<<","<<IVec[i]<<"), ";		
  std::cout<<std::endl<<std::endl;		
}
inline void IMPACT_Vector::PrintArraytoFile(std::ofstream &outfile){    
  for(int i=0; i<size; ++i)     
    outfile<<IVec[i] <<" ";   
  //outfile<<std::endl<<std::endl;    
}
inline void IMPACT_Vector::PrintNonZero(){		
    for(int i=0; i<size; ++i)		  
      if (fabs(IVec[i])>zerotolerance::zerothreshold){std::cout<<"("<<i+1<<","<<IVec[i]<<"), ";}		
  std::cout<<std::endl<<std::endl;		
  }
inline void IMPACT_Vector::Duplicate(IMPACT_Vector *vec)
{
  if (size!=vec->length()){
  size=vec->length();
  delete[] IVec;		
  IVec = new double[size];}
  for (int i=0;i<size;++i)
    IVec[i]=vec->Get(i+1);
}
inline int IMPACT_Vector::length()
{
  return size; 
}
inline void IMPACT_Vector::reset() // resets vector so all elements are zero
{
  for (int i=0;i<size;++i) IVec[i]=0.0;
}
inline void IMPACT_Vector::unitary() // resets vector so all elements are one
{
  for (int i=0;i<size;++i) IVec[i]=1.0+(double)i/(double)size;
}
inline void IMPACT_Vector::checkerr(int *element)
  {
    Vec_check_err(element,&size);
  }

inline int IMPACT_Vector::chkzero(int element) // checks if a particular element is 0
{
  int answer=0;
  if (fabs(Get(element))>zerotolerance::zerothreshold)
    answer=1;
  return answer;
}
inline double IMPACT_Vector::VecSum() // tells me that the vector is zero
{
  double answer=0.0;
  for (int i=0;i<size;++i)
    answer+=IVec[i]*IVec[i];
  return answer;
}
inline double ** IMPACT_Vector::GetVec() // returns entire vector
{
  return &IVec;
}
inline double * IMPACT_Vector::Get_add(int element)
{
  checkerr(&element);
  return &IVec[element-1];
}
