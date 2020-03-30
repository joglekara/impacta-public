#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <iostream>
using namespace std;

#define     amu      1.66053873E-27
#define     kB       1.380658E-23
#define     qe       1.60217733E-19

//extern double Te_ext,ne_ext,AN,ZA;

double eos(double, double, double, double, double) ;


int main(void)
{
  double Te=0.0, ne, ZA, AN, Z;
  
  
  cout << "Temperature (eV) = " ;
  cin >> Te ;

  while (Te>=0.0)
    {
      cout << "ne (m^-3) = " ;
      cin >> ne ;
      cout << "ZA = " ;
      cin >> ZA ;
      cout << "AN = " ;
      cin >> AN ;
      cout << "Z = " ;
      cin >> Z ;
      
      cout << eos(Z,Te,ne,ZA,AN) << endl << endl ;

      cout << "Temperature (eV) = " ;
      cin >> Te ;
    }

}



/* EOS 
   Calculation of Zfree using the Thomas-Fermi model (Salzmann, 1998) */

double eos(double Z, double Te, double ne, double ZA, double AN)
{
  //	double Te,ne;
	double rho,R,T0,alpha,beta;
	double TF,a1,a2,a3,a4,A,b0,b1,b2,B,c1,c2,CC,Q1,Q,x,ans;
	
//	Te = Te_ext*kB/qe;
//      ne = ne_ext*1.0E-6;
	
	rho   = ne/Z*AN*amu*1.0E3;
	R     = rho/(ZA*AN);
	T0    = Te/pow(ZA,4.0/3.0);
	alpha = 14.3139;
	beta  = 0.6624;
	TF    = T0/(1.0+T0);
	a1    = 3.323E-3;
	a2    = 0.971832;
	a3    = 9.26148E-5;
	a4    = 3.10165;
	A     = a1*pow(T0,a2)+a3*pow(T0,a4);
	b0    = -1.7630;
	b1    = 1.43175;
	b2    = 0.315463;
	B     = -exp(b0+b1*TF+b2*pow(TF,7.0));
	c1    = -0.366667;
	c2    = 0.983333;
	CC    = c1*TF+c2;
	Q1    = A*pow(R,B);
	Q     = pow(pow(R,CC)+pow(Q1,CC),1.0/CC);
	x     = alpha*pow(Q,beta);
	
	ans = Z-ZA*x/(1.0+x+sqrt(1.0+2.0*x));
	
	return ans;
}
