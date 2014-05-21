#include <iostream>
#include <math.h>
using namespace std;


double H(const double *xx)
{
  const double M     = xx[0];
  const double K1    = xx[1];
  const double K2    = xx[2];
  const double Phi   = xx[3];
  const double Alpha = xx[4];

  std::cout << "M:     " <<  M     << std::endl;
  std::cout << "K1:    " <<  K1    << std::endl;
  std::cout << "K2:    " <<  K2    << std::endl;
  std::cout << "Phi:   " <<  Phi   << std::endl;
  std::cout << "Alpha: " <<  Alpha << std::endl;

  const double Eq1   = -1/(M*sin(Phi-Alpha));
  const double Eq2   = (4*M_PI*M*M-2*K1-4*K2)*cos(Phi)*sin(Phi);
  const double Eq3   = 4*K2*cos(Phi)*sin(Phi)*sin(Phi)*sin(Phi);

  return Eq1*(Eq2+Eq3);
}

double H1(const double *xx)
{
  const double M     = xx[0];
  const double K1    = xx[1];
  const double K2    = xx[2];
  const double Phi   = xx[3];
  const double Alpha = xx[4];

  const double Eq1   = -1/(M*sin(Phi-Alpha));
  const double Eq2   = (4*M_PI*M*M-2*K1-4*K2)*cos(Phi)*sin(Phi);
  const double Eq3   = 4*K2*cos(Phi)*sin(Phi)*sin(Phi)*sin(Phi);
  const double Eq4   = M*cos(Alpha-Phi);
  const double Eq5   = (4*M_PI*M*M-2*K1-4*K2)*cos(2*Phi);
  const double Eq6   = 4*K2*(3*sin(Phi)*sin(Phi)*(cos(Phi)*cos(Phi)-sin(Phi)*sin(Phi)*sin(Phi)*sin(Phi)));

  return Eq1*(Eq2+Eq3)*Eq4+Eq5+Eq6;
}

double H2(const double *xx)
{
  const double M     = xx[0];
  const double K1    = xx[1];
  const double K2    = xx[2];
  const double Phi   = xx[3];
  const double Alpha = xx[4];

  const double Eq1   = -1/(M*sin(Phi-Alpha));
  const double Eq2   = (4*M_PI*M*M-2*K1-4*K2)*cos(Phi)*sin(Phi);
  const double Eq3   = 4*K2*cos(Phi)*sin(Phi)*sin(Phi)*sin(Phi);
  const double Eq4   = M*cos(Alpha-Phi);
  const double Eq5   = (4*M_PI*M*M-2*K1-4*K2)*sin(Phi)*sin(Phi);
  const double Eq6   = 4*K2*sin(Phi)*sin(Phi)*sin(Phi)*sin(Phi);


  return Eq1*(Eq2+Eq3)*Eq4-Eq5-Eq6;
}


double F(const double *xx)
{
  const double M     = xx[0];
  return 1.76*10000000/(2*M_PI*M)*sqrt(H1(xx)*H2(xx));	
}
 
int main()
{
  double xx[5] ={2.0,2.0,2.0,2.0,1.0};
  std::cout << "H1: " << H1(xx) << std::endl;
  std::cout << "H2: " << H2(xx) << std::endl;
  std::cout << "H: "  << H(xx)  << std::endl;
  std::cout << "F: "  << F(xx)  << std::endl;
}
