#include <iostream>
#include <map>
#include <math.h>
#include "TGraph.h"
using namespace std;


double H(const double *xx)
{
  const double M     = xx[0];
  const double K1    = xx[1];
  const double K2    = xx[2];
  const double Phi   = xx[3];
  const double Alpha = xx[4];

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

  double m     = 1257.;
  double alpha = 0.1651;
  double k1    = 5000000.;
  double k2    = 5000.;
  double phi   = alpha;

  double xx[5] = {m, k1, k2, phi, alpha};
  typedef std::map<double, double> HFMap;
  typename HFMap::iterator hfmap_it;
  HFMap hfmap;
  
  double nstep    = 100;
  double phi_step = (M_PI/2 - alpha)/nstep;
  for(int i = 0; i < nstep; i++){
    double h = H(xx);    
    double f = F(xx);
    xx[3]    = xx[3] + phi_step;
    hfmap.insert(std::pair<double,double>(h,f));
  }

  for(hfmap_it=hfmap.begin(); hfmap_it != hfmap.end(); hfmap_it++){
    std::cout << "H vs F: "  << hfmap_it->first  << " " << hfmap_it->second << std::endl;
  }

}






