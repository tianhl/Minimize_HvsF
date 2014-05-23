#include <iostream>
#include <stdlib.h>   
#include <map>
#include <math.h>
#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "TRandom2.h"
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
  const double Eq6   = 4*K2*(3*sin(Phi)*sin(Phi)*cos(Phi)*cos(Phi)-sin(Phi)*sin(Phi)*sin(Phi)*sin(Phi));

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
  return 1.76*1e7/(2*M_PI*M)*sqrt(H1(xx)*H2(xx));	
}


void HvsF(const double *xx, double *h, double *f){
  *h = H(xx);
  *f = F(xx);
}

double AbsDiffToH(const double *pars)
{
  double xx[5] = {pars[0], pars[1], pars[2], pars[3], pars[4]};
  double hpoint = pars[5];

  return abs(H(xx)-hpoint);
}
 
int main()
{

  const double dataPointH[1] = {3041.38127};
  const double dataPointF[1] = {7.09871e9};

  double m     = 1257.;
  double alpha = 0.1651;
  double k1    = 1e7;
  double k2    = 5e5;
  double phi   = alpha;

  double xx[5] = {m, k1, k2, phi, alpha};

  // find H point

  const double hpoint  = 3041.38;
  const string minName = "Minuit2";
  const string algName = "";

  ROOT::Math::Functor funH(&AbsDiffToH,6); 

  // initialize minimizer
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer(minName.c_str(), algName.c_str());
  min->SetMaxFunctionCalls(1000000); 
  min->SetMaxIterations(10000);  
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  double minstep[6]    = {0.0,  0.0,  0.0,  0.1,  0.0,    0.0};
  double startpoint[6] = {m,    k1,   k2,   phi,  alpha,  hpoint};
  int    randomSeed = 100;

  TRandom2 r(randomSeed);
  startpoint[3] = r.Uniform(alpha,M_PI/2);
  
  min->SetFunction(funH);

  // Set the free variables to be minimized!
  min->SetVariable(0,"m",     startpoint[0], minstep[0]);
  min->SetVariable(1,"k1",    startpoint[1], minstep[1]);
  min->SetVariable(2,"k2",    startpoint[2], minstep[2]);
  min->SetVariable(3,"phi",   startpoint[3], minstep[3]);
  min->SetVariable(4,"alpha", startpoint[4], minstep[4]);
  min->SetVariable(5,"hpoint",startpoint[5], minstep[5]);

  // do the minimization
  min->Minimize();

  const double *xs = min->X();
  std::cout << "m:      " << xs[0] << std::endl; 
  std::cout << "k1:     " << xs[1] << std::endl; 
  std::cout << "k2:     " << xs[2] << std::endl; 
  std::cout << "phi:    " << xs[3] << std::endl; 
  std::cout << "alpha:  " << xs[4] << std::endl; 
  std::cout << "hpoint: " << xs[5] << std::endl; 
  std::cout << "Diff2point: " << min->MinValue()  << std::endl;

  xx[0]=xs[0]; 
  xx[1]=xs[1]; 
  xx[2]=xs[2]; 
  xx[3]=xs[3]; 
  xx[4]=xs[4]; 

  std::cout << "H point: " << H(xx) << std::endl;
  std::cout << "F point: " << F(xx) << std::endl;

  // Draw Plot 

  typedef std::map<double, double> HFMap;
  typename HFMap::iterator hfmap_it;
  HFMap hfmap;

  const int Nstep    = 100;
  double phi_step = (M_PI/2 - alpha)/(Nstep+2);

  double harray[Nstep];
  double farray[Nstep];

  int i = 0;
  for(i = 0; i < Nstep; i++){
    xx[3]    = xx[3] + phi_step;
    double h, f;    
    HvsF(xx, &h, &f);
    hfmap.insert(std::pair<double,double>(h,f));
    std::cout << "H: " << h << " F: " << f << std::endl;
  }

  i = 0;
  for(hfmap_it=hfmap.begin(); hfmap_it != hfmap.end(); hfmap_it++){
    harray[i] = hfmap_it->first;
    farray[i] = hfmap_it->second;
    i++;
  }


  TCanvas *cav  = new TCanvas("c1");
  TGraph  *hvsf = new TGraph(Nstep, harray, farray);
  TMarker *data = new TMarker(dataPointH[0], dataPointF[0], 20);
  data->SetMarkerColor(2);
  hvsf->SetMarkerStyle(21);
  hvsf->Draw("APL");
  data->Draw();
  cav->Print("FunGraph.eps");


}






