// gas.hpp — модель газа: совершенный газ с переменной cp(T) по полиномам NASA.
#pragma once
#include <cmath>
namespace gd {
constexpr double Ru    = 8.314462618;
constexpr double M_AIR = 0.0289647;
struct Nasa { double lo[7], hi[7]; double Tmid; };
inline const Nasa& N2(){ static const Nasa s = {
  {3.298677e0,1.4082404e-3,-3.963222e-6,5.641515e-9,-2.444854e-12,-1.0208999e3,3.950372e0},
  {2.92664e0,1.4879768e-3,-5.68476e-7,1.0097038e-10,-6.753351e-15,-9.227977e2,5.980528e0},1000.0}; return s; }
inline const Nasa& O2(){ static const Nasa s = {
  {3.212936e0,1.127486e-3,-5.75615e-7,1.3138773e-9,-8.768554e-13,-1.005249e3,6.034738e0},
  {3.697578e0,6.135197e-4,-1.258842e-7,1.775281e-11,-1.136435e-15,-1.23393e3,3.189166e0},1000.0}; return s; }
inline const double* coef(const Nasa& s,double T){ return T<s.Tmid? s.lo:s.hi; }
inline double cp_mol(const Nasa& s,double T){const double*a=coef(s,T);
  return Ru*(a[0]+a[1]*T+a[2]*T*T+a[3]*T*T*T+a[4]*T*T*T*T);}
inline double h_mol(const Nasa& s,double T){const double*a=coef(s,T);
  return Ru*(a[0]*T+a[1]/2*T*T+a[2]/3*T*T*T+a[3]/4*T*T*T*T+a[4]/5*T*T*T*T*T+a[5]);}
inline double s0_mol(const Nasa& s,double T){const double*a=coef(s,T);
  return Ru*(a[0]*log(T)+a[1]*T+a[2]/2*T*T+a[3]/3*T*T*T+a[4]/4*T*T*T*T+a[6]);}
class Gas {
public:
  bool variable; double R; double cp_c{},gamma_c{}; double xN2{0.79},xO2{0.21};
  static Gas Air(){ Gas g; g.variable=true; g.R=Ru/M_AIR; return g; }
  static Gas Perfect(double gamma=1.4,double R=287.06){
    Gas g; g.variable=false; g.R=R; g.gamma_c=gamma; g.cp_c=gamma*R/(gamma-1.0); return g; }
  double cp(double T) const { if(!variable) return cp_c;
    return (xN2*cp_mol(N2(),T)+xO2*cp_mol(O2(),T))/M_AIR; }
  double h(double T) const { if(!variable) return cp_c*T;
    return (xN2*h_mol(N2(),T)+xO2*h_mol(O2(),T))/M_AIR; }
  double phiR(double T) const { if(!variable) return (cp_c/R)*log(T);
    return (xN2*s0_mol(N2(),T)+xO2*s0_mol(O2(),T))/Ru; }
  double gamma(double T) const { if(!variable) return gamma_c; double c=cp(T); return c/(c-R); }
  double a(double T) const { return sqrt(gamma(T)*R*T); }
  double T_from_h(double htarget,double Tguess=300.0) const { if(!variable) return htarget/cp_c;
    double T=(Tguess>50?Tguess:300.0);
    for(int i=0;i<100;i++){ double f=h(T)-htarget,fp=cp(T),dT=-f/fp;
      if(T+dT<10) dT=10-T; T+=dT; if(std::fabs(dT)<1e-9*T) break; } return T; }
};
} // namespace gd
