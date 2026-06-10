// riemann.hpp — точное решение задачи о распаде разрыва для двух сверхзвуковых потоков.
#pragma once
#include "gas.hpp"
#include <algorithm>
namespace gd {
struct State { double rho,u,v,p,T; };
inline State make_state(const Gas& g,double p,double T,double V,double theta){
  State s; s.p=p; s.T=T; s.rho=p/(g.R*T); s.u=V*cos(theta); s.v=V*sin(theta); return s; }
inline double Vmag(const State& s){ return std::sqrt(s.u*s.u+s.v*s.v); }
inline double theta_of(const State& s){ return std::atan2(s.v,s.u); }

struct ShockOut{ double T2,p2,V2,M2,sigma,turn; bool ok; };
inline ShockOut oblique_shock(const Gas& g,double T1,double p1,double V1,double X){
  ShockOut o; o.ok=false; const double R=g.R; const double p2=X*p1; const double h1=g.h(T1);
  double T2=T1*X;
  for(int i=0;i<100;i++){ double f=g.h(T2)-h1-0.5*R*(X-1.0)*(T1+T2/X);
    double fp=g.cp(T2)-0.5*R*(X-1.0)/X; double dT=-f/fp; T2+=dT; if(std::fabs(dT)<1e-10*T2) break; }
  double rho1=p1/(R*T1),rho2=p2/(R*T2),r=rho1/rho2,denom=1.0-r; if(denom<=0) return o;
  double un1sq=(p2-p1)/(rho1*denom); if(un1sq<=0||un1sq>V1*V1) return o;
  double un1=std::sqrt(un1sq),sigma=std::asin(un1/V1),ut=V1*std::cos(sigma),un2=rho1*un1/rho2;
  o.T2=T2;o.p2=p2;o.V2=std::sqrt(un2*un2+ut*ut);o.M2=o.V2/g.a(T2);
  o.sigma=sigma;o.turn=sigma-std::atan2(un2,ut);o.ok=true; return o; }

struct ExpOut{ double T2,p2,V2,M2,dnu; bool ok; };
inline ExpOut prandtl_meyer(const Gas& g,double T1,double p1,double V1,double X){
  ExpOut o;o.ok=false; const double R=g.R; const double lnX=std::log(X);
  double T2=T1*std::pow(X,0.2857);
  for(int i=0;i<100;i++){ double f=g.phiR(T2)-g.phiR(T1)-lnX,fp=g.cp(T2)/(R*T2),dT=-f/fp;
    T2+=dT; if(std::fabs(dT)<1e-10*T2) break; }
  const double H=g.h(T1)+0.5*V1*V1; double V2sq=2*(H-g.h(T2)); if(V2sq<=0) return o;
  double V2=std::sqrt(V2sq); const int n=30; double dV=(V2-V1)/n,s=0; double Tcur=T1;
  auto integ=[&](double V){ double T=g.T_from_h(H-0.5*V*V,Tcur); Tcur=T;
    double M=V/g.a(T),m2=M*M-1; return (m2>0?std::sqrt(m2):0.0)/V; };
  for(int i=0;i<=n;i++){ double V=V1+i*dV,w=(i==0||i==n)?1:(i%2?4:2); s+=w*integ(V); }
  o.T2=T2;o.p2=X*p1;o.V2=V2;o.M2=V2/g.a(T2);o.dnu=s*dV/3.0;o.ok=true; return o; }

inline bool wave_turn(const Gas& g,const State& st,double pstar,double& G,double& Tpost,double& Vpost){
  double V1=Vmag(st);
  if(pstar>=st.p){ ShockOut o=oblique_shock(g,st.T,st.p,V1,pstar/st.p);
    if(!o.ok) return false; G=o.turn;Tpost=o.T2;Vpost=o.V2; return true; }
  ExpOut o=prandtl_meyer(g,st.T,st.p,V1,pstar/st.p);
  if(!o.ok) return false; G=-o.dnu;Tpost=o.T2;Vpost=o.V2; return true; }

struct RiemOut{ double pstar,thetaStar; State sL,sR; bool ok; int waveL,waveR; };
inline RiemOut steady_riemann(const Gas& g,const State& L,const State& R){
  RiemOut o;o.ok=false; const double thL=theta_of(L),thR=theta_of(R);
  auto Phi=[&](double pstar,bool& good){ double GL,GR,t,v;
    bool a=wave_turn(g,L,pstar,GL,t,v),b=wave_turn(g,R,pstar,GR,t,v);
    good=a&&b; return (thL-GL)-(thR+GR); };
  // поиск слабого корня: первая смена знака Phi на лог-сетке давлений
  double pmin=1e-4*std::min(L.p,R.p), pmax=1e3*std::max(L.p,R.p);
  const int NS=160; double pa=-1,pb=-1,fa=0;
  bool pok; double pp=pmin, pf=Phi(pp,pok);
  for(int i=1;i<=NS;i++){
    double p=pmin*std::pow(pmax/pmin,(double)i/NS); bool ok; double f=Phi(p,ok);
    if(pok&&ok&&pf*f<0){ pa=pp; pb=p; fa=pf; break; }
    pp=p; pf=f; pok=ok; }
  if(pa<0) return o;                       // решение отсутствует (вакуум или режим Маха)
  double pstar=0.5*(pa+pb);
  for(int it=0;it<80;it++){ double pm=0.5*(pa+pb); bool gm; double fm=Phi(pm,gm);
    if(!gm){ pb=pm; continue; } if(fa*fm<=0) pb=pm; else { pa=pm; fa=fm; }
    pstar=0.5*(pa+pb); if(std::fabs(pb-pa)<1e-10*pstar) break; }
  double GL,GR,TLs,VLs,TRs,VRs; bool gL,gR;
  gL=wave_turn(g,L,pstar,GL,TLs,VLs); gR=wave_turn(g,R,pstar,GR,TRs,VRs);
  if(!gL||!gR) return o; double thStar=thL-GL;
  o.pstar=pstar;o.thetaStar=thStar; o.sL=make_state(g,pstar,TLs,VLs,thStar);
  o.sR=make_state(g,pstar,TRs,VRs,thStar);
  o.waveL=(pstar>=L.p)?1:0; o.waveR=(pstar>=R.p)?1:0; o.ok=true; return o; }
} // namespace gd
