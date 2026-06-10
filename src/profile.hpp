// profile.hpp — построение автомодельного профиля решения по полярному углу.
#pragma once
#include "riemann.hpp"
#include "../input/params.hpp"
#include <algorithm>
#include <cmath>
#include <ostream>
#include <vector>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace gd {

struct ProfilePoint { double phi, p, V, T, M, theta; };

struct ProfileResult {
  bool ok;
  RiemOut riem;
  std::vector<ProfilePoint> pts;
};

inline ProfileResult build_profile(const Gas& g, const State& L, const State& R) {
  using namespace gd::params;
  const double DEG = M_PI/180.0;
  ProfileResult res; res.ok = false;
  RiemOut o = steady_riemann(g, L, R);
  if (!o.ok) return res;
  res.riem = o;

  const double thL = theta_of(L), thR = theta_of(R);
  const double VL  = Vmag(L),    VR  = Vmag(R);
  const double ML  = VL/g.a(L.T), MR = VR/g.a(R.T);
  const double muL = std::asin(1.0/ML), muR = std::asin(1.0/MR);
  const double thS = o.thetaStar, pstar = o.pstar;
  const double MLs = Vmag(o.sL)/g.a(o.sL.T), MRs = Vmag(o.sR)/g.a(o.sR.T);
  const double muLs = std::asin(1.0/MLs), muRs = std::asin(1.0/MRs);

  auto add = [&](double phi, const State& s){
    double V = Vmag(s);
    res.pts.push_back({phi, s.p, V, s.T, V/g.a(s.T), theta_of(s)});
  };
  auto addRaw = [&](double phi, double p, double V, double T, double M, double th){
    res.pts.push_back({phi, p, V, T, M, th});
  };

  // нижняя волна
  double loWaveLo, loWaveHi;
  if (o.waveL == 1) { // косой скачок
    ShockOut so = oblique_shock(g, L.T, L.p, VL, pstar/L.p);
    double phis = thL - so.sigma; loWaveLo = phis; loWaveHi = phis;
  } else {            // волна разрежения
    loWaveLo = thL - muL; loWaveHi = thS - muLs;
  }
  add(loWaveLo - PROFILE_PAD_DEG*DEG, L);
  add(loWaveLo - 1e-4, L);
  if (o.waveL == 1) {
    add(loWaveLo - 1e-5, L); add(loWaveLo + 1e-5, o.sL);
  } else {
    const int N = PROFILE_FAN_N;
    for (int i = 0; i <= N; i++) {
      double X = 1.0 + (pstar/L.p - 1.0)*i/N;
      ExpOut e = prandtl_meyer(g, L.T, L.p, VL, X);
      double th = thL + e.dnu, mu = std::asin(1.0/e.M2);
      addRaw(th - mu, e.p2, e.V2, e.T2, e.M2, th);
    }
  }
  // зона *L
  add(loWaveHi + 1e-4, o.sL); add(thS - 1e-4, o.sL);
  // контактный разрыв
  add(thS - 1e-5, o.sL); add(thS + 1e-5, o.sR);
  // зона *R
  double upWaveLo, upWaveHi;
  if (o.waveR == 1) {
    ShockOut so = oblique_shock(g, R.T, R.p, VR, pstar/R.p);
    double phis = thR + so.sigma; upWaveLo = phis; upWaveHi = phis;
  } else {
    upWaveLo = thS + muRs; upWaveHi = thR + muR;
  }
  add(thS + 1e-4, o.sR); add(upWaveLo - 1e-4, o.sR);
  // верхняя волна
  if (o.waveR == 1) {
    add(upWaveLo - 1e-5, o.sR); add(upWaveLo + 1e-5, R);
  } else {
    const int N = PROFILE_FAN_N;
    for (int i = 0; i <= N; i++) {
      double X = pstar/R.p + (1.0 - pstar/R.p)*i/N;
      ExpOut e = prandtl_meyer(g, R.T, R.p, VR, X);
      double th = thR - e.dnu, mu = std::asin(1.0/e.M2);
      addRaw(th + mu, e.p2, e.V2, e.T2, e.M2, th);
    }
  }
  add(upWaveHi + 1e-4, R);
  add(upWaveHi + PROFILE_PAD_DEG*DEG, R);

  std::sort(res.pts.begin(), res.pts.end(),
            [](const ProfilePoint& a, const ProfilePoint& b){ return a.phi < b.phi; });
  res.ok = true;
  return res;
}

inline void write_profile_csv(std::ostream& f, const ProfileResult& r) {
  const double DEG = M_PI/180.0;
  f << "phi,p,V,T,M,theta\n";
  for (const auto& q : r.pts)
    f << q.phi/DEG << ',' << q.p << ',' << q.V << ',' << q.T << ','
      << q.M << ',' << q.theta/DEG << '\n';
}

} // namespace gd
