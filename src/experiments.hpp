// experiments.hpp — численные эксперименты: профили, карты, углы, сравнение моделей газа.
#pragma once
#include "riemann.hpp"
#include "profile.hpp"
#include "io.hpp"
#include "../input/params.hpp"
#include <cstdio>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
namespace gd::exper {

inline constexpr double DEG = M_PI/180.0;

inline void run_profiles_basic(const Gas& g) {
  using namespace gd::params;
  const double V = PROFILE_M * g.a(T0);
  auto write = [&](const State& L, const State& R, const char* name) {
    auto r = build_profile(g, L, R);
    if (!r.ok) { std::printf("  %s: решение отсутствует\n", name); return; }
    auto f = io::open_csv(name); write_profile_csv(f, r);
    std::printf("  %s: p*/p_L=%.4f theta*=%+0.3f deg (waves L%d R%d)\n",
                name, r.riem.pstar/L.p, r.riem.thetaStar/DEG, r.riem.waveL, r.riem.waveR);
  };
  std::printf("=== Автомодельные профили (gamma=1.4, M=%g) ===\n", PROFILE_M);
  write(make_state(g, P0, T0, V,  SS_ALPHA_DEG*DEG),
        make_state(g, P0, T0, V, -SS_ALPHA_DEG*DEG), "prof_SS.csv");
  write(make_state(g, P0, T0, V, -RR_ALPHA_DEG*DEG),
        make_state(g, P0, T0, V,  RR_ALPHA_DEG*DEG), "prof_RR.csv");
  write(make_state(g, SR_PRATIO*P0, T0, V, 0.0),
        make_state(g, P0,            T0, V, 0.0), "prof_SR.csv");
}

inline void run_existence_map(const Gas& g) {
  using namespace gd::params;
  auto f = io::open_csv("exist_map.csv");
  f << "M,dtheta,code\n";
  for (double M = MAP_M_MIN; M <= MAP_M_MAX + 1e-9; M += MAP_M_STEP) {
    double V = M * g.a(T0);
    for (double dth = MAP_DTH_MIN; dth <= MAP_DTH_MAX + 1e-9; dth += MAP_DTH_STEP) {
      State L = make_state(g, P0, T0, V,  0.5*dth*DEG);
      State R = make_state(g, P0, T0, V, -0.5*dth*DEG);
      RiemOut o = steady_riemann(g, L, R); int code;
      if (o.ok) code = (o.waveL==1 && o.waveR==1) ? 3
                     : (o.waveL==0 && o.waveR==0) ? 1 : 2;
      else      code = (dth >= 0) ? 4 : 0;
      f << M << ',' << dth << ',' << code << '\n';
    }
  }
  std::printf("=== Карта областей существования записана ===\n");
}

inline void run_existence_bounds(const Gas& pg, const Gas& air) {
  using namespace gd::params;
  auto ok = [&](const Gas& g, double M, double dt) {
    double V = M * g.a(T0);
    State L = make_state(g, P0, T0, V, 0.5*dt*DEG);
    State R = make_state(g, P0, T0, V,-0.5*dt*DEG);
    return steady_riemann(g, L, R).ok;
  };
  auto bisect = [&](const Gas& g, double M, double hi) {
    double lo = 0.0;
    if (!ok(g, M, lo)) return 0.0;
    for (int i = 0; i < BNDS_BISECT_ITER; i++) {
      double m = 0.5*(lo + hi);
      if (ok(g, M, m)) lo = m; else hi = m;
    }
    return lo;
  };
  auto f = io::open_csv("exist_bounds.csv");
  f << "M,dth_mach,dth_vac,dth_mach_air,dth_vac_air\n";
  for (double M = BNDS_M_MIN; M <= BNDS_M_MAX + 1e-9; M += BNDS_M_STEP) {
    f << M << ',' << bisect(pg, M, BNDS_MACH_HI) << ','
                  << bisect(pg, M, BNDS_VAC_HI)  << ','
                  << bisect(air,M, BNDS_MACH_HI) << ','
                  << bisect(air,M, BNDS_VAC_HI)  << '\n';
  }
  std::printf("=== Границы существования записаны ===\n");
}

inline void run_angles_symmetric(const Gas& pg, const Gas& air) {
  using namespace gd::params;
  auto f = io::open_csv("ang_sym.csv");
  f << "alpha,pstar_c,sigma_c,T2_c,pstar_v,sigma_v,T2_v\n";
  const double Vp = ANG_SYM_M * pg.a(T0);
  const double Va = ANG_SYM_M * air.a(T0);
  for (double a = ANG_SYM_ALPHA_MIN; a <= ANG_SYM_ALPHA_MAX + 1e-9; a += ANG_SYM_ALPHA_STEP) {
    State Lc = make_state(pg, P0, T0, Vp, a*DEG),  Rc = make_state(pg, P0, T0, Vp, -a*DEG);
    State Lv = make_state(air,P0, T0, Va, a*DEG),  Rv = make_state(air,P0, T0, Va, -a*DEG);
    RiemOut oc = steady_riemann(pg, Lc, Rc);
    RiemOut ov = steady_riemann(air,Lv, Rv);
    double sc=0, sv=0, Tc=0, Tv=0, pc=0, pv=0;
    if (oc.ok) { ShockOut s = oblique_shock(pg, T0, P0, Vp, oc.pstar/P0);
                 sc = s.sigma/DEG; Tc = oc.sL.T; pc = oc.pstar/P0; }
    if (ov.ok) { ShockOut s = oblique_shock(air,T0, P0, Va, ov.pstar/P0);
                 sv = s.sigma/DEG; Tv = ov.sL.T; pv = ov.pstar/P0; }
    f << a << ',' << pc << ',' << sc << ',' << Tc << ','
                  << pv << ',' << sv << ',' << Tv << '\n';
  }
  std::printf("=== Зависимости sigma, p*/p от alpha записаны ===\n");
}

inline void run_angles_pratio(const Gas& g) {
  using namespace gd::params;
  auto f = io::open_csv("ang_pratio.csv");
  f << "ratio,pstar,thetaStar,sigma_shock\n";
  const double V = ANG_PR_M * g.a(T0);
  for (double r = ANG_PR_MIN; r <= ANG_PR_MAX + 1e-9; r *= ANG_PR_FACTOR) {
    State L = make_state(g, r*P0, T0, V, 0.0);
    State R = make_state(g, P0,   T0, V, 0.0);
    RiemOut o = steady_riemann(g, L, R); if (!o.ok) continue;
    double sig = 0;
    if      (o.waveR == 1) { ShockOut s = oblique_shock(g, R.T, R.p, Vmag(R), o.pstar/R.p); sig = s.sigma/DEG; }
    else if (o.waveL == 1) { ShockOut s = oblique_shock(g, L.T, L.p, Vmag(L), o.pstar/L.p); sig = s.sigma/DEG; }
    f << r << ',' << o.pstar << ',' << o.thetaStar/DEG << ',' << sig << '\n';
  }
  std::printf("=== Зависимости theta*, sigma от p_L/p_R записаны ===\n");
}

inline void run_angles_vsM(const Gas& pg, const Gas& air) {
  using namespace gd::params;
  auto f = io::open_csv("ang_vsM.csv");
  f << "M,pstar_c,sigma_c,T2_c,pstar_v,sigma_v,T2_v\n";
  const double a = ANG_VSM_ALPHA_DEG;
  for (double M = ANG_VSM_M_MIN; M <= ANG_VSM_M_MAX + 1e-9; M += ANG_VSM_M_STEP) {
    double Vp = M * pg.a(T0), Va = M * air.a(T0);
    State Lc = make_state(pg, P0, T0, Vp,  a*DEG), Rc = make_state(pg, P0, T0, Vp, -a*DEG);
    State Lv = make_state(air,P0, T0, Va,  a*DEG), Rv = make_state(air,P0, T0, Va, -a*DEG);
    RiemOut oc = steady_riemann(pg, Lc, Rc);
    RiemOut ov = steady_riemann(air,Lv, Rv);
    double sc=0, sv=0, pc=0, pv=0, Tc=0, Tv=0;
    if (oc.ok) { ShockOut s = oblique_shock(pg, T0, P0, Vp, oc.pstar/P0);
                 sc = s.sigma/DEG; pc = oc.pstar/P0; Tc = oc.sL.T; }
    if (ov.ok) { ShockOut s = oblique_shock(air,T0, P0, Va, ov.pstar/P0);
                 sv = s.sigma/DEG; pv = ov.pstar/P0; Tv = ov.sL.T; }
    f << M << ',' << pc << ',' << sc << ',' << Tc << ','
                  << pv << ',' << sv << ',' << Tv << '\n';
  }
  std::printf("=== Зависимости sigma, p*/p, T2 от M записаны ===\n");
}

inline void run_comparison_profiles(const Gas& pg, const Gas& air) {
  using namespace gd::params;
  const double Vp = CMP_M * pg.a(T0), Va = CMP_M * air.a(T0);
  const double a  = CMP_ALPHA_DEG;
  auto write_one = [&](const Gas& g, double V, const char* name) {
    State L = make_state(g, P0, T0, V,  a*DEG);
    State R = make_state(g, P0, T0, V, -a*DEG);
    auto r = build_profile(g, L, R);
    if (!r.ok) { std::printf("  %s: решение отсутствует\n", name); return; }
    auto f = io::open_csv(name); write_profile_csv(f, r);
    std::printf("  %s: p*/p=%.4f T2=%.1f K\n", name, r.riem.pstar/P0, r.riem.sL.T);
  };
  std::printf("=== Сравнение профилей при M=%g, alpha=%g deg ===\n", CMP_M, CMP_ALPHA_DEG);
  write_one(pg,  Vp, "prof_cmp_const.csv");
  write_one(air, Va, "prof_cmp_air.csv");
}

inline void run_custom_case(const Gas& g) {
  using namespace gd::params;
  if (!CUSTOM_ENABLED) return;
  std::printf("=== Произвольные асимметричные конфигурации (модель cp(T)) ===\n");
  for (const auto& c : CUSTOM_CASES) {
    const double VL = c.ML * g.a(c.TL);
    const double VR = c.MR * g.a(c.TR);
    State L = make_state(g, c.pL, c.TL, VL, c.thL_deg*DEG);
    State R = make_state(g, c.pR, c.TR, VR, c.thR_deg*DEG);
    std::printf("  [%s]\n", c.name);
    std::printf("    L: p=%.0f Pa  T=%.1f K  M=%.2f  theta=%+0.2f deg\n",
                c.pL, c.TL, c.ML, c.thL_deg);
    std::printf("    R: p=%.0f Pa  T=%.1f K  M=%.2f  theta=%+0.2f deg\n",
                c.pR, c.TR, c.MR, c.thR_deg);
    auto r = build_profile(g, L, R);
    if (!r.ok) {
      std::printf("    -> решение отсутствует (вакуум или режим Маха)\n");
      continue;
    }
    auto f = io::open_csv(c.name);
    write_profile_csv(f, r);
    std::printf("    -> p*/p_L=%.4f  theta*=%+0.3f deg  waves L%d R%d  T*_L=%.1f K  T*_R=%.1f K\n",
                r.riem.pstar/c.pL, r.riem.thetaStar/DEG,
                r.riem.waveL, r.riem.waveR, r.riem.sL.T, r.riem.sR.T);
  }
}

} // namespace gd::exper
