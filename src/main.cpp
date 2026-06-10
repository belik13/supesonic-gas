// main.cpp — драйвер: создаёт модели газа и запускает все эксперименты.
#include "experiments.hpp"
#include "../input/params.hpp"
#include <cstdio>

int main() {
  using namespace gd;
  Gas pg  = Gas::Perfect(params::GAMMA_CONST, params::R_AIR_CONST);
  Gas air = Gas::Air();

  exper::run_profiles_basic(pg);
  exper::run_existence_map(pg);
  exper::run_existence_bounds(pg, air);
  exper::run_angles_symmetric(pg, air);
  exper::run_angles_pratio(pg);
  exper::run_angles_vsM(pg, air);
  exper::run_comparison_profiles(pg, air);
  exper::run_custom_case(air);

  std::printf("\nГотово. Результаты в %s\n", io::RESULTS_DIR);
  return 0;
}
