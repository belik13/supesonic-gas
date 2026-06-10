// params.hpp — входные параметры задачи и численные настройки решателя.
#pragma once
namespace gd::params {

// ---------- термодинамика воздуха ----------
constexpr double X_N2 = 0.79;
constexpr double X_O2 = 0.21;
constexpr double R_AIR_CONST = 287.06;
constexpr double GAMMA_CONST = 1.4;

// ---------- базовые условия набегающих потоков ----------
constexpr double T0 = 300.0;
constexpr double P0 = 1.0e5;

// ---------- автомодельные профили трёх конфигураций ----------
constexpr double PROFILE_M    = 3.0;
constexpr double SS_ALPHA_DEG = 12.0;
constexpr double RR_ALPHA_DEG = 10.0;
constexpr double SR_PRATIO    = 2.0;

// ---------- карта областей существования решения ----------
constexpr double MAP_M_MIN     = 1.6,  MAP_M_MAX   = 8.0,  MAP_M_STEP   = 0.1;
constexpr double MAP_DTH_MIN   = -80.0, MAP_DTH_MAX = 60.0, MAP_DTH_STEP = 1.0;
constexpr double BNDS_M_MIN    = 1.8,  BNDS_M_MAX  = 7.0,  BNDS_M_STEP  = 0.2;
constexpr double BNDS_MACH_HI  = 179.0;
constexpr double BNDS_VAC_HI   = -200.0;

// ---------- зависимости углов от начальных параметров ----------
constexpr double ANG_SYM_M         = 3.0;
constexpr double ANG_SYM_ALPHA_MIN = 1.0;
constexpr double ANG_SYM_ALPHA_MAX = 32.0;
constexpr double ANG_SYM_ALPHA_STEP = 1.0;
constexpr double ANG_PR_M         = 3.0;
constexpr double ANG_PR_MIN       = 0.2;
constexpr double ANG_PR_MAX       = 5.0;
constexpr double ANG_PR_FACTOR    = 1.07;
constexpr double ANG_VSM_ALPHA_DEG = 20.0;
constexpr double ANG_VSM_M_MIN    = 3.0;
constexpr double ANG_VSM_M_MAX    = 10.0;
constexpr double ANG_VSM_M_STEP   = 0.25;

// ---------- сравнение профилей при высокой T (cp(T) vs gamma=const) ----------
constexpr double CMP_M         = 4.0;
constexpr double CMP_ALPHA_DEG = 18.0;

// ---------- произвольные асимметричные конфигурации ----------
// Левый и правый потоки каждой строки задаются полностью независимо.
// Добавляй/меняй строки в таблице CUSTOM_CASES — пересобирай — запускай:
// для каждой строки получишь профиль results/<name> (модель cp(T)).
constexpr bool CUSTOM_ENABLED = true;
struct CustomCase {
  const char* name;                 // имя выходного CSV
  double pL, TL, ML, thL_deg;       // левый поток:  давление[Па], T[К], M, угол[град]
  double pR, TR, MR, thR_deg;       // правый поток: давление[Па], T[К], M, угол[град]
};
inline constexpr CustomCase CUSTOM_CASES[] = {
  // name                  pL     TL   ML   thL    pR     TR   MR   thR
  {"prof_custom1.csv",   1.0e6, 300, 5.0,  +20,  2.0e5, 300, 3.0,  -10}, // разрежение–разрежение (расхождение ±18°)
  {"prof_custom2.csv",   3.0e5, 400, 3.0,   +5,  1.0e5, 300, 3.5,   -2}, // смешанная: скачок + разрежение, наклонный контакт
  {"prof_custom3.csv",   1.0e5, 300, 3.0,  +40,  1.0e5, 300, 3.0,  -40}, // режим Маха (регулярного решения нет)
};

// ---------- численные параметры решателя ----------
constexpr int    SCAN_NODES      = 160;
constexpr int    BISECT_ITER     = 80;
constexpr double PSTAR_REL_TOL   = 1.0e-10;
constexpr int    NEWTON_ITER     = 100;
constexpr double NEWTON_REL_TOL  = 1.0e-10;
constexpr int    SIMPSON_N       = 30;
constexpr int    BNDS_BISECT_ITER = 24;

// ---------- технические ----------
constexpr int    PROFILE_FAN_N  = 30;
constexpr double PROFILE_PAD_DEG = 12.0;

} // namespace gd::params
