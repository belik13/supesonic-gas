// io.hpp — открытие CSV-файлов в директории результатов.
#pragma once
#include <fstream>
#include <string>
namespace gd::io {
constexpr const char* RESULTS_DIR = "results/";
inline std::ofstream open_csv(const std::string& name) {
  return std::ofstream(std::string(RESULTS_DIR) + name);
}
} // namespace gd::io
