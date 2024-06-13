#pragma once

#include <random>
#include <thread>

static thread_local std::mt19937 engine(
    *INI.Seed() + std::hash<std::thread::id>{}(std::this_thread::get_id()));
static thread_local uint32_t rand_r_state = engine();

inline double UegoRand() {
  static thread_local std::uniform_real_distribution<double> dist(0.0, 1.0);
  return dist(engine);
}

inline long UegoLongRand(long low, long high) {
  static thread_local std::uniform_int_distribution<int> dist(low, high);
  return dist(engine);
}

inline double UegoDoubleRand(double low, double high) {
  static thread_local std::uniform_real_distribution<double> dist(low, high);
  return dist(engine);
}

inline double UegoDoubleRandNoNormalized(double low, double high) {
  return (low + ((double)rand_r(&rand_r_state) / RAND_MAX) * (high - low));
}
