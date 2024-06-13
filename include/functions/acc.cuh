#pragma once

#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int32_t

#include "uego/uegoini.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#define CUDA_CHECK(call)                                                       \
  if ((call) != cudaSuccess) {                                                 \
    cudaError_t err = cudaGetLastError();                                      \
    cerr << "CUDA error calling \"" #call "\", code is "                       \
         << cudaGetErrorName(err) << endl;                                     \
    exit(err);                                                                 \
  }

typedef Eigen::Vector3d vec3d;
typedef Eigen::Quaterniond quatd;

class Acc {
private:
  Acc();

  // Props. for enabling or not CUDA
  bool use_cuda;
  uint32_t max_CUDA_devices;
  uint32_t max_threads;

  // CUDA query molecule
  vec3d *query_atoms;
  double *query_weight;
  double *query_radius;
  // CUDA target/variable molecule
  vec3d *target_atoms;
  double *target_weight;
  double *target_radius;

public:
  Acc(Acc const &) = delete;
  void operator=(Acc const &) = delete;

  static Acc &get_self() {
    static Acc self{};
    return self;
  }

  void initialize(Ini *ini);
  void upload_molecules(Ini *ini);
  bool is_using_CUDA() { return this->use_cuda; }

  double ndim_value(const double *x);
};