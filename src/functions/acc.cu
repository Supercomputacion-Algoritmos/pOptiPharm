#include "functions/acc.cuh"

#include "functions/Tanimoto.h"
#include "uego/uego.h"

Acc::Acc() {
  this->use_cuda = false;

  this->query_atoms = nullptr;
  this->query_weight = nullptr;
  this->query_radius = nullptr;

  this->target_atoms = nullptr;
  this->target_weight = nullptr;
  this->target_radius = nullptr;
}

void Acc::initialize(Ini *ini) {
  // Query CUDA props.
  int max_devs;
  CUDA_CHECK(cudaGetDeviceCount(&max_devs));

  this->max_CUDA_devices = (uint32_t)max_devs;
  this->max_threads = ini->getThreadNumber();

  uint32_t query_size = ini->getMolQuery()->getNumAtoms();
  uint32_t target_size = ini->getMolVariable()->getNumAtoms();
  uint32_t combined_size = query_size + target_size;

  if (ini->isForcingCUDA() || combined_size >= ini->getCUDAMinCombinedSize()) {
    this->use_cuda = true;

    printf("cuda     : -- GPU support enabled. Molecules combined size (%d, "
           "Query: %d + Target: %d) greater or equal to limit (%d).CPU\n",
           combined_size, query_size, target_size,
           ini->getCUDAMinCombinedSize());

    if (this->max_threads > this->max_CUDA_devices) {
      printf("cuda     : -- There are more CPU threads (%d) than CUDA devices "
             "(%d). "
             "Some CPUthreads will use the same CUDA device. \n"
             "              It MAY incur in a performance penalty. Keep "
             "watching!\n",
             this->max_threads, this->max_CUDA_devices);
    }

    for (uint32_t i = 0; i < this->max_CUDA_devices; i++) {
      cudaDeviceProp props;
      CUDA_CHECK(cudaGetDeviceProperties(&props, i));

      printf("cuda     : -- Using CUDA device %d: %s - %.2lf GiB\n", i,
             props.name,
             ceil(props.totalGlobalMem / 1024.0f / 1024.0f / 1024.0f));
    }
  } else {
    printf("cuda     : -- GPU support disabled. Molecules combined size (%d, "
           "Query: %d + Target: %d) is smaller than the limit (%d).\n",
           combined_size, query_size, target_size,
           ini->getCUDAMinCombinedSize());
    printf(
        "cuda     : -- GPU support can be forced using -fgpu 1 or by changing "
        "the min. combined size via -mc <limit>.\n");
  }
}

void Acc::upload_molecules(Ini *ini) {
  if (!this->is_using_CUDA()) {
    return;
  }

  // Query sizes
  double *h_query_atoms = ini->getMolQuery()->getAtomsXYZ();
  double *h_query_weight = ini->getMolQuery()->getWeightAtoms();
  double *h_query_radius = ini->getMolQuery()->getRadiusAtoms();
  uint32_t query_bytes = ini->getMolQuery()->getNumAtoms() * sizeof(vec3d);
  uint32_t query_bytes_rem = ini->getMolQuery()->getNumAtoms() * sizeof(double);

  // Async stream
  cudaStream_t stream;
  CUDA_CHECK(cudaStreamCreate(&stream));

  // Query allocations
  CUDA_CHECK(cudaMallocAsync(&this->query_atoms, query_bytes, stream));
  CUDA_CHECK(cudaMallocAsync(&this->query_weight, query_bytes_rem, stream));
  CUDA_CHECK(cudaMallocAsync(&this->query_radius, query_bytes_rem, stream));
  // Query copies
  CUDA_CHECK(cudaMemcpyAsync(this->query_atoms, h_query_atoms, query_bytes,
                             cudaMemcpyHostToDevice, stream));
  CUDA_CHECK(cudaMemcpyAsync(this->query_weight, h_query_weight,
                             query_bytes_rem, cudaMemcpyHostToDevice, stream));
  CUDA_CHECK(cudaMemcpyAsync(this->query_radius, query_radius, query_bytes_rem,
                             cudaMemcpyHostToDevice, stream));

  // Target sizes
  double *h_target_atoms = ini->getMolVariable()->getAtomsXYZ();
  double *h_target_weight = ini->getMolVariable()->getWeightAtoms();
  double *h_target_radius = ini->getMolVariable()->getRadiusAtoms();
  uint32_t target_bytes = ini->getMolVariable()->getNumAtoms() * sizeof(vec3d);
  uint32_t target_bytes_rem =
      ini->getMolVariable()->getNumAtoms() * sizeof(double);

  // Target allocations
  CUDA_CHECK(cudaMallocAsync(&this->target_atoms, target_bytes, stream));
  CUDA_CHECK(cudaMallocAsync(&this->target_weight, target_bytes_rem, stream));
  CUDA_CHECK(cudaMallocAsync(&this->target_radius, target_bytes_rem, stream));

  // Target upload
  CUDA_CHECK(cudaMemcpyAsync(this->target_atoms, h_target_atoms, target_bytes,
                             cudaMemcpyHostToDevice, stream));
  CUDA_CHECK(cudaMemcpyAsync(this->target_weight, h_target_weight,
                             target_bytes_rem, cudaMemcpyHostToDevice, stream));
  CUDA_CHECK(cudaMemcpyAsync(this->target_radius, h_target_radius,
                             target_bytes_rem, cudaMemcpyHostToDevice, stream));

  CUDA_CHECK(cudaStreamSynchronize(stream));
  CUDA_CHECK(cudaStreamDestroy(stream));
}

// TODO: Keep reducing register usage (currently is 40 for sm_70 - CUDA 11.7)
// 32 or less would be ideal
// use set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} --ptxas-options=-v) with cmake
__global__ void kernelRotateMol(const vec3d *atoms, uint32_t atoms_len,
                                const vec3d *x, vec3d *new_atoms) {
  uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= atoms_len) {
    return;
  }

  vec3d vector = (x[2] - x[1]).normalized() * sin(x[0].x() * 0.5);
  vec3d temp = atoms[idx] - x[1];

  quatd q1(cos(x[0].x() * 0.5), vector.x(), vector.y(), vector.z());
  quatd atom_position(0.0, temp.x(), temp.y(), temp.z());

  quatd part2 = q1 * atom_position;
  quatd part3 = part2 * q1.conjugate();

  new_atoms[idx] = x[1] + part3.vec();
}

__global__ void kernelMolToNewPosition(vec3d *new_atoms, uint32_t atoms_len,
                                       const vec3d *x) {
  uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= atoms_len) {
    return;
  }

  new_atoms[idx] += x[3];
}

__global__ void kernelPreciseOverlapSameVDW(
    const vec3d *atoms_query, const double *weight_query,
    const double *radius_query, uint32_t query_size, const vec3d *new_atoms,
    const double *weight_target, const double *radius_target,
    uint32_t target_size, double *results) {
  uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
  uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;

  if (idx >= query_size || idy >= target_size) {
    return;
  }

  vec3d temp = atoms_query[idx] - new_atoms[idy];
  double Rij2 = temp.squaredNorm();

  double Kij = exp(Rij2 * -0.3731438999881213);
  double Vij = 24.428790199 * Kij;

  results[idy * query_size + idx] =
      weight_query[idx] * weight_target[idy] * Vij;
}

__global__ void kernelPreciseOverlapNotSameVDW(
    const vec3d *atoms_query, const double *weight_query,
    const double *radius_query, uint32_t query_size, const vec3d *new_atoms,
    const double *weight_target, const double *radius_target,
    uint32_t target_size, double *results) {
  uint32_t idx = blockDim.x * blockIdx.x + threadIdx.x;
  uint32_t idy = blockDim.y * blockIdx.y + threadIdx.y;

  if (idx >= query_size || idy >= target_size) {
    return;
  }

  vec3d temp = atoms_query[idx] - new_atoms[idy];
  double Rij2 = temp.squaredNorm();

  double alphai = 2.417972471923026 / (radius_query[idx] * radius_target[idy]);
  double alphaj = 2.417972471923026 / (pow(radius_target[idy], 2.0));

  double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
  double Vij = 7.99984656 * Kij * pow((M_PI / (alphai + alphaj)), 1.5);

  results[idy * query_size + idx] =
      weight_query[idx] * weight_target[idy] * Vij;
}

double Acc::ndim_value(const double *x) {
  // Set device
  static std::atomic_uint32_t current_dev = 0;
  CUDA_CHECK(cudaSetDevice(current_dev));
  current_dev = (current_dev + 1) % this->max_CUDA_devices;

  // Pre-data
  uint32_t query_atoms_len = INI.getMolQuery()->getNumAtoms();
  uint32_t target_atoms_len = INI.getMolVariable()->getNumAtoms();
  bool SVDW = INI.getSameVanDerWaalsRadius();
  double query_tanimoto = INI.getMolQuery()->tanimoto;
  double target_tanimoto = INI.getMolVariable()->tanimoto;

  uint32_t new_x_bytes = 4 * sizeof(vec3d);
  uint32_t new_atoms_bytes = target_atoms_len * sizeof(vec3d);
  uint32_t VABs_bytes = query_atoms_len * target_atoms_len * sizeof(double);

  // Host
  static thread_local cudaStream_t stream;
  static thread_local bool initialized = false;

  static thread_local double *h_new_x;
  static thread_local double *h_VABs;
  // Device
  static thread_local vec3d *d_x;
  static thread_local vec3d *d_new_atoms;
  static thread_local double *d_VABs;

  if (!initialized) {
    initialized = true;
    // Host
    CUDA_CHECK(cudaStreamCreate(&stream));
    h_new_x = new double[12];
    h_VABs = new double[query_atoms_len * target_atoms_len];

    // Device
    CUDA_CHECK(cudaMallocAsync(&d_x, new_x_bytes, stream));
    CUDA_CHECK(cudaMallocAsync(&d_new_atoms, new_atoms_bytes, stream));
    CUDA_CHECK(cudaMallocAsync(&d_VABs, VABs_bytes, stream));
  }

  // Data copies host
  h_new_x[0] = x[0];
  h_new_x[1] = x[0];
  memcpy(h_new_x + 2, x, 10 * sizeof(double));

  // Data copies device
  CUDA_CHECK(cudaMemcpyAsync(d_x, h_new_x, new_x_bytes, cudaMemcpyHostToDevice,
                             stream));

  // Kernel launches
  // Rotate mol.
  uint32_t rotate_block_size = 512;
  uint32_t rotate_actual_grid_size =
      (target_atoms_len + rotate_block_size - 1) / rotate_block_size;

  kernelRotateMol<<<rotate_actual_grid_size, rotate_block_size, 0, stream>>>(
      this->target_atoms, target_atoms_len, d_x, d_new_atoms);

  // MolToNewPosition
  uint32_t position_block_size = 1024;
  uint32_t position_actual_grid_size =
      (target_atoms_len + rotate_block_size - 1) / rotate_block_size;

  kernelMolToNewPosition<<<position_actual_grid_size, position_block_size, 0,
                           stream>>>(d_new_atoms, target_atoms_len, d_x);

  // Overlap
  // ! Big molecule must be always be the query one
  // Kernel are splitted because combined they were using 32 registers.
  // A bigger register usage would have had tanked some performance.
  // Splitting kernels, we assure when same Van der Waals radius are used
  // with other CUDA versions / other devices / whatever, that kernel remains
  // well under 32 registers.

  dim3 precise_block_size(32, 16);
  dim3 precise_actual_grid_size(
      (query_atoms_len + precise_block_size.x - 1) / precise_block_size.x,
      (target_atoms_len + precise_block_size.y - 1) / precise_block_size.y);
  if (SVDW) {
    kernelPreciseOverlapSameVDW<<<precise_actual_grid_size, precise_block_size,
                                  0, stream>>>(
        this->query_atoms, this->query_weight, this->query_radius,
        query_atoms_len, d_new_atoms, this->target_weight, this->target_radius,
        target_atoms_len, d_VABs);
  } else {

    kernelPreciseOverlapNotSameVDW<<<precise_actual_grid_size,
                                     precise_block_size, 0, stream>>>(
        this->query_atoms, this->query_weight, this->query_radius,
        query_atoms_len, d_new_atoms, this->target_weight, this->target_radius,
        target_atoms_len, d_VABs);
  }

  CUDA_CHECK(cudaMemcpyAsync(h_VABs, d_VABs, VABs_bytes, cudaMemcpyDeviceToHost,
                             stream));
  CUDA_CHECK(cudaStreamSynchronize(stream));

  double VAB = 0.0;
  for (uint32_t i = 0; i < query_atoms_len * target_atoms_len; i++) {
    VAB += h_VABs[i];
  }

  return (VAB / (target_tanimoto + query_tanimoto - VAB));
}