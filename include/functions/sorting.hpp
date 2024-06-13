#pragma once

#include "uego/uego.h"

#include "functions/pool.hpp"

#include <cinttypes>

class Sorting {
private:
  static void sec_sort(SpeciesList *start, uint32_t len);
  static void inner_parallel_merge_sort(Pool &thread_pool, SpeciesList *start,
                                        uint32_t concurrency, uint32_t len);

public:
  static void parallel_merge_sort(Pool &thread_pool, SpeciesList *root,
                                  uint32_t len);
};
