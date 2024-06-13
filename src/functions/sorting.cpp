#include "functions/sorting.hpp"

#include "uego/uegoini.h"

#include <chrono>
#include <thread>

void Sorting::parallel_merge_sort(Pool &thread_pool, SpeciesList *root,
                                  uint32_t len) {
  inner_parallel_merge_sort(thread_pool, root->next,
                            thread_pool.work_threads_count(), len);
}

void Sorting::sec_sort(SpeciesList *start, uint32_t len) {
  SearchSpElement *tmp = INI.Prototype()->RandNew();
  SpeciesList *iter = start;

  for (uint32_t i = 0; i < len; i++) {
    for (uint32_t j = 0; j < len && iter->next != nullptr; j++) {
      double left = iter->center->CurrValue();
      double right = iter->next->center->CurrValue();
      if (left < right) {
        // Sorting
        tmp->UpdateFrom(iter->center);
        iter->center->UpdateFrom(iter->next->center);
        iter->next->center->UpdateFrom(tmp);
      }
      iter = iter->next;
    }
    iter = start;
  }
}

void Sorting::inner_parallel_merge_sort(Pool &thread_pool, SpeciesList *start,
                                        uint32_t concurrency, uint32_t len) {
  if (concurrency == 1) {
    Sorting::sec_sort(start, len);
    return;
  }

  uint32_t split = len / concurrency;
  int32_t rem = len % concurrency;

  SpeciesList *iter = start;
  for (uint32_t id = 0; id < concurrency; id++) {
    uint32_t my_work = split + ((rem-- > 0) ? 1 : 0);

    thread_pool.queue_job(
        [iter, my_work] { Sorting::sec_sort(iter, my_work); });

    for (uint32_t i = 0; i < my_work && iter != nullptr; i++) {
      iter = iter->next;
    }
  }

  thread_pool.wait();

  inner_parallel_merge_sort(thread_pool, start, concurrency / 2, len);
}
