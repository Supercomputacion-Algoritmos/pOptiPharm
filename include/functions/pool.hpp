#pragma once

#include "uego/searchsp.h"
#include "uego/speclist.h"

#include <atomic>
#include <cinttypes>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

// Implementation mainly from https://stackoverflow.com/a/32593825
class Pool {
private:
  // Threading
  uint32_t max_threads;
  std::vector<std::thread> threads;

  // Queue
  std::mutex queue_mutex;
  std::condition_variable mutex_condvar;
  std::atomic_uint32_t current_jobs;
  std::queue<std::function<void()>> job_queue;

  // Other
  bool should_terminate;

  void thread_loop();

public:
  Pool(uint32_t thread_num = 0, bool pin_threads = false);
  Pool(Pool &) = delete;
  ~Pool();

  uint32_t work_threads_count() { return this->max_threads; }

  void queue_job(const std::function<void()> &job);
  void wait();
};
