#include "functions/pool.hpp"

Pool::~Pool() {
  {
    std::unique_lock<std::mutex> lock(this->queue_mutex);
    this->should_terminate = true;
  }

  this->mutex_condvar.notify_all();
  for (auto &&thread : this->threads) {
    thread.join();
  }
  threads.clear();
  cout << "threadpool -- Goodbye." << endl;
}

Pool::Pool(uint32_t thread_num, bool pin_threads) {
  this->current_jobs.store(0);
  this->should_terminate = false;

  this->max_threads =
      thread_num == 0 ? std::thread::hardware_concurrency() : thread_num;

  threads.resize(this->max_threads);
  for (uint32_t i = 0; i < this->max_threads; i++) {
    threads.at(i) = std::thread(&Pool::thread_loop, this);

    // TODO Proper affinity when system has multiple NUMA nodes
    // For max. performance, disable SMT on node
    if (pin_threads) {
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(i, &cpuset);

      if (pthread_setaffinity_np(threads.at(i).native_handle(),
                                 sizeof(cpu_set_t), &cpuset) != 0) {
        cout << "threadpool -- Could not set affinity for thread " << i << endl;
      }
    }
  }

  cout << "threadpool -- " << this->max_threads << " work threads available."
       << endl;

  if (pin_threads) {
    cout << "threadpool -- Pool's threads are pinned!" << endl;
  }
}

void Pool::thread_loop() {
  while (true) {
    std::function<void()> job;
    {
      std::unique_lock<std::mutex> lock(this->queue_mutex);
      this->mutex_condvar.wait(lock, [this] {
        return !job_queue.empty() || this->should_terminate;
      });

      if (this->should_terminate) {
        return;
      }

      job = job_queue.front();
      job_queue.pop();
    }

    job();
    this->current_jobs.fetch_sub(1);
  }
}

void Pool::queue_job(const std::function<void()> &job) {
  {
    std::unique_lock<std::mutex> lock(this->queue_mutex);
    this->job_queue.push(job);
  }
  this->current_jobs.fetch_add(1);
  this->mutex_condvar.notify_one();
}

// TODO I dont like the spin-loop
void Pool::wait() {
  while (this->current_jobs.load() > 0)
    ;
}