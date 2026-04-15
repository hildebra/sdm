#pragma once

#include <vector>
#include <thread>
#include <queue>
#include <future>
#include <functional>
#include <condition_variable>
#include <memory>
#include <atomic>
#include <mutex>
#include <stdexcept>
#include "Benchmark.h"


class ThreadPool {
public:
    static ThreadPool& instance();
    static void configure(size_t threads);
    static size_t configured_thread_count();
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    template<class F, class... Args>
    auto submit(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
        using return_type = decltype(f(args...));
        auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            if (stop_.load()) throw std::runtime_error("submit on stopped ThreadPool");
            // wait while queue is full
            while (max_queue_size_ > 0 && tasks_.size() >= max_queue_size_) {
                // block until a worker pops a task (and notifies)
                condition_not_full_.wait(lock);
                if (stop_.load()) throw std::runtime_error("submit on stopped ThreadPool");
            }
            // capture submission timestamp for queue wait instrumentation (time of enqueue)
            long long submit_ts = Instr::now_us();
            tasks_.emplace([task, submit_ts]() {
                // when the worker starts this wrapper, record queue wait
                long long start_ts = Instr::now_us();
                long long wait_us = start_ts - submit_ts;
                Instr::add_pool_queue_wait_us(wait_us);
                // run the actual packaged task
                (*task)();
            });
            // record sampled queue length (after push)
            Instr::record_pool_queue_len((long long)tasks_.size());
        }
        condition_.notify_one();
        return res;
    }

    void shutdown();
    ~ThreadPool();

private:
    ThreadPool(size_t threads = 0);
    static std::atomic<size_t>& configured_threads_();
    static std::atomic<bool>& instance_created_();

    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;

    std::mutex queue_mutex_;
    std::condition_variable condition_;
    std::condition_variable condition_not_full_;
    std::atomic<bool> stop_;
    size_t max_queue_size_ = 0;
};
// Inline implementations to make ThreadPool header-only and avoid linker issues
inline std::atomic<size_t>& ThreadPool::configured_threads_() {
    static std::atomic<size_t> configured(0);
    return configured;
}

inline std::atomic<bool>& ThreadPool::instance_created_() {
    static std::atomic<bool> created(false);
    return created;
}

inline ThreadPool& ThreadPool::instance() {
    static ThreadPool pool(configured_threads_().load());
    instance_created_().store(true);
    return pool;
}

inline void ThreadPool::configure(size_t threads) {
    if (threads == 0) {
        return;
    }
    if (instance_created_().load()) {
        return;
    }
    configured_threads_().store(threads);
}

inline size_t ThreadPool::configured_thread_count() {
    size_t configured = configured_threads_().load();
    if (configured > 0) {
        return configured;
    }
    size_t hw = std::thread::hardware_concurrency();
    return hw > 0 ? hw : 2;
}

inline ThreadPool::ThreadPool(size_t threads) : stop_(false) {
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 2;
    }
    // set a bounded queue size to provide backpressure; a small multiple of threads
    max_queue_size_ = threads * 4;
    for (size_t i = 0; i < threads; ++i) {
        workers_.emplace_back([this]() {
            for (;;) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex_);
                    this->condition_.wait(lock, [this]() { return this->stop_.load() || !this->tasks_.empty(); });
                    if (this->stop_.load() && this->tasks_.empty()) return;
                    task = std::move(this->tasks_.front());
                    this->tasks_.pop();
                    // notify one blocked submitter that there is space
                    this->condition_not_full_.notify_one();
                }
                task();
            }
        });
    }
}

inline ThreadPool::~ThreadPool() {
    shutdown();
}

inline void ThreadPool::shutdown() {
    bool expected = false;
    if (!stop_.compare_exchange_strong(expected, true)) return;
    condition_.notify_all();
    condition_not_full_.notify_all();
    for (auto &worker : workers_) {
        if (worker.joinable()) worker.join();
    }
}

