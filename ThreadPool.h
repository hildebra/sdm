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
    static ThreadPool& gzip_instance();
    static void configure(size_t threads);
    static void configure_gzip(size_t threads);
    static size_t configured_thread_count();
    static size_t configured_gzip_thread_count();
    static void set_gzip_mode(bool enabled);
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    template<class F, class... Args>
    auto submit(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
        return submit_impl_(false, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template<class F, class... Args>
    auto submit_gzip(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
        return gzip_instance().submit_impl_(false, std::forward<F>(f), std::forward<Args>(args)...);
    }

private:
    template<class F, class... Args>
    auto submit_impl_(bool isGzipTask, F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
        using return_type = decltype(f(args...));
        auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            if (stop_.load()) throw std::runtime_error("submit on stopped ThreadPool");
            // wait while queue is full
            while (max_queue_size_ > 0 && (tasks_.size() + gzip_tasks_.size()) >= max_queue_size_) {
                // block until a worker pops a task (and notifies)
                condition_not_full_.wait(lock);
                if (stop_.load()) throw std::runtime_error("submit on stopped ThreadPool");
            }
            // capture submission timestamp for queue wait instrumentation (time of enqueue)
            long long submit_ts = Instr::now_us();
            auto wrappedTask = [task, submit_ts]() {
                // when the worker starts this wrapper, record queue wait
                long long start_ts = Instr::now_us();
                long long wait_us = start_ts - submit_ts;
                Instr::add_pool_queue_wait_us(wait_us);
                // run the actual packaged task
                (*task)();
            };
            if (isGzipTask) {
                gzip_tasks_.emplace(std::move(wrappedTask));
            }
            else {
                tasks_.emplace(std::move(wrappedTask));
            }
            // record sampled queue length (after push)
            Instr::record_pool_queue_len((long long)(tasks_.size() + gzip_tasks_.size()));
        }
        condition_.notify_one();
        return res;
    }

public:
    void shutdown();
    ~ThreadPool();

private:
    ThreadPool(size_t threads = 0);
    static std::atomic<size_t>& configured_threads_();
    static std::atomic<size_t>& configured_gzip_threads_();
    static std::atomic<bool>& instance_created_();
    static std::atomic<bool>& gzip_instance_created_();

    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::queue<std::function<void()>> gzip_tasks_;

    std::mutex queue_mutex_;
    std::condition_variable condition_;
    std::condition_variable condition_not_full_;
    std::atomic<bool> stop_;
    bool gzip_mode_enabled_ = false;
    size_t gzip_reserved_target_ = 0;
    size_t gzip_inflight_ = 0;
    size_t tune_counter_ = 0;
    size_t worker_count_ = 0;
    size_t max_queue_size_ = 0;

    bool has_work_() const {
        return !tasks_.empty() || !gzip_tasks_.empty();
    }
    bool should_take_gzip_() const;
    void maybe_tune_gzip_reservation_();
};
// Inline implementations to make ThreadPool header-only and avoid linker issues
inline std::atomic<size_t>& ThreadPool::configured_threads_() {
    static std::atomic<size_t> configured(0);
    return configured;
}

inline std::atomic<size_t>& ThreadPool::configured_gzip_threads_() {
    static std::atomic<size_t> configured(0);
    return configured;
}

inline void ThreadPool::set_gzip_mode(bool enabled) {
    (void)enabled;
}

inline std::atomic<bool>& ThreadPool::instance_created_() {
    static std::atomic<bool> created(false);
    return created;
}

inline std::atomic<bool>& ThreadPool::gzip_instance_created_() {
    static std::atomic<bool> created(false);
    return created;
}

inline ThreadPool& ThreadPool::instance() {
    static ThreadPool pool(configured_threads_().load());
    instance_created_().store(true);
    return pool;
}

inline ThreadPool& ThreadPool::gzip_instance() {
    static ThreadPool pool(configured_gzip_thread_count());
    gzip_instance_created_().store(true);
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
    if (!gzip_instance_created_().load() && configured_gzip_threads_().load() == 0) {
        size_t gzipThreads = threads / 4;
        if (gzipThreads == 0) {
            gzipThreads = 1;
        }
        configured_gzip_threads_().store(gzipThreads);
    }
}

inline void ThreadPool::configure_gzip(size_t threads) {
    if (threads == 0) {
        return;
    }
    if (gzip_instance_created_().load()) {
        return;
    }
    configured_gzip_threads_().store(threads);
}

inline size_t ThreadPool::configured_thread_count() {
    size_t configured = configured_threads_().load();
    if (configured > 0) {
        return configured;
    }
    size_t hw = std::thread::hardware_concurrency();
    return hw > 0 ? hw : 2;
}

inline size_t ThreadPool::configured_gzip_thread_count() {
    size_t configured = configured_gzip_threads_().load();
    if (configured > 0) {
        return configured;
    }
    size_t compute = configured_thread_count();
    size_t derived = compute / 4;
    if (derived == 0) {
        derived = 1;
    }
    return derived;
}

inline ThreadPool::ThreadPool(size_t threads) : stop_(false) {
    if (threads == 0) {
        threads = std::thread::hardware_concurrency();
        if (threads == 0) threads = 2;
    }
    worker_count_ = threads;
    // set a bounded queue size to provide backpressure; a small multiple of threads
    max_queue_size_ = threads * 4;
    for (size_t i = 0; i < threads; ++i) {
        workers_.emplace_back([this]() {
            for (;;) {
                std::function<void()> task;
                bool tookGzip = false;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex_);
                    this->condition_.wait(lock, [this]() { return this->stop_.load() || this->has_work_(); });
                    if (this->stop_.load() && !this->has_work_()) return;

                    if (should_take_gzip_()) {
                        task = std::move(this->gzip_tasks_.front());
                        this->gzip_tasks_.pop();
                        ++this->gzip_inflight_;
                        tookGzip = true;
                    }
                    else if (!this->tasks_.empty()) {
                        task = std::move(this->tasks_.front());
                        this->tasks_.pop();
                    }
                    else if (!this->gzip_tasks_.empty()) {
                        task = std::move(this->gzip_tasks_.front());
                        this->gzip_tasks_.pop();
                        ++this->gzip_inflight_;
                        tookGzip = true;
                    }
                    else {
                        continue;
                    }
                    // notify one blocked submitter that there is space
                    this->condition_not_full_.notify_one();
                }
                task();
                if (tookGzip) {
                    std::lock_guard<std::mutex> lock(this->queue_mutex_);
                    if (this->gzip_inflight_ > 0) {
                        --this->gzip_inflight_;
                    }
                    this->maybe_tune_gzip_reservation_();
                }
            }
        });
    }
}

inline bool ThreadPool::should_take_gzip_() const {
    if (gzip_tasks_.empty()) {
        return false;
    }
    if (!gzip_mode_enabled_) {
        return true;
    }
    if (tasks_.empty()) {
        return true;
    }
    return gzip_inflight_ < gzip_reserved_target_;
}

inline void ThreadPool::maybe_tune_gzip_reservation_() {
    if (!gzip_mode_enabled_ || worker_count_ < 2) {
        gzip_reserved_target_ = 0;
        return;
    }
    ++tune_counter_;
    if (tune_counter_ < 32) {
        return;
    }
    tune_counter_ = 0;

    const size_t maxReserve = worker_count_ > 2 ? (worker_count_ - 1) : 1;
    const size_t backlog = gzip_tasks_.size() + gzip_inflight_;

    if (backlog > gzip_reserved_target_ * 2 && gzip_reserved_target_ < maxReserve) {
        ++gzip_reserved_target_;
        return;
    }
    if (backlog + 1 < gzip_reserved_target_ && gzip_reserved_target_ > 1) {
        --gzip_reserved_target_;
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

