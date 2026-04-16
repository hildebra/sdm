//
// Created by joachim on 05/06/2020.
//

#pragma once

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include <string>
#include <chrono>
#include <iostream>
#include <array>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include "DNAconsts.h"

using namespace std::chrono;
using namespace std;

class Benchmark {
    string name;
    long event_count = 0;
    long time_sum = 0;
    long dummy_time;
    time_point<high_resolution_clock> start_time;
    time_point<high_resolution_clock> last_time;

public:
    Benchmark(string name) : name(name) {
        auto start_time_dummy = high_resolution_clock::now();

        start_time = high_resolution_clock::now();
        last_time = start_time;
        auto stop_time = high_resolution_clock::now();

        auto stop_time_dummy = high_resolution_clock::now();
        auto duration_dummy = duration_cast<microseconds>(stop_time_dummy - start_time_dummy);

        dummy_time = (long)duration_dummy.count();
    }

    uint now_total_secs() {
        auto now_time = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(now_time - start_time);
        long time_loc = (long)duration.count();
        uint secs = time_loc / 1000000;
        return secs;
    }

    float now_interval_secs(int decimals=0) {
        auto now_time = high_resolution_clock::now();       
        auto duration = duration_cast<microseconds>(now_time - last_time);
        float time_loc = (float)duration.count();
        float secs = time_loc / 1000000.f;
        last_time = now_time;
        //convert decimal places
        if (decimals > 0) {
            const float multiplier = (float)std::pow(10.0, decimals);
            secs = std::ceil(secs * multiplier) / multiplier;
        }
        return secs;
    }

    void start() {
        event_count++;
        start_time = high_resolution_clock::now();
    }

    void stop() {
        auto stop_time = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop_time-start_time);
        time_sum += (long)duration.count();
    }

    void printResults(ostream &os = std::cout) {

        uint hours = time_sum/3600000000;
        time_sum -= hours * 3600000000;
        uint mins = time_sum/60000000;
        time_sum -= mins * 60000000;
        uint secs = time_sum/1000000;
        time_sum -= secs * 1000000;
        uint msecs = time_sum/1000;

        os << name << ": ";
        if (hours) os << hours << "h ";
        if (mins) os << mins << "m ";
        if (secs) os << secs << "s ";
        if (msecs) os << msecs << "ms";
        os << std::endl;
    }
};

// Lightweight instrumentation placed here to avoid adding new files in this workspace run.
namespace Instr {
#if defined(DEBUG) || defined(_DEBUG)
    enum Point {
        GETDNALINES = 0,
        STR2DNA,
        MULTI_READ_READY,
        PROCESS_DNA,
        OUTPUT_WRITE,
        ASYNC_SUBMIT,
        ASYNC_LAG,
        POOL_QUEUE_WAIT,
        POOL_QUEUE_LEN,
        LOCK_WAIT,
        LOCK_HOLD,
        COUNT_POINTS
    };

    static std::array<std::atomic<long long>, COUNT_POINTS> time_us = {};
    static std::array<std::atomic<long long>, COUNT_POINTS> counts = {};

    inline long long now_us() {
        using namespace std::chrono;
        return duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count();
    }

    struct ScopedTimer {
        Point p;
        long long start;
        ScopedTimer(Point pt) : p(pt), start(now_us()) {}
        ~ScopedTimer() {
            long long dur = now_us() - start;
            time_us[p].fetch_add(dur, std::memory_order_relaxed);
            counts[p].fetch_add(1, std::memory_order_relaxed);
        }
    };

    inline void add_us(Point p, long long us) {
        time_us[p].fetch_add(us, std::memory_order_relaxed);
        counts[p].fetch_add(1, std::memory_order_relaxed);
    }

    // Record measured queue wait time (microseconds) for pool tasks
    inline void add_pool_queue_wait_us(long long us) {
        add_us(POOL_QUEUE_WAIT, us);
    }

    // Record a sampled queue length. We store the length as a value in time_us and
    // increment the sample count in counts so avg = total_len / samples.
    inline void record_pool_queue_len(long long len) {
        add_us(POOL_QUEUE_LEN, len);
    }

    inline void printResults(std::ostream &os = std::cerr) {
        static const char* names[COUNT_POINTS] = { "GETDNALINES","STR2DNA","MULTI_READ_READY","PROCESS_DNA","OUTPUT_WRITE","ASYNC_SUBMIT","ASYNC_LAG","POOL_QUEUE_WAIT","POOL_QUEUE_LEN","LOCK_WAIT","LOCK_HOLD" };
        os << "\n--- Instrumentation results (microseconds total / count / avg) ---\n";
        for (int i = 0; i < COUNT_POINTS; ++i) {
            long long t = time_us[i].load(std::memory_order_relaxed);
            long long c = counts[i].load(std::memory_order_relaxed);
            double avg = (c > 0) ? (double)t / (double)c : 0.0;
            os << names[i] << ": total=" << t << " us, count=" << c << ", avg=" << avg << " us\n";
        }
        os << "----------------------------------------------------------\n";
    }

    // Timed exclusive lock for std::mutex and similar
    class TimedLockGuard {
        std::mutex* mtx_ = nullptr;
        long long hold_start_ = 0;
    public:
        TimedLockGuard(std::mutex& m) : mtx_(&m) {
            long long t0 = now_us();
            mtx_->lock();
            long long t1 = now_us();
            add_us(LOCK_WAIT, t1 - t0);
            hold_start_ = now_us();
        }
        ~TimedLockGuard() {
            long long hold = now_us() - hold_start_;
            add_us(LOCK_HOLD, hold);
            mtx_->unlock();
        }
        // non-copyable
        TimedLockGuard(const TimedLockGuard&) = delete;
        TimedLockGuard& operator=(const TimedLockGuard&) = delete;
    };

    // Timed exclusive lock for std::shared_mutex using unique_lock
    class TimedSharedLockGuard {
        std::shared_mutex* smtx_ = nullptr;
        std::unique_lock<std::shared_mutex> lk_;
        long long hold_start_ = 0;
    public:
        TimedSharedLockGuard(std::shared_mutex& m) : smtx_(&m), lk_(*smtx_) {
            // unique_lock already locked by ctor; estimate wait time as negligible here
            // To measure wait time more precisely we'd need try_lock loop; keep simple
            hold_start_ = now_us();
        }
        ~TimedSharedLockGuard() {
            long long hold = now_us() - hold_start_;
            add_us(LOCK_HOLD, hold);
        }
        TimedSharedLockGuard(const TimedSharedLockGuard&) = delete;
        TimedSharedLockGuard& operator=(const TimedSharedLockGuard&) = delete;
    };

#else
    // Minimal no-op instrumentation when DEBUG is not set to avoid runtime overhead and logging.
enum Point {
    GETDNALINES = 0,
    STR2DNA,
    MULTI_READ_READY,
    PROCESS_DNA,
    OUTPUT_WRITE,
    ASYNC_SUBMIT,
    ASYNC_LAG,
    POOL_QUEUE_WAIT,
    POOL_QUEUE_LEN,
    LOCK_WAIT,
    LOCK_HOLD,
    COUNT_POINTS
};

    inline long long now_us() { return 0; }

    struct ScopedTimer { ScopedTimer(Point) {} };

    inline void add_us(Point, long long) {}

    inline void printResults(std::ostream &os = std::cerr) {}

    // No-op helpers for pool instrumentation (present in debug build)
    inline void add_pool_queue_wait_us(long long) {}
    inline void record_pool_queue_len(long long) {}

    // Timed exclusive lock for std::mutex and similar (no timing)
    class TimedLockGuard {
        std::mutex* mtx_ = nullptr;
    public:
        TimedLockGuard(std::mutex& m) : mtx_(&m) { mtx_->lock(); }
        ~TimedLockGuard() { mtx_->unlock(); }
        TimedLockGuard(const TimedLockGuard&) = delete;
        TimedLockGuard& operator=(const TimedLockGuard&) = delete;
    };

    // Timed exclusive lock for std::shared_mutex using unique_lock (no timing)
    class TimedSharedLockGuard {
        std::shared_mutex* smtx_ = nullptr;
        std::unique_lock<std::shared_mutex> lk_;
    public:
        TimedSharedLockGuard(std::shared_mutex& m) : smtx_(&m), lk_(*smtx_) {}
        ~TimedSharedLockGuard() {}
        TimedSharedLockGuard(const TimedSharedLockGuard&) = delete;
        TimedSharedLockGuard& operator=(const TimedSharedLockGuard&) = delete;
    };

#endif
}
