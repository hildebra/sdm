// SPDX-License-Identifier: GPL-3.0-or-later
#include "PartitionWorkPool.h"
#include <array>
#include <atomic>
#include <chrono>
#include <iostream>
#include <thread>

static void check(bool ok, const char* message) {
    if (!ok) { throw std::runtime_error(message); }
}
static void testHelping() {
    PartitionWorkPool pool(2, 1);
    std::mutex mutex;
    std::condition_variable changed;
    bool blocked = false, helped = false;
    std::array<std::atomic<int>, 2> active{};
    std::array<std::exception_ptr, 2> errors{};
    std::thread::id helper;
    std::jthread first([&] {
        try {
            const std::array<PartitionWorkPool::Work, 2> work{{{0}, {1}}};
            pool.run(work, [&](size_t i) {
                check(active[i].fetch_add(1) == 0, "Concurrent updates to one partition");
                if (i == 0) {
                    std::unique_lock lock(mutex);
                    blocked = true; changed.notify_all();
                    check(changed.wait_for(lock, std::chrono::seconds(5), [&] { return helped; }),
                        "Ready partition did not make progress while another was blocked");
                } else {
                    std::lock_guard lock(mutex);
                    helper = std::this_thread::get_id();
                    helped = true; changed.notify_all();
                }
                --active[i];
            });
        } catch (...) { errors[0] = std::current_exception(); }
    });
    {
        std::unique_lock lock(mutex);
        check(changed.wait_for(lock, std::chrono::seconds(5), [&] { return blocked; }), "First worker did not start");
    }
    std::jthread second([&] {
        try {
            const std::array<PartitionWorkPool::Work, 1> work{{{0}}};
            pool.run(work, [&](size_t) {
                check(active[0].fetch_add(1) == 0, "Busy partition was claimed twice");
                --active[0];
            });
        } catch (...) { errors[1] = std::current_exception(); }
    });
    const auto secondId = second.get_id();
    first.join(); second.join();
    for (const auto& error : errors) { if (error) { std::rethrow_exception(error); } }
    check(helper == secondId, "Worker waiting for partition zero must execute another batch's ready partition");
    const auto stats = pool.stats();
    check(stats.tasks == 3 && stats.batches == 2 && stats.helpedTasks >= 1 && stats.peakActive == 2,
        "Pool records completed work and simultaneous partition owners");
}
static void testStressAndExceptions() {
    constexpr size_t threads = 8, partitions = 16, batches = 128, tasks = 32;
    PartitionWorkPool pool(partitions, 16);
    std::array<std::atomic<unsigned>, partitions> active{};
    std::array<std::atomic<uint64_t>, partitions> counts{};
    std::array<uint64_t, partitions> expected{};
    std::array<std::exception_ptr, threads> errors{};
    auto partition = [](size_t t, size_t b, size_t i) { return i % 2 ? 0 : (t + b + i) % partitions; };
    for (size_t t = 0; t < threads; ++t) for (size_t b = 0; b < batches; ++b)
        for (size_t i = 0; i < tasks; ++i) { ++expected[partition(t, b, i)]; }
    std::vector<std::jthread> workers;
    for (size_t t = 0; t < threads; ++t) {
        workers.emplace_back([&, t] {
            try {
                for (size_t b = 0; b < batches; ++b) {
                    std::array<PartitionWorkPool::Work, tasks> work;
                    for (size_t i = 0; i < tasks; ++i) { work[i] = {partition(t, b, i), 1 + i % 7}; }
                    pool.run(work, [&](size_t i) {
                        const size_t id = work[i].partition;
                        check(active[id].fetch_add(1) == 0, "Partition exclusivity violated under contention");
                        ++counts[id];
                        --active[id];
                    });
                }
            } catch (...) { errors[t] = std::current_exception(); }
        });
    }
    workers.clear(); // jthread joins; every callback must have finished.
    for (const auto& error : errors) { if (error) { std::rethrow_exception(error); } }
    for (size_t i = 0; i < partitions; ++i) { check(counts[i] == expected[i] && active[i] == 0, "Lost or duplicated work"); }
    check(pool.stats().tasks == threads * batches * tasks, "Incomplete work was reported as completed");
    unsigned finished = 0;
    const std::array<PartitionWorkPool::Work, 3> work{{{0}, {1}, {0}}};
    bool caught = false;
    try {
        pool.run(work, [&](size_t i) { ++finished; if (i == 0) { throw std::runtime_error("sentinel"); } });
    } catch (const std::runtime_error& error) { caught = std::string(error.what()) == "sentinel"; }
    check(caught && finished == 3, "Callback failure must drain the batch before propagation");
    pool.run(work, [&](size_t) { ++finished; });
    check(finished == 6, "Pool could not be reused after a callback failure");
    const std::array<PartitionWorkPool::Work, 2> invalid{{{0}, {partitions}}};
    caught = false;
    try { pool.run(invalid, [&](size_t) { ++finished; }); }
    catch (const std::invalid_argument&) { caught = true; }
    check(caught && finished == 6, "Invalid batch was partially published");
    pool.resetStats();
    check(pool.stats().tasks == 0 && pool.stats().peakActive == 0, "Pool statistics did not reset");
}
int main() {
    try { testHelping(); testStressAndExceptions(); }
    catch (const std::exception& error) { std::cerr << "FAIL: " << error.what() << '\n'; return 1; }
    std::cout << "Partition work pool tests passed\n";
}
