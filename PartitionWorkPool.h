#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <span>
#include <stdexcept>
#include <vector>

// A cooperative work pool: existing callers execute ready partitions while
// waiting for their own batch. No extra threads, detached jobs, or nested
// submissions to the processing ThreadPool. One executor per partition.
class PartitionWorkPool {
public:
    struct Work { size_t partition; size_t weight = 1; };
    struct Stats {
        uint64_t batches = 0, tasks = 0, helpedTasks = 0;
        size_t peakActive = 0;
    };
    explicit PartitionWorkPool(size_t partitions, size_t quantum = 256)
        : partitions_(partitions), quantum_(quantum) {
        if (!partitions || !quantum) { throw std::invalid_argument("Empty partition work pool"); }
    }
    PartitionWorkPool(const PartitionWorkPool&) = delete;
    PartitionWorkPool& operator=(const PartitionWorkPool&) = delete;

    // Work and callback remain alive until ALL of this batch has completed,
    // including when a callback throws. Callbacks must not recursively enter
    // this pool or wait for a lifecycle operation which waits for this batch.
    template<class F> void run(std::span<const Work> work, F&& execute) {
        run(work, std::forward<F>(execute), [] {});
    }
    // As run(); afterEnqueue() is called once this batch's jobs are queued.
    // Jobs of one partition execute in queue order, so callers that serialize
    // their enqueues (and release that order in afterEnqueue) get a fixed
    // per-partition execution order across batches.
    template<class F, class G> void run(std::span<const Work> work, F&& execute, G&& afterEnqueue) {
        if (work.empty()) { return; }
        auto callback = [&execute](size_t index) { execute(index); };
        Batch batch{work.size(), &callback, [](void* context, size_t index) {
            (*static_cast<decltype(callback)*>(context))(index);
        }, {}};
        // All potentially throwing allocations/validation precede publication.
        std::vector<Job> jobs(work.size());
        for (size_t i = 0; i < work.size(); ++i) {
            if (work[i].partition >= partitions_.size() || !work[i].weight) {
                throw std::invalid_argument("Invalid partition work item");
            }
            jobs[i] = {&batch, i, work[i].weight, nullptr, {}};
        }
        std::unique_lock<std::mutex> lock(mutex_);
        ++stats_.batches;
        for (size_t i = 0; i < work.size(); ++i) {
            auto& partition = partitions_[work[i].partition];
            if (partition.tail) { partition.tail->next = &jobs[i]; }
            else { partition.head = &jobs[i]; }
            partition.tail = &jobs[i];
            if (!partition.running && !partition.ready) { makeReady(work[i].partition); }
        }
        changed_.notify_all();
        // The jobs point into this frame and are already published: a throwing
        // callback must not unwind before they complete (dangling Job*/Batch*).
        std::exception_ptr enqueueError;
        try { afterEnqueue(); }
        catch (...) { enqueueError = std::current_exception(); }
        while (batch.remaining) {
            changed_.wait(lock, [&] { return !batch.remaining || readyHead_ != none; });
            if (!batch.remaining) { break; }
            const size_t id = readyHead_;
            auto& partition = partitions_[id];
            readyHead_ = partition.nextReady;
            if (readyHead_ == none) { readyTail_ = none; }
            partition.ready = false;
            partition.running = true;
            ++active_;
            if (active_ > stats_.peakActive) { stats_.peakActive = active_; }

            // Coalesce queued batches for this partition to amortize dispatch,
            // then return it to the back of the ready queue for fairness.
            Job* first = partition.head;
            Job* last = first;
            size_t weight = first->weight;
            while (last->next && weight < quantum_ && last->next->weight <= quantum_ - weight) {
                last = last->next;
                weight += last->weight;
            }
            partition.head = last->next;
            if (!partition.head) { partition.tail = nullptr; }
            last->next = nullptr;
            lock.unlock();
            for (Job* job = first; job; job = job->next) {
                try { job->batch->invoke(job->batch->context, job->index); }
                catch (...) { job->error = std::current_exception(); }
            }
            lock.lock();
            for (Job* job = first; job; job = job->next) {
                if (job->error && !job->batch->error) { job->batch->error = job->error; }
                --job->batch->remaining;
                ++stats_.tasks;
                stats_.helpedTasks += job->batch != &batch;
            }
            --active_;
            partition.running = false;
            if (partition.head) { makeReady(id); }
            changed_.notify_all();
        }
        lock.unlock();
        if (enqueueError) { std::rethrow_exception(enqueueError); }
        if (batch.error) { std::rethrow_exception(batch.error); }
    }
    Stats stats() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return stats_;
    }
    // Called only after producers have drained (e.g. under LifecyclePause).
    void resetStats() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (active_ || readyHead_ != none) { throw std::logic_error("Reset of active partition work pool"); }
        stats_ = {};
    }
private:
    static constexpr size_t none = std::numeric_limits<size_t>::max();
    struct Batch {
        size_t remaining;
        void* context;
        void (*invoke)(void*, size_t);
        std::exception_ptr error;
    };
    struct Job {
        Batch* batch;
        size_t index, weight;
        Job* next;
        std::exception_ptr error;
    };
    struct Partition {
        Job* head = nullptr;
        Job* tail = nullptr;
        size_t nextReady = none;
        bool running = false, ready = false;
    };
    void makeReady(size_t id) {
        auto& partition = partitions_[id];
        partition.ready = true;
        partition.nextReady = none;
        if (readyTail_ != none) { partitions_[readyTail_].nextReady = id; }
        else { readyHead_ = id; }
        readyTail_ = id;
    }
    std::vector<Partition> partitions_;
    size_t quantum_, readyHead_ = none, readyTail_ = none, active_ = 0;
    mutable std::mutex mutex_;
    std::condition_variable changed_;
    Stats stats_;
};
