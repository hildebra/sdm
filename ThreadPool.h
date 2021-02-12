//
// Created by fritsche on 04/11/2020.
//

#pragma once

#include <condition_variable>
#include <future>
#include <functional>
#include <vector>
#include <thread>
#include <queue>


class ThreadPool {
public:
    bool active = true;
    using Task = std::function<void()>;

    const size_t num_threads;

    explicit ThreadPool(std::size_t numThreads) : num_threads(numThreads) {
        start(numThreads);
#ifdef DEBUG
		cerr << "Starting multithreading with "<< numThreads << "threads\n";
#endif // DEBUG

    }

    ~ThreadPool() {
        stop();
    }

    void enqueue(Task task) {
        {
            std::unique_lock<std::mutex> {mEventMutex};
            mTasks.emplace(std::move(task));
        }
        mEventVar.notify_one();
    }

private:
    std::vector<std::thread> mThreads;
    std::queue<Task> mTasks;
    std::condition_variable mEventVar;
    std::mutex mEventMutex;
    bool mStopping = false;

    void start(std::size_t numThreads) {
        for (auto i = 0u; i < numThreads; ++i) {
            mThreads.emplace_back([=] {
                while (true) {
                    Task task;
                    {
                        std::unique_lock<std::mutex> lock{mEventMutex};

                        mEventVar.wait(lock, [=] { return mStopping || !mTasks.empty(); });

                        if (mStopping && mTasks.empty()) break;

                        task = std::move(mTasks.front());
                        mTasks.pop();
                    }
                    task();
                }
            });
        }
    }

    void stop() noexcept {
        {
            std::unique_lock<std::mutex> lock{mEventMutex};
            mStopping = true;
        }
        mEventVar.notify_all();

        for (auto &thread : mThreads)
            thread.join();
    }
};