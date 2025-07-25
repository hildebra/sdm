//
// Created by joachim on 05/06/2020.
//

#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#include <string>
#include <chrono>
#include <iostream>
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
