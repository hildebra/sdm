#pragma once
// Shared rough conserved/variable-block simulator for dereplication experiments.
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <string_view>
#include <vector>
class CoarseSimulation {
public:
    CoarseSimulation(size_t meanLength, unsigned quality, size_t jitter = 0)
        : jitter(jitter), errorRate(std::pow(10., -double(quality) / 10.)) {
    const size_t length = meanLength + jitter;
    std::string ancestral(length, 'A');
    for (char& base : ancestral) { base = "ACGT"[templateRng() % 4]; }
    auto variable = [meanLength](size_t position) {
        position = position * 1500 / meanLength;
        for (size_t block = 0; block < 9; ++block) {
            if (position >= 80 + block * 155 && position < 140 + block * 155) { return true; }
        }
        return false;
    };
    for (size_t lineage = 0; lineage < 5; ++lineage) {
        auto parent = ancestral;
        for (size_t i = 0; i < length; ++i) {
            if (variable(i) && templateRng() % 10000 < 2000) { randomMutation(parent[i], templateRng); }
        }
        for (size_t variant = 0; variant < 10; ++variant) {
            auto sequence = parent;
            for (size_t i = 0; i < length; ++i) {
                if (variable(i) && templateRng() % 10000 < 200) { randomMutation(sequence[i], templateRng); }
            }
            templates.push_back(std::move(sequence));
        }
    }
    }
    std::string next() {
        lastTemplate = orderRng() % templates.size();
        auto sequence = templates[lastTemplate];
        if (jitter) { sequence.resize(sequence.size() - orderRng() % (2 * jitter + 1)); }
        for (char& base : sequence) {
            if (double(readRng() % 10000) < errorRate * 10000) { randomMutation(base, readRng); ++errors; }
        }
        bases += sequence.size(); return sequence;
    }
    uint64_t errors = 0, bases = 0;
    size_t lastTemplate = 0;
private:
    static void randomMutation(char& base, std::mt19937_64& rng) {
        const size_t code = std::string_view("ACGT").find(base);
        base = "ACGT"[(code + 1 + rng() % 3) % 4];
    }
    size_t jitter;
    double errorRate;
    std::mt19937_64 templateRng{2026091201ULL}, readRng{2026091202ULL}, orderRng{2026091220ULL};
    std::vector<std::string> templates;
};
