#pragma once

#include <memory>
#include <string>
#include <vector>

// HTSlib stays behind this interface so the portable build has no HTS dependency.
// Records use the existing FASTQ batch representation; optional element 4 is BC.
class HtsInput {
public:
    HtsInput(const std::string& filename, int pairs, const std::string& onlyPair,
             const std::string& reference, bool barcodeTags);
    ~HtsInput();
    HtsInput(const HtsInput&) = delete;
    HtsInput& operator=(const HtsInput&) = delete;

    static std::string detectFormat(const std::string& filename);
    static std::string version();
    bool read(std::vector<std::string>& record, int position);
    bool hasData() const;
    int pairCount() const;
    void reset();
    const std::string& filename() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
