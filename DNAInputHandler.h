//
// Created by joachim on 29/07/2020.
//
#pragma once

#include <queue>
#include <thread>
#include "InputStream.h"

class DNAInputHandler {
    shared_ptr<InputStreamer> is_;
    shared_ptr<DNA>* chunk;
    
public:
    DNAInputHandler(const shared_ptr<InputStreamer> &is);
    int chunk_size_ = 100000;

private:
    std::queue<shared_ptr<DNA>> dna_queue_;
    int buffer_ = 10000000;
    bool read_ = true;
    bool sync_ = false;
    
    void supply();
    void supply2();
public:
    bool has_next = true;
    shared_ptr<DNA> getDNAObject();
    shared_ptr<DNA>* getDNAChunk();
    std::thread runAsThread();
    bool hasNext();
};
