//
// Created by joachim on 29/07/2020.
//

#include "DNAInputHandler.h"
#include <thread>
#include <mutex>

std::mutex queue_mutex;

void DNAInputHandler::supply2() {
    cout << "start supply" << endl;
    cout << "supply thread: " << std::this_thread::get_id() << endl;
    int count = 0;
    while (has_next) {
        //cout << "read chunk" << endl;
        chunk = new shared_ptr<DNA>[chunk_size_];
        //cout << "init chunk" << endl;
        for (int i = 0; i < chunk_size_; i++) {
            //cout << i << ", ";
            chunk[i] = is_->getDNA(has_next, 0, sync_);
            if (!has_next) {
                read_ = false;
                chunk[i+1] = nullptr;
                return;
            }
        }
        read_ = false;
        
        while (!read_)
            std::this_thread::sleep_for(1ms);
    }
}

void DNAInputHandler::supply() {
    cout << "start supply" << endl;
    cout << "supply thread: " << std::this_thread::get_id() << endl;
    int count = 0;
    while (has_next) {
        count++;
        shared_ptr<DNA> dna_ptr = is_->getDNA(has_next, 0, sync_);
        
        //std::lock_guard<std::mutex> queue_lock(queue_mutex);
        queue_mutex.lock();
        dna_queue_.push(dna_ptr);
        queue_mutex.unlock();
        
        //cout << "too big? " << (dna_queue_.size() >= buffer_) << endl;
        while (dna_queue_.size() >= buffer_) {
            //cout << "wait till queue gets smaller: " << std::this_thread::get_id() << endl;
            std::this_thread::sleep_for(1ms);
        }
    }
}

std::thread DNAInputHandler::runAsThread() {
    return std::thread(&DNAInputHandler::supply2, this);
}

DNAInputHandler::DNAInputHandler(const shared_ptr<InputStreamer> &is) : is_(is) {}

shared_ptr<DNA> DNAInputHandler::getDNAObject() {
    shared_ptr<DNA> dna;
    {
        queue_mutex.lock();
        if (!dna_queue_.empty()) {
            dna = dna_queue_.front();
            dna_queue_.pop();
        }
        queue_mutex.unlock();
    }
    
    return dna;
}


bool DNAInputHandler::hasNext() {
    return has_next;
}

shared_ptr<DNA> *DNAInputHandler::getDNAChunk() {
    while (read_) {
        //cout << ".";
        std::this_thread::sleep_for(0.1ms);
    }
    
    shared_ptr<DNA>* chunk_ptr = chunk;
    read_ = true;
    
    return chunk_ptr;
}
