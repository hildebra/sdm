#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
//
// Created by joachim on 13/08/2020.
//

#pragma once

#include "FastxReader.h"
#include "containers.h"
#include "Benchmark.h"
#include "FastxReader.h"
#include <mutex>
#include <thread>
#include "ThreadPool.h"
#include <condition_variable>


static qual_score determineFastqOffset(std::ifstream* is, qual_score fastQver) {
    BufferedFastxReader reader;
    FastxRecord record;
    reader.LoadBatch(*is, (size_t) 20);
    int maxQScore = SCHAR_MIN;
    int minQScore = SCHAR_MAX;
    while (reader.NextSequence(record)) {
        if (record.format == FORMAT_FASTA) return 0;
        for (auto c : record.quality) {
            maxQScore = std::max(maxQScore, (int) c);
            minQScore = std::max(minQScore, (int) c);
        }
    }
    
    int fqDiff(0);
    if (minQScore >= 100 || maxQScore < 2) {
        return fqDiff;
    }
    auto fqSolexaFmt = false;
    if (minQScore >= 59 && maxQScore > 74){
        fqDiff = (fastQver - 64); fastQver = 64;
        if (minQScore < 64) { //set to illumina1.0 (solexa)
            fqSolexaFmt = true;
            cerr << "\nSetting to illumina 1.0-1.3 (solexa) fastq version (q offset = 64, min Q=-5).\n\n";
        } else {
            cerr << "\nSetting to illumina 1.3-1.8 fastq version (q offset = 64).\n\n";
        }
    } else if (minQScore >= 33 && maxQScore <= 74) {
        fqDiff = (fastQver - 33); fastQver = 33;
        cerr << "\nSetting to Sanger fastq version (q offset = 33).\n\n";
    } else {
        cerr << "\nUndecided fastq version..\n";
        fqDiff = (fastQver - 33); fastQver = 0;
        //exit(53);
    }
    is->clear();
    is->seekg(0);
    
    return fqDiff;
}

void setupStreams(shared_ptr<InputStreamer> is, istream*& read1_is, istream*& read2_is, istream*& quality1_is, istream*& quality2_is) {
    if (is->fasta_istreams[0] != nullptr) {
        read1_is = is->fasta_istreams[0];
        cout << "reset" << endl;
        read1_is->seekg(0, ios::beg);

        if (is->quality_istreams[0] != nullptr) {
            quality1_is = is->quality_istreams[0];
            //quality1_is->seekg(0);
        }
        if (is->fasta_istreams[1] != nullptr) {
            read2_is = is->fasta_istreams[1];
            //read2_is->seekg(0);
            if (is->quality_istreams[1] != nullptr) {
                quality2_is = is->quality_istreams[1];
                //quality2_is->seekg(0);
            }
        }
    } else if (is->fastq_istreams[0] != nullptr) {
        read1_is = is->fastq_istreams[0];
        read1_is->seekg(0);
        if (is->fastq_istreams[1] != nullptr) {
            read2_is = is->fastq_istreams[1];
            //read2_is->seekg(0);
        }
    }
}



static std::mutex input_mtx;
void readSingleTpWorker(const shared_ptr<OutputStreamer> &md, const shared_ptr<InputStreamer>& is,
                        uint thread_id, istream *read1_is, istream *read2_is, istream *quality1_is,
                        istream *quality2_is, bool *finished, std::condition_variable &check_quit) {
	shared_ptr<Filters> curFil = md->getFilters();
    bool checkReversedRead = curFil->checkRevRd();

    //Init readers for single_end_read1, single_end_read2 and potential quality files
    BufferedFastxReader reader_read1, reader_read2, reader_qual1, reader_qual2;
    FastxRecord read1, read2, quality1, quality2;

    size_t batchSize = 10000;

    uint64_t processed_reads = 0;

    while (true) {
        bool ok_r1 = false, ok_r2 = false, ok_q1 = false, ok_q2 = false;
        bool valid_fragment_r1 = false, valid_fragment_r2 = false, valid_fragment_q1 = false, valid_fragment_q2 = false;

        {
            // Lock the input area
            std::lock_guard<std::mutex> guard(input_mtx);
            if (read1_is)
                ok_r1 = reader_read1.LoadBatch(*read1_is, batchSize);
            if (read2_is)
                ok_r2 = reader_read2.LoadBatch(*read2_is, batchSize);
            //ok_r2 = reader_read1.LoadBatch(*read2_is, batchSize);
            if (quality1_is) {
                ok_q1 = reader_qual1.LoadBatch(*quality1_is, batchSize);
                reader_qual1.auto_detect();
            }
            if (quality2_is)
                ok_q1 = reader_qual2.LoadBatch(*quality2_is, batchSize);
            //ok_q2 = reader_read1.LoadBatch(*quality2_is, batchSize);
        }

        if (!(ok_r1 || ok_r2 || ok_q1 || ok_q2)) break;

        while (true) {
            valid_fragment_r1 = reader_read1.NextSequence(read1);
            if (quality1_is != nullptr) {
                valid_fragment_q1 = reader_qual1.NextSequence(quality1);
            }

            if (!valid_fragment_r1) break;
            //if (!valid_fragment_q1) break;

            shared_ptr<DNA> dna1, dna2;

            if (quality1_is) {
                is->getDNA(&read1, nullptr, &quality1, nullptr,nullptr, &dna1, &dna2, nullptr);
            }
            else
                is->getDNA(&read1, nullptr, nullptr, nullptr,nullptr, &dna1, &dna2, nullptr);

//            curFil->preFilterSeqStat(dna1, 0);

            curFil->preFilterSeqStatMT(dna1, 0, thread_id);
            curFil->sTotalPlusMT(0);//mutex

            int tagIdx(-2);
            if (checkReversedRead) {
                string presentBC;
                int c_err(0);
                int chkRev(1);
                // Does not alter the Filter object therefore threadsafe
                tagIdx = curFil->findTag(dna1, presentBC, c_err, true, chkRev);

                if (chkRev==0) {//no? undo revTranscr
                    dna1->reverse_transcribe();
                }
            }
            tagIdx = -2;
            md->analyzeDNA(dna1, -1, -1, tagIdx);
            md->dereplicateDNA(dna1, nullptr); //run in extra thread?
            md->write2Demulti(dna1, 0, curFil->getBCoffset());

            // thread_id:: Joachim, there is no saveForWrite that takes thread_id
            if (!md->saveForWrite(dna1)) {
                break;
            }
//            if (dna1->isGreenQual()) {
//                chkDerep++;
//            }
            //

            processed_reads++;
        }
    }
    *finished = true;
    check_quit.notify_all();

//    std::cout << "processed: " << processed_reads << std::endl;
}

void readSingleTp(OptContainer& cmdArgs, shared_ptr<OutputStreamer> md, shared_ptr<InputStreamer> is, ThreadPool *pool) {
	shared_ptr<Filters> curFil = md->getFilters();
    curFil->singReadBC2();
    int chkDerep(0);
    bool checkReversedRead = curFil->checkRevRd();
    bool cont(true); bool sync(false);

    std::istream* read1_is = nullptr;
    std::istream* read2_is = nullptr;
    std::istream* quality1_is = nullptr;
    std::istream* quality2_is = nullptr;

    static const unsigned int batchSize = 10000;
    //static const unsigned int batchSize = 1;
    uint64_t sum = 0;



//    for (int i = 0; i < 3; i++) {
//        cout << "fasta" << i << ": " << (is->fasta_istreams[i] != nullptr) << "    ";
//        cout << "fastq" << i << ": " << (is->fastq_istreams[i] != nullptr) << "    ";
//        cout << "quality" << i << ": " << (is->quality_istreams[i] != nullptr) << endl;
//    }

    setupStreams(is, read1_is, read2_is, quality1_is, quality2_is);

    size_t processing_thread_count = pool->num_threads - 1;
	bool* finished_array = NULL;
	finished_array = new bool[processing_thread_count];
	for (int i = 0; i < processing_thread_count; i++) {		finished_array[i] = false;	}

    std::condition_variable check_quit;


    for (auto i = 0; i < processing_thread_count; i++) {
        bool *finished = &finished_array[i];
        pool->enqueue([md, is, i, read1_is, read2_is, quality1_is, quality2_is, finished, &check_quit] {
            readSingleTpWorker(md, is, i, read1_is, read2_is, quality1_is, quality2_is, finished, check_quit);
        });
    }

    // Has to wait until every thread is finished!
    std::mutex wait_mtx;
    {
        std::unique_lock<std::mutex> lock(wait_mtx);
        check_quit.wait(lock, [&finished_array, processing_thread_count] {
            //for (auto f : finished_array)  if (!f) return false;
			for (size_t x = 0; x < processing_thread_count; x++) { if (!finished_array[x]) return false; }
            return true;
        });
    }
    pool->active = false;

    // automatically merges the filterobjects
    md->closeOutStreams();
}

