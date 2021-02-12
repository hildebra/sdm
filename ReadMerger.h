//
// Created by fritsche on 09/12/2020.
//

#pragma once

#include <cstdio>
#include <vector>
#include <string_view>
#include <iostream>
#include <unordered_map>
#include "include/robin_hood.h"
#include "include/robin_map.h"
#include "Benchmark.h"
#include "InputStream.h"


typedef std::string_view SeedStr;
//typedef std::string SeedStr;

typedef robin_hood::unordered_map<SeedStr, int> SeedMap;
//typedef std::unordered_map<SeedStr, int> SeedMap;
//typedef tsl::robin_map<SeedStr, int> SeedMap;

struct Seed {
    size_t pos1 = 0;
    size_t pos2 = 0;
    size_t length = 0;
    bool is2reversed = false;
};

struct MergeResult {
    Seed seed;
	int overlap = -1;
    int offset = -1;
	int offset1 = -1;
	int offset2 = -1;
    double percent_identity = -1.f;
};

struct mergeStats {
	mergeStats():mismatchFwd(0), mismatchRev(0), matchFwd(0), matchRev(0),
		percFwd(0), percRev(0){}
	vector<float> mismatchFwd, mismatchRev, matchFwd, matchRev;
	vector<float> percFwd, percRev;
	float maxF, maxR;
	void prepStats();
	void printLogs();
	void printHisto(string File);
	void logDistri(int p1, int p2, int overlap, bool same);
};
struct qualStats {
	qualStats() :r1(0), r2(0),N1(0),N2(0) {}
	vector<long> r1, r2;
	vector<long> N1,N2;
	void logQuals(vector<qual_score> q1, vector<qual_score> q2);
	void printHisto(string File);
};

class ReadMerger {
private:
    // initial distance from beginning/end of the read to start finding seeds
    const size_t seed_margin_ = 0;

    const size_t seed_dist_ = 3;

    // Number of seeds per read per end (if also search of reverse take this number by 4 to get the total count of extacted seeds.)
    const size_t seed_count_ = 23;

    // Lengh of the seed
    const size_t seed_size_ = 5;

    //Minimum overlap to consider read mergable
    const size_t min_overlap_ = 12; // used to be 16

    // threshold of percent identity necessary to allow merge
    const float percent_identity_threshold_ = 0.9f;

    // check reverse complement
    bool check_reverse_complement_ = true;

    //
    int seed_positions_[30] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,5,7,10,13,17,21,27,34,44,55 };

    // temp for rev transcribed;
    std::string reverse_complement_tmp_ = "";

    // Map for seeds.
//    tsl::robin_map<std::string_view, int> seedmap_;
    SeedMap seedmap_;

//    static std::string generate(char c, size_t times) {
//        if (times == 0) return "";
//        string str = "";
//        str.resize(times);
//        for (auto i = 0; i < times; i++)
//            str[i] = c;
//        return str;
//    }

	static char complement(char c) {
		static int base_map['t' + 1];
		base_map['A'] = 'T';
		base_map['C'] = 'G';
		base_map['G'] = 'C';
		base_map['T'] = 'A';
		base_map['N'] = 'N';

		base_map['a'] = 't';
		base_map['c'] = 'g';
		base_map['g'] = 'c';
		base_map['t'] = 'a';
		base_map['n'] = 'n';

		return base_map[c];
	}

	mergeStats mS;
	qualStats qS;
	bool b_takeStats;



public:
    MergeResult result;

	//functions
	ReadMerger(bool collStats=false):b_takeStats(collStats){}
	~ReadMerger() { cout << "Destroyed Read Merger\n"; }

	void printMergeHisto() { mS.printLogs(); }
	void printQualHisto(string File,string SRb="") { qS.printHisto(File); }
	void printMergeHisto(string File, string SRb = "") { mS.printHisto(File); }
	void restStats() { mS = mergeStats(); qS = qualStats(); }
	void setStatsTaken(bool b) { b_takeStats = b; }

	void reverseStringInPlace(char *str, int len);

	void reverseComplement(char *str, int len);

	std::string reverseComplement(std::string str);

	double percentIdentity(std::string_view sequence1, std::string_view sequence2, int mismatches = -1);

	double percentIdentity(const char* sequence1, const char* sequence2, 
		int length, int mismatches = -1);
   


	/***
     * Find seed in two paired end reads. Store the result in the parameter result.
     * @param sequence1
     * @param sequence2
     * @param result
     * @return
     */
	bool findSeed(std::string sequence1, std::string sequence2, MergeResult &result);

    static double qualToProb(char c) {
        return pow(10, -1 * c/10.f);
    }

    static int probToQual(double prob) {
        return int(log10(prob) * -10 +0.5);
    }

	inline double mergeQProbabilities(double p1, double p2);
	//new qual for mismatches
	inline double mismatchQProbability(double p1, double p2);

    static qual_score inline fastMin(qual_score a, qual_score b) {
        return (a >= b) * b + (b > a) * a;
    }


	qual_score getMergedQual(qual_score q1, qual_score q2,bool same);
	shared_ptr<DNA> merge(shared_ptr<DNA> read1, shared_ptr<DNA> read2);

    bool findSeed(std::string &sequence1, std::string &sequence2) {
        return findSeed(sequence1, sequence2, result);
    }


	void testMergeWithReads(std::istream &is1, std::istream &is2);


    

    static void highlight(std::string str, size_t pos, size_t len = 1, char c = '^', size_t offset = 0) {
        std::cout << std::string(' ', offset) << str << std::endl;
        std::cout << std::string(' ', pos+offset) << std::string(c, len) << std::endl;
    }

    static void highlight(shared_ptr<DNA> dna) {
        std::cout << std::string(dna->merge_offset_, ' ') << dna->getSequence() << std::endl;
        std::cout << std::string(dna->merge_seed_pos_ + dna->merge_offset_, ' ') << std::string(5, '^') << std::endl;
    }
};
