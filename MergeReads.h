//
// Created by fritsche on 02/12/2020.
//

#pragma once


#include <string>
#include "FastxReader.h"
#include "containers.h"

namespace Merge {

    enum Orientation {
        AUTO_DETECT,  OFFSET_FIRST, OFFSET_SECOND
    };

    struct MergeParams {
        Orientation orientation = AUTO_DETECT;
        int seed_margin = 5;
        size_t seed_dist = 10;
        size_t seed_size = 10;
        size_t seed_count = 2;
        size_t max_search = 150;
        size_t mismatches = 0;
        std::vector<size_t> offsets;
    };

    struct MergeResult {
        int anchor1_seq1;
        int anchor1_seq2;
        int anchor2_seq1;
        int anchor2_seq2;
        int overlap;
        double percent_identity;
    };

    static std::string generate(char c, size_t times) {
        if (times == 0) return "";
        string str = "";
        str.resize(times);
        for (auto i = 0; i < times; i++)
            str[i] = c;
        return str;
    }

    // utility functions
    static void highlight(std::string str, size_t pos, size_t len = 1, char c = '^', size_t offset = 0) {
        std::cout << generate(' ', offset) << str << std::endl;
        std::cout << generate(' ', pos+offset) << generate(c, len) << std::endl;
    }

    // utility functions
    static void highlight(std::string str, std::vector<size_t> &pos, size_t len = 1, char c = '^', size_t offset = 0) {
        std::cout << generate(' ', offset) << str << std::endl;

        string highlight_str(str.length() + offset, ' ');
        for (auto hpos : pos) {
            for (auto i = 0; i < len; i++)
                highlight_str[hpos + offset + i] = c;
        }
        std::cout << highlight_str << std::endl;
    }

    static void reverseStringInPlace(char *str, int len) {
        char *p1 = str;
        char *p2 = str + len - 1;

        while (p1 < p2) {
            char tmp = *p1;
            *p1++ = *p2;
            *p2-- = tmp;
        }
    }

    static char complement(char c) {
        static int base_map['t'+1];
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

    static void reverseComplement(char *str, int len) {
        char *p1 = str;
        char *p2 = str + len - 1;

        while (p1 < p2) {
            char tmp = complement(*p1);
            *p1++ = complement(*p2);
            *p2-- = tmp;
        }
    }

    static double percentIdentity (std::string_view sequence1, std::string_view sequence2, bool search_forward=false) {
        int match_count = 0;
        int base_count = 0;

//        std::cout << sequence1 << std::endl;
//        std::cout << sequence2 << std::endl;

        if (sequence1.length() != sequence2.length()) {
            std::cerr << "Sequences must be of the same length." << std::endl;
            return -1.f;
        }

        if (search_forward) {
            for (auto i = 0; i < sequence1.length(); ++i) {
                match_count += sequence1[i] == sequence2[i];
                ++base_count;
            }
        } else {
            for (auto i = sequence1.length()-1; i != 0; --i) {
                match_count += sequence1[i] == sequence2[i];
                ++base_count;
            }
        }

        return (double) match_count / base_count;
    }

    static int getAnkerDistDifference(int anchor1_s1, int anchor2_s1, int anchor1_s2, int anchor2_s2) {
        int diff_seq1 = anchor1_s1 - anchor2_s1;
        int diff_seq2 = anchor1_s2 - anchor2_s2;
        if (diff_seq1 < 0 && diff_seq2 > 0 || diff_seq1 > 0 && diff_seq2 < 0)
            return -1;

        return abs(diff_seq2 - diff_seq1);
    }


    static inline int findSeed(std::string& sequence, std::string_view& seed, int start = 0, int mismatches = 0) {
        int seed_pos;
        bool found = false;

//        std::cout << "findseed: " << seed << std::endl;

        for (auto seq_pos = start; seq_pos < sequence.length() - seed.length() + 1; seq_pos++) {
            seed_pos = 0;
            auto mismatch_count = 0;
            while (true) {
                if (sequence[seq_pos + seed_pos] != seed[seed_pos] && ++mismatch_count > mismatches) {
//                    std::cout << "abort " << seq_pos << " "  << std::endl;
                    break;
                }
                if (++seed_pos == seed.length()) return seq_pos;
            }
        }
        return -1;
    }

    static inline int findSeed2(std::string& sequence, std::string_view& seed, int start_margin = 5, int mismatches = 0, bool search_forward=true) {
        int seed_pos;
        int seq_pos;

        size_t mismatch_count;
        if (search_forward) {
            for (seq_pos = start_margin; seq_pos < sequence.length(); ++seq_pos) {
//                std::cout << seq_pos << " ";
                seed_pos = -1;
                while (seed.length() != ++seed_pos && sequence[seq_pos + seed_pos] == seed[seed_pos]) ;
//                std::cout << seed_pos << " (" << seed.length() << ") seq_pos: " << seq_pos << " ";
                if (seed_pos == seed.length()) {
//                    std::cout << " return: " << seq_pos << std::endl;
                    return seq_pos;
                }
            }
        } else {
            for (seq_pos = sequence.length() - start_margin - 1; seq_pos != 0; --seq_pos) {
                seed_pos = -1;
                while (seed.length() != ++seed_pos && sequence[seq_pos + seed_pos] == seed[seed_pos]) ;
                if (seed_pos == seed.length()) {
                    return seq_pos;
                }
            }
        }

//        std::cout << " return: " << -1 << std::endl;
        return -1;
    }

    static inline int findSeed3(std::string& sequence, std::string_view& seed, int start_margin = 5, int mismatches = 0, bool search_forward=true) {
        int seed_pos;
        int seq_pos;

        size_t mismatch_count;
        if (search_forward) {
            for (seq_pos = start_margin; seq_pos < sequence.length(); ++seq_pos) {
//                std::cout << seq_pos << " ";
                seed_pos = -1;
                while (seed.length() != ++seed_pos && sequence[seq_pos + seed_pos] == seed[seed_pos]) ;
//                std::cout << seed_pos << " (" << seed.length() << ") seq_pos: " << seq_pos << " ";
                if (seed_pos == seed.length()) {
//                    std::cout << " return: " << seq_pos << std::endl;
                    return seq_pos;
                }
            }
        } else {
            for (seq_pos = sequence.length() - start_margin - 1; seq_pos != 0; --seq_pos) {
                seed_pos = -1;
                while (seed.length() != ++seed_pos && sequence[seq_pos + seed_pos] == seed[seed_pos]) ;
                if (seed_pos == seed.length()) {
                    return seq_pos;
                }
            }
        }

//        std::cout << " return: " << -1 << std::endl;
        return -1;
    }

    static double merge2(std::string& sequence1, std::string& sequence2, MergeResult &stats, MergeParams &params) {
        switch (params.orientation) {
            case OFFSET_FIRST: {
                std::string_view seed_s1a(sequence1.c_str() + params.seed_margin, params.seed_size);
                std::string_view seed_s1b(sequence1.c_str() + params.seed_margin, params.seed_size);
                break;
            }
            case OFFSET_SECOND: {
                std::string_view seed_s2a(sequence2.c_str() + params.seed_margin, params.seed_size);
                std::string_view seed_s2b(sequence2.c_str() + params.seed_margin, params.seed_size);
                break;
            }
            case AUTO_DETECT: {
                std::string_view seed_s1a(sequence1.c_str() + params.seed_margin, params.seed_size);
                std::string_view seed_s1b(sequence1.c_str() + params.seed_margin, params.seed_size);
                std::string_view seed_s2a(sequence2.c_str() + params.seed_margin, params.seed_size);
                std::string_view seed_s2b(sequence2.c_str() + params.seed_margin, params.seed_size);


                for (int seq_pos = sequence2.length() - params.seed_size; seq_pos >= 0; --seq_pos) {
                    auto seed_pos = -1;
                    auto mismatch_count = 0;
                    while ((mismatch_count -= sequence2[seq_pos] != seed_s1a[++seed_pos]) != params.mismatches);
                    if (seed_pos == params.seed_size) {

                    }
                }
                for (int seq_pos = sequence2.length(); seq_pos >= params.seed_margin; --seq_pos) {
                    auto seed_pos = 0;
                    auto mismatch_count = 0;
                    while ((mismatch_count -= !(sequence2[seq_pos])) != params.mismatches ) {

                    }
                }

                break;
            }
        }

        if (params.orientation == OFFSET_SECOND) {


        }
        if (params.orientation == OFFSET_FIRST) {
        }


    }

    static double merge(std::string& sequence1, std::string& quality1, std::string& sequence2, std::string& quality2) {
        const int min_seed_num = 2;
        const int min_overlap = 30;
        const int seed_length = 10;
        const int mismatches = 2;
        const int margin = 19;
        const int tries_per_side = 20;
        const int try_multiplier = 1;

        // seed positions are relative to sequence1
        int seed1_pos = sequence1.length() - seed_length - margin;
        std::string_view seed1(sequence1.c_str() + seed1_pos, seed_length);
        auto seed1_seq2 = findSeed(sequence2, seed1);

        std::vector<std::size_t> anchors1;
        std::vector<std::size_t> anchors2;

        int anchor1_s1;
        int anchor2_s1;
        int anchor1_s2;
        int anchor2_s2;
        int offset_s1 = 0;
        int offset_s2 = 0;
        double percent_identity;

        // If true sequence1 is offset relativ to s2 and not the other way around
        // s1             GGTGTGACGATCGATCGG
        // s2 ACGTGACTGCATGGTGTGAC

        bool s1_is_offset = false;

//        anchor1_s1 = margin;
        anchor1_s1 = 150;
//        std::cout << "start: " << anchor1_s1 << std::endl;
//        std::cout << "len: " << sequence1.length() << std::endl;

        // FIND FIRST ANCHOR
        for (auto pos = 0; ; pos++) {
//            std::cout << anchor1_s1 << std::endl;
            if (anchor1_s1 < 0 || anchor1_s1 > sequence1.length() - seed_length) break;

            std::string_view seed(sequence1.c_str() + anchor1_s1, seed_length);
//            anchor1_s2 = findSeed(sequence2, seed, 0, mismatches);
            anchor1_s2 = findSeed2(sequence2, seed, 5, 0, false);

//            std::cout << anchor1_s1 << " " << anchor1_s2 << " " << seed << std::endl;

            if (anchor1_s2 != -1) {
//                std::cout << "found ";
                if (anchor1_s1 > anchor1_s2) {
                    s1_is_offset = false;
                    offset_s2 = anchor1_s1 - anchor1_s2;
//                    std::cout << " s2 offset " << anchor1_s2 << " " << anchor1_s1;
                } else {
//                    std::cout << " s1 offset " << anchor1_s2 << " " << anchor1_s1;
                    s1_is_offset = true;
                    offset_s1 = anchor1_s2 - anchor1_s1;
                }
//                std::cout << std::endl;
                anchors1.push_back(anchor1_s1);
                anchors2.push_back(anchor1_s2);
                break;
            }

            anchor1_s1 += try_multiplier;
        }

        if (anchor1_s2 == -1) {
//            std::cerr << "could not find first anchor" << std::endl;
            return 0.f;
        }

        // extend seed


        auto overlap = s1_is_offset ? sequence2.length() - offset_s1 : sequence1.length() - offset_s2;
        percent_identity = percentIdentity(std::string_view(sequence1.c_str() + offset_s2, overlap), std::string_view(sequence2.c_str() + offset_s1, overlap));
//        percent_identity = 1.;

        // FIND SECOND ANCHOR

////        std::cout << "s1 is offset: " << s1_is_offset << std::endl;
//        anchor2_s2 = s1_is_offset ? sequence2.length() - seed_length : 0;
////        std::cout << "starting anchor2: " << anchor2_s1 << std::endl;
//        while (true) {
//            if (anchor2_s2 < 0 || anchor2_s2 > sequence2.length() - seed_length + 1)
//                break;
//
//            std::string_view seed(sequence2.c_str() + anchor2_s2, seed_length);
////            std::cout << anchor2_s2 << " " << seed << std::endl;
//            anchor2_s1 = findSeed(sequence1, seed, 0, mismatches);
//
//            if (anchor2_s1 != -1) {
////                std::cout << anchor2_s1 << std::endl;
//                anchors1.push_back(anchor2_s1);
//                anchors2.push_back(anchor2_s2);
//                break;
//            }
//
//            if (s1_is_offset)
//                anchor2_s2--;
//            else
//                anchor2_s2++;
//        }

//        if (anchor2_s2 == -1) {
////            std::cout << s1_is_offset << std::endl;
//            highlight(sequence1, anchor1_s1, seed_length, '^', offset_s1);
//            highlight(sequence2, anchor1_s2, seed_length, '^', offset_s2);
//            std::cerr << "could not find second anchor" << std::endl;
//            return false;
//        }
//        highlight(sequence1, anchor2_s1, seed_length);
//        highlight(sequence2, anchor2_s2, seed_length, '^', anchor2_s1 - anchor2_s2);

        size_t min_dist = 40;

//        int dist = getAnkerDistDifference(anchor1_s1, anchor2_s1, anchor1_s2, anchor2_s2);
//        cout << "dist: " << dist << std::endl;

//        if (percent_identity > 0.95) {
//
//            std::cout << "anchor1:" << anchor1_s1 << ", " << anchor1_s2 << std::endl;
//            std::cout << "anchor2:" << anchor2_s1 << ", " << anchor2_s2 << std::endl;
//            highlight(sequence1, anchors1, seed_length, '^', offset_s1);
//            highlight(sequence2, anchors2, seed_length, '^', offset_s2);
//            std::cout << "overlap: " << overlap << std::endl;
//            std::cout << "percent identity: " << percent_identity << std::endl;
//        }
        return percent_identity;
    }



    static void testMerge() {
        std::string sequence = "ACGTACGTCAGTCATGTCGATCGTAGTAGCTAGCATGTACTACGTAGCTGACTGACGACGTACTACACGGTGTGTAGCTGCTGACTG";
        std::string sequence2 = "TAGCATGTACTACGTAGCTGACTGACGACGTACTACACGGTGTGTAGCTGCTGACTGGCTGTTGACTGATCGTCGTATTATATGTA";
        std::string qual1 = "";
        std::string qual2 = "";

        int seed_length = 10;
        int seed_pos = 0;
        std::string_view seed(sequence.c_str() + seed_pos, seed_length);

        for (seed_pos = 0; seed_pos < sequence.length() - seed.length() + 1; seed_pos++) {
            std::string_view seed(sequence.c_str() + seed_pos, seed_length);
            auto result = findSeed(sequence, seed);
            std::cout << (string)seed << ": " << result << std::endl;
        }
        std::cout << sequence << std::endl;
        std::cout << "0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890" << std::endl;


        seed = "ATGGCGCC";
        std::cout << (string)seed << ": " << findSeed(sequence, seed) << std::endl;

        seed = "ACTG";

        int result = -1;
        while ((result = findSeed(sequence, seed, ++result)) >= 0)
            std::cout << (string) seed << ": " << result << std::endl;


        std::cout << "\nTEST MERGE" << std::endl;

        merge(sequence, qual1, sequence2, qual2);
    }

    void testMergeWithReads(std::istream &is1, std::istream &is2) {
        Benchmark merge_bm("Merge");
        merge_bm.start();
        BufferedFastxReader r1;
        BufferedFastxReader r2;

        FastxRecord record1;
        FastxRecord record2;

        const size_t batch_size = 100000;

        bool check_r1 = false;
        bool check_r2 = false;

        size_t counter = 0;
        size_t merge_counter = 0;

        size_t greater95 = 0;
        size_t smaller95 = 0;
        size_t zero = 0;

        while (true) {
            check_r1 = r1.LoadBatch(is1, batch_size);
            check_r2 = r2.LoadBatch(is2, batch_size);
            if (!check_r1 || !check_r2) break;

            while (r1.NextSequence(record1) && r2.NextSequence(record2)) {
                ++counter;
//                std::cout << record1.header << std::endl;
                //std::cout << record2.to_string();

                //rev second
                reverseComplement(const_cast<char *>(record2.sequence.c_str()), record2.sequence.length());

                auto percent_identity = merge(record1.sequence, record1.quality, record2.sequence, record2.quality);


                if (percent_identity >= 0.90) {
                    greater95++;
                } else if (percent_identity == 0) {
                    zero++;
                } else {
                    smaller95++;
                }
//                if (success) {
//                    std::string stop;
//                    std::cin >> stop;
//                }
            }
        }
        merge_bm.stop();
        std::cout << counter << " reads"  << std::endl;
        std::cout << greater95 << " > 95% id reads"  << std::endl;
        std::cout << smaller95 << " < 95% id reads"  << std::endl;
        std::cout << zero << " no seed found"  << std::endl;
        std::cout << counter << " reads"  << std::endl;
        std::cout << (double)greater95/counter << " merged"  << std::endl;

        merge_bm.printResults();
    }

}