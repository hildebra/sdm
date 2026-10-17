#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <fstream>
#include <cstdio>

#include "FastxReader.h"
#include "InputStream.h"

using namespace std;

int test_empty_input() {
    stringstream ss("");
    FastxRecord rec;
    string buf;
    bool ok = BufferedFastxReader::ReadNextSequence(ss, rec, buf);
    if (ok) {
        cerr << "test_empty_input: expected false, got true\n";
        return 1;
    }

int test_ifbufstream_crlf_split_boundary() {
    const char* path = "TestFastxReader_tmp_crlf_split.txt";
    {
        ofstream out(path, ios::binary);
        out << "AAAA\r\nBBBB\r\nCCCC"; // no trailing newline
    }

    ifbufstream in(path, 5, false, false); // force chunk boundary between '\r' and '\n'
    string l1, l2, l3;
    bool ok1 = in.getline(l1);
    bool ok2 = in.getline(l2);
    bool ok3 = in.getline(l3);
    bool ok4 = in.getline(l3);

    std::remove(path);

    if (!ok1 || !ok2 || !ok3 || ok4) {
        cerr << "test_ifbufstream_crlf_split_boundary: unexpected getline status\n";
        return 1;
    }
    if (l1 != "AAAA" || l2 != "BBBB" || l3 != "CCCC") {
        cerr << "test_ifbufstream_crlf_split_boundary: expected AAAA/BBBB/CCCC, got '"
             << l1 << "'/'" << l2 << "'/'" << l3 << "'\n";
        return 1;
    }
    return 0;
}

int test_ifbufstream_getlines_stops_on_header_after_crlf() {
    const char* path = "TestFastxReader_tmp_getlines.txt";
    {
        ofstream out(path, ios::binary);
        out << "seq1\r\nseq2\r\n>next\r\n";
    }

    ifbufstream in(path, 5, false, false);
    string merged;
    int linesRead = 0;
    bool ok = in.getlines(merged, linesRead, true);
    std::remove(path);

    if (!ok) {
        cerr << "test_ifbufstream_getlines_stops_on_header_after_crlf: expected true\n";
        return 1;
    }
    if (linesRead != 2) {
        cerr << "test_ifbufstream_getlines_stops_on_header_after_crlf: expected 2 lines, got " << linesRead << "\n";
        return 1;
    }
    if (merged != "seq1 seq2 ") {
        cerr << "test_ifbufstream_getlines_stops_on_header_after_crlf: unexpected merged string '" << merged << "'\n";
        return 1;
    }
    return 0;
}
    return 0;
}

int test_malformed_header() {
    stringstream ss(">not_fastq\nFRAME\n+\n!!!!\n");
    FastxRecord rec;
    string buf;
    try {
        // auto-detect will see '>' and treat as FASTA; but header starting with '>' is OK
        // To trigger malformed FASTQ, provide a header that is not '@' but force FORMAT_FASTQ
        bool ok = BufferedFastxReader::ReadNextSequence(ss, rec, buf, FORMAT_FASTQ);
        (void)ok;
    } catch (const runtime_error &e) {
        return 0; // expected
    } catch (...) {
        cerr << "test_malformed_header: unexpected exception type\n";
        return 1;
    }
    cerr << "test_malformed_header: expected exception, none thrown\n";
    return 1;
}

int test_fastq_mismatched_qual() {
    // FASTQ with sequence length 4 and quality length 3 should trigger validation
    stringstream ss("@id\nACGT\n+\n!!!\n");
    FastxRecord rec;
    string buf;
    try {
        bool ok = BufferedFastxReader::ReadNextSequence(ss, rec, buf, FORMAT_FASTQ);
        (void)ok;
    } catch (const runtime_error &e) {
        // expected
        return 0;
    } catch (...) {
        cerr << "test_fastq_mismatched_qual: unexpected exception type\n";
        return 1;
    }
    cerr << "test_fastq_mismatched_qual: expected exception, none thrown\n";
    return 1;
}

int test_remove_paired_info_short() {
    string s = "A";
    remove_paired_info(s, 0);
    if (s != "A") {
        cerr << "test_remove_paired_info_short: expected 'A', got '" << s << "'\n";
        return 1;
    }
    return 0;
}

int test_header_stem_bad() {
    string h = "badheader";
    try {
        string out = header_stem(h);
        (void)out;
    } catch (const runtime_error &e) {
        return 0;
    } catch (...) {
        cerr << "test_header_stem_bad: unexpected exception type\n";
        return 1;
    }
    cerr << "test_header_stem_bad: expected exception, none thrown\n";
    return 1;
}

int main() {
    int failures = 0;
    failures += test_empty_input();
    failures += test_malformed_header();
    failures += test_fastq_mismatched_qual();
    failures += test_remove_paired_info_short();
    failures += test_header_stem_bad();
    failures += test_ifbufstream_crlf_split_boundary();
    failures += test_ifbufstream_getlines_stops_on_header_after_crlf();

    if (failures == 0) {
        cout << "All tests passed\n";
        return 0;
    } else {
        cerr << failures << " test(s) failed\n";
        return 1;
    }
}
