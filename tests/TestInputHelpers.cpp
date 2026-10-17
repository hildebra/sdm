#include <iostream>
#include <stdexcept>
#include <fstream>
#include <cstdio>

#include "InputStream.h"

using namespace std;

int test_ifbufstream_crlf_split_boundary() {
    const char* path = "TestInputHelpers_tmp_crlf_split.txt";
    {
        ofstream out(path, ios::binary);
        out << "AAAAAAAAA\r\nBBBB\r\nCCCC"; // no trailing newline
    }

    ifbufstream in(path, 10, false, false); // force chunk boundary between '\r' and '\n'
    string l1, l2, l3;
    bool ok1 = in.getlines(l1);
    bool ok2 = in.getlines(l2);
    bool ok3 = in.getlines(l3);
    bool ok4 = in.getlines(l3);

    std::remove(path);

    if (!ok1 || !ok2 || !ok3 || ok4) {
        cerr << "test_ifbufstream_crlf_split_boundary: unexpected getline status\n";
        return 1;
    }
    if (l1 != "AAAAAAAAA" || l2 != "BBBB" || l3 != "CCCC") {
        cerr << "test_ifbufstream_crlf_split_boundary: expected AAAAAAAAA/BBBB/CCCC, got '"
             << l1 << "'/'" << l2 << "'/'" << l3 << "'\n";
        return 1;
    }
    return 0;
}

int test_ifbufstream_getlines_stops_on_header_after_crlf() {
    const char* path = "TestInputHelpers_tmp_getlines.txt";
    {
        ofstream out(path, ios::binary);
        out << "seq1\r\nseq2\r\n>next\r\n";
    }

    ifbufstream in(path, 10, false, false);
    string merged;
    bool ok = in.getlines(merged, true, true);
    std::remove(path);

    if (!ok) {
        cerr << "test_ifbufstream_getlines_stops_on_header_after_crlf: expected true\n";
        return 1;
    }
    if (merged != "seq1 seq2 ") {
        cerr << "test_ifbufstream_getlines_stops_on_header_after_crlf: unexpected merged string '" << merged << "'\n";
        return 1;
    }
    return 0;
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
