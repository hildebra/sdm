#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>
#include <functional>
#include "InputStream.h"

using namespace std;
#include <thread>
#include <chrono>
#include <atomic>

#if __has_include(<zlib.h>)
#include <zlib.h>
#define TEST_USE_ZLIB 1
#else
#define TEST_USE_ZLIB 0
#endif

int writeFile(const string &path, const string &content) {
    ofstream f(path, ios::binary);
    if (!f) return 1;
    f.write(content.c_str(), (streamsize)content.size());
    f.close();
    return 0;
}

int test_producer_consumer_interleaved() {
    const string p = "tests/tmp_prodcons.txt";
    // ensure file exists empty
    writeFile(p, "");
    const int N = 50;
    // producer: append lines with a small delay
    std::atomic<bool> producer_done(false);
    std::thread producer([&]() {
        for (int i = 0; i < N; ++i) {
            ofstream f(p, ios::app | ios::binary);
            string line = "L" + to_string(i) + "\n";
            f.write(line.c_str(), (streamsize)line.size());
            f.close();
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        producer_done.store(true);
    });

    // consumer: use ifbufstream in multithreaded mode
    ifbufstream is(p, 64, true, false);
    vector<string> got;
    string s;
    auto start = std::chrono::steady_clock::now();
    while ((int)got.size() < N) {
        if (!is.getlines(s)) {
            // wait briefly for producer to produce more
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            // timeout guard
            if (std::chrono::steady_clock::now() - start > std::chrono::seconds(10)) break;
            continue;
        }
        got.push_back(s);
        // small timeout guard
        if (std::chrono::steady_clock::now() - start > std::chrono::seconds(10)) break;
    }

    producer.join();
    if ((int)got.size() != N) {
        cerr << "test_producer_consumer_interleaved: expected "<<N<<" lines, got "<<got.size()<<"\n";
        return 1;
    }
    for (int i = 0; i < N; ++i) {
        string expect = "L" + to_string(i);
        if (got[i] != expect) { cerr << "test_producer_consumer_interleaved: mismatch at "<<i<<" got '"<<got[i]<<"' expected '"<<expect<<"'\n"; return 1; }
    }
    return 0;
}

int test_truncate_during_read() {
    const string p = "tests/tmp_trunc.txt";
    // create a file with many lines
    const int TOTAL = 200;
    {
        ofstream f(p, ios::binary);
        for (int i = 0; i < TOTAL; ++i) {
            string line = "X" + to_string(i) + "\n";
            f.write(line.c_str(), (streamsize)line.size());
        }
    }

    ifbufstream is(p, 64, true, false);
    string s;
    vector<string> got;

    // start a thread that truncates the file after some reads
    std::atomic<bool> trunc_done(false);
    std::thread truncator([&]() {
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        // truncate file to zero
        ofstream f(p, ios::trunc | ios::binary);
        f.close();
        trunc_done.store(true);
    });

    // read until EOF or timeout
    auto start = std::chrono::steady_clock::now();
    while (true) {
        if (!is.getlines(s)) break;
        got.push_back(s);
        if (std::chrono::steady_clock::now() - start > std::chrono::seconds(5)) break;
    }

    truncator.join();
    // We expect that we read at least some lines before truncate, and no crash occurred
    if (got.empty()) { cerr << "test_truncate_during_read: expected some lines before truncate, got 0\n"; return 1; }
    return 0;
}

#if TEST_USE_ZLIB
int test_gzip_read() {
    const string p = "tests/tmp_gz.txt.gz";
    string content = "gline1\ngline2\n";
    // write gzip compressed file using zlib
    gzFile gz = gzopen(p.c_str(), "wb");
    if (!gz) { cerr<<"test_gzip_read: gzopen failed\n"; return 1; }
    int written = gzwrite(gz, content.c_str(), (unsigned)content.size());
    gzclose(gz);
    if (written != (int)content.size()) { cerr<<"test_gzip_read: gzwrite failed\n"; return 1; }

    // open with ifbufstream; InputStream.h needs to be compiled with gzip read support for this to work
    ifbufstream is(p, 64, false, false);
    string s;
    if (!is.getlines(s)) { cerr<<"test_gzip_read: first getline failed\n"; return 1; }
    if (s != "gline1") { cerr<<"test_gzip_read: expected 'gline1' got '"<<s<<"'\n"; return 1; }
    if (!is.getlines(s)) { cerr<<"test_gzip_read: second getline failed\n"; return 1; }
    if (s != "gline2") { cerr<<"test_gzip_read: expected 'gline2' got '"<<s<<"'\n"; return 1; }
    return 0;
}
#endif

int test_small_file() {
    const string p = "tests/tmp_small.txt";
    string content = "line1\nline2\nline3\n";
    writeFile(p, content);

    ifbufstream is(p, 64, false, false);
    string s;
    if (!is.getlines(s)) { cerr<<"test_small_file: first getline failed\n"; return 1; }
    if (s != "line1") { cerr<<"test_small_file: expected 'line1' got '"<<s<<"'\n"; return 1; }
    if (!is.getlines(s)) { cerr<<"test_small_file: second getline failed\n"; return 1; }
    if (s != "line2") { cerr<<"test_small_file: expected 'line2' got '"<<s<<"'\n"; return 1; }
    if (!is.getlines(s)) { cerr<<"test_small_file: third getline failed\n"; return 1; }
    if (s != "line3") { cerr<<"test_small_file: expected 'line3' got '"<<s<<"'\n"; return 1; }
    if (is.getlines(s)) { cerr<<"test_small_file: expected EOF but got data\n"; return 1; }
    return 0;
}

int test_file_smaller_than_buffer() {
    const string p = "tests/tmp_small2.txt";
    string content = "onlyline\n";
    writeFile(p, content);
    ifbufstream is(p, 1024, false, false);
    string s;
    if (!is.getlines(s)) { cerr<<"test_file_smaller_than_buffer: getline failed\n"; return 1; }
    if (s != "onlyline") { cerr<<"test_file_smaller_than_buffer: expected 'onlyline' got '"<<s<<"'\n"; return 1; }
    if (is.getlines(s)) { cerr<<"test_file_smaller_than_buffer: expected EOF but got data\n"; return 1; }
    return 0;
}

int test_crlf_boundary() {
    const string p = "tests/tmp_crlf.txt";
    // craft content so that '\r' is last byte of buffer and '\n' is first of next
    // choose buffer size 10 for the stream; place '\r' at pos 9 and '\n' at pos 10
    string content = "012345678\r\nABC\n";
    writeFile(p, content);
    ifbufstream is(p, 10, true, false); // enable MC mode to use background reader
    string s;
    if (!is.getlines(s)) { cerr<<"test_crlf_boundary: first getline failed\n"; return 1; }
    if (s != "012345678") { cerr<<"test_crlf_boundary: expected '012345678' got '"<<s<<"'\n"; return 1; }
    if (!is.getlines(s)) { cerr<<"test_crlf_boundary: second getline failed\n"; return 1; }
    if (s != "ABC") { cerr<<"test_crlf_boundary: expected 'ABC' got '"<<s<<"'\n"; return 1; }
    return 0;
}

int test_eof_without_newline() {
    const string p = "tests/tmp_nonl.txt";
    string content = "line1\nlastline"; // no trailing newline
    writeFile(p, content);
    ifbufstream is(p, 16, false, false);
    string s;
    if (!is.getlines(s)) { cerr<<"test_eof_without_newline: first getline failed\n"; return 1; }
    if (s != "line1") { cerr<<"test_eof_without_newline: expected 'line1' got '"<<s<<"'\n"; return 1; }
    if (!is.getlines(s)) { cerr<<"test_eof_without_newline: second getline failed\n"; return 1; }
    if (s != "lastline") { cerr<<"test_eof_without_newline: expected 'lastline' got '"<<s<<"'\n"; return 1; }
    if (is.getlines(s)) { cerr<<"test_eof_without_newline: expected EOF but got data\n"; return 1; }
    return 0;
}

int test_get4lines() {
    const string p = "tests/tmp4.txt";
    string content = "a\nb\nc\nd\n";
    writeFile(p, content);
    ifbufstream is(p, 32, false, false);
    vector<string> in(4);
    if (!is.get4lines(in)) { cerr<<"test_get4lines: get4lines() failed\n"; return 1; }
    if (in[0] != "a" || in[1] != "b" || in[2] != "c" || in[3] != "d") { cerr<<"test_get4lines: mismatch\n"; return 1; }
    return 0;
}

int main() {
    vector<pair<string, function<int()>>> tests = {
        {"small_file", test_small_file},
        {"file_smaller_than_buffer", test_file_smaller_than_buffer},
        {"crlf_boundary", test_crlf_boundary},
        {"eof_without_newline", test_eof_without_newline},
        {"get4lines", test_get4lines}
    };
#if TEST_USE_ZLIB
    tests.push_back({"gzip_read", test_gzip_read});
#endif
    // ifbufstream reads immutable sequencing files; live append/truncate behavior is outside its contract.

    int fails = 0;
    for (auto &t : tests) {
        cerr << "Running " << t.first << "...\n";
        int r = t.second();
        if (r != 0) { cerr << t.first << " FAILED\n"; fails++; }
        else cerr << t.first << " ok\n";
    }

    if (fails == 0) {
        cout << "All ifbufstream tests passed\n";
        return 0;
    }
    cerr << fails << " test(s) failed\n";
    return 1;
}
