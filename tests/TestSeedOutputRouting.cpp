#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "Filters.h"
#include "InputStream.h"
#include "OutputStreamer.h"
#include "containers.h"

using namespace std;

namespace {

void writeFile(const string& path, const string& contents) {
	ofstream out(path, ios::binary | ios::trunc);
	if (!out) { throw runtime_error("cannot create test file: " + path); }
	out << contents;
}

string readFile(const string& path) {
	ifstream in(path, ios::binary);
	if (!in) { return ""; }
	ostringstream contents;
	contents << in.rdbuf();
	return contents.str();
}

size_t countRecords(const string& fastq) {
	size_t count = 0;
	for (size_t pos = 0; (pos = fastq.find("\n+\n", pos)) != string::npos; pos += 3) {
		++count;
	}
	return count;
}

OptContainer seedOptions(const string& stem, const string& paf) {
	OptContainer options;
	options["-sample_sep"] = "___";
	options["-paired"] = "2";
	options["-onlyPair"] = "";
	options["-i_qual_offset"] = "33";
	options["-o_qual_offset"] = "33";
	options["-5PR1cut"] = "0";
	options["-5PR2cut"] = "0";
	options["-logLvsQ"] = "";
	options["-GoldenAxe"] = "0";
	options["-o_fastq"] = stem + ".1.fq," + stem + ".2.fq";
	options["-o_fna"] = "";
	options["-o_qual"] = "";
	options["-o_fastq2"] = "";
	options["-o_fastq_noBC"] = "";
	options["-suppressOutput"] = "0";
	options["-threadIO"] = "0";
	options["-merge_pairs_seed"] = "1";
	options["-uparseVer"] = "";
	options["-count_chimeras"] = "F";
	options["-optimalRead2Cluster"] = paf;
	return options;
}

void removeOutputs(const string& stem) {
	std::remove((stem + ".1.fq").c_str());
	std::remove((stem + ".2.fq").c_str());
	std::remove((stem + ".1.singl.fq").c_str());
	std::remove((stem + ".2.singl.fq").c_str());
	std::remove((stem + ".merg.fq").c_str());
}

void require(bool condition, const string& message) {
	if (!condition) { throw runtime_error(message); }
}

void testFailedMergeUsesOnlySingletonOutput() {
	const string stem = "tests/seed_route_failed";
	const string fallback = stem + ".fna";
	const string paf = stem + ".paf";
	removeOutputs(stem);
	writeFile(fallback, ">OTU_failed\nACGTACGT\n");
	writeFile(paf, "");

	OptContainer options = seedOptions(stem, paf);
	Filters filters(&options);
	{
		UClinks links(&options);
		auto fallbackInput = make_shared<InputStreamer>(true, 33, "0", "1", 1);
		fallbackInput->setupFna(fallback);
		links.setupDefSeeds(fallbackInput, {});
		OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
		output.activateReadMerger(1);
		links.writeNewSeeds(&output, &filters, false);
	}

	require(readFile(stem + ".merg.fq").empty(),
		"a failed merge was incorrectly written to merged output");
	require(countRecords(readFile(stem + ".1.singl.fq")) == 1,
		"a failed merge was not written exactly once to singleton output");

	std::remove(fallback.c_str());
	std::remove(paf.c_str());
	removeOutputs(stem);
}

} // namespace

int main() {
	ini_DNAconstants();
	OptContainer globalOptions;
	setGlobalOptContainer(&globalOptions);
	try {
		testFailedMergeUsesOnlySingletonOutput();
	}
	catch (const exception& error) {
		cerr << "Seed-output routing regression failed: " << error.what() << '\n';
		return 1;
	}
	cout << "Seed-output routing regression passed\n";
	return 0;
}
