#pragma once
/* sdm: simple demultiplexer
Copyright (C) 2026  Falk Hildebrand

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string>
#include <vector>
#include <sstream>
#include <array>
#include <algorithm>
#include <cctype>
#include <istream>
#include <memory>
#include <iostream>

using std::string;
using std::vector;

// Utility functions used across DNA and InputStream
string spaceX(unsigned int k);
int digitsInt(int x);
int digitsFlt(float x);
string intwithcommas(int value);
// itos was removed but kept as a thin wrapper to std::to_string for backward compatibility.
// Replace calls to itos(...) with std::to_string(...) across the codebase when convenient.
inline std::string itos(int number) { return std::to_string(number); }
std::string ftos(float number, int digits = 4);
bool isGZfile(const string fileS);

// trim leading and trailing whitespace
void trim(std::string& str);

// check if string contains only digits
bool is_digits(const std::string &str);

// trim whitespace from end
void rtrim(std::string& s);

// fast integer parse from C-string pointer: advances pointer
int parseInt(const char** p1);

std::istream& safeGetline(std::istream& is, std::string& t);

std::vector<std::string> header_string_split(const std::string& str, const std::string& sep);
void remove_paired_info(string&, short = -1);
std::string header_stem(const std::string& header);

string reverseTS2(const std::string& Seq);
void reverseTS(std::string& Seq);

bool any_lowered(const string& is);
std::ptrdiff_t len_common_prefix_base(char const a[], char const b[]);

string applyFileIT(string x, int it, const string xtr = "");
bool fileExists(const std::string& name, int i = -1, bool extiffail = true);

inline vector<string> splitByComma(const string& fileS, bool requireTwo, char SrchStr = ',') {
    string::size_type pos = fileS.find(SrchStr);
    if (pos == string::npos) {
        if (requireTwo) {
            std::cerr << fileS << std::endl;
            std::cerr << "Could not find '" << SrchStr << "' in input file (required for paired end sequences)\n";
            exit(14);
        }
        else {
            return vector<string>(1, fileS);
        }
    }
    vector<string> tfas(2, "");
    tfas[0] = fileS.substr(0, pos);
    tfas[1] = fileS.substr(pos + 1);
    return tfas;
}

inline vector<string> splitByCommas(const string& fileS, char SrchStr = ',') {
    if (fileS.find(SrchStr) == string::npos) { return vector<string>(1, fileS); }
    vector<string> res = splitByComma(fileS, true, SrchStr);
    vector<string> ret; ret.push_back(res[0]);
    while (res[1].find(SrchStr) != string::npos) {
        res = splitByComma(res[1], true, SrchStr);
        ret.push_back(res[0]);
    }
    ret.push_back(res[1]);
    return ret;
}

template<class TYPE>
TYPE calc_median2(vector<TYPE>& in, float perc) {
    size_t sum = in.size();
    if (sum == 0) return TYPE();
    size_t tar = (size_t)(((float)(sum - 1)) * perc);
    return in[tar];
}

template<class TYPE>
TYPE median(vector<TYPE>& v) {
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin() + n, v.end());
    return v[n];
}

string detectSeqFmt(const string);

