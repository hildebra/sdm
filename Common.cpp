#include "Common.h"
#include <algorithm>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <cctype>
#ifdef _gzipread
#include "include/gzstream.h"
#include "include/zstr.h"
#endif
// DNA_trans is defined in InputStream.cpp; declare extern here so utilities can use it
extern char DNA_trans[256];

void rtrim(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim leading and trailing whitespace
void trim(std::string& str){
    // trim trailing spaces
    size_t endpos = str.find_last_not_of(" \t");
    if( string::npos != endpos )
    {
        str = str.substr( 0, endpos+1 );
    }

    // trim leading spaces
    size_t startpos = str.find_first_not_of(" \t");
    if( string::npos != startpos )
    {
        str = str.substr( startpos );
    }
}

//from http://stackoverflow.com/questions/8888748/how-to-check-if-given-c-string-or-char-contains-only-digits
bool is_digits(const std::string &str)
{
    return std::all_of(str.begin(), str.end(), ::isdigit); // C++11
}



string spaceX(unsigned int k){
    string ret = "";
    for (unsigned int i = 0; i < k; i++){
        ret += ' ';
    }
    return ret;
}

int digitsInt(int x){
    int length = 1;
    while (x /= 10)
        length++;
    return length;
}

int digitsFlt(float x){
    std::stringstream s;
    s << x;
    return (int)s.str().length();
}

string intwithcommas(int value) {
    string numWithCommas = std::to_string((long long)value);
    int insertPosition = (int)numWithCommas.length() - 3;
    while (insertPosition > 0) {
        numWithCommas.insert(insertPosition, ",");
        insertPosition -= 3;
    }
    return (numWithCommas);
}

// itos removed - use std::to_string everywhere

std::string ftos(float number,int digits) {
    std::stringstream ss;
    ss << std::setprecision(digits)<< number;
    return ss.str();
}

static string normalizePathToken(string fi) {
    trim(fi);
    if (fi.size() >= 2) {
        const char first = fi.front();
        const char last = fi.back();
        if ((first == '"' && last == '"') || (first == '\'' && last == '\'')) {
            fi = fi.substr(1, fi.size() - 2);
            trim(fi);
        }
    }
    return fi;
}

bool isGZfile(const string fi) {
    string f = normalizePathToken(fi);
    if (f.size() < 3) {
        return false;
    }
    const size_t p = f.size() - 3;
    const char c0 = static_cast<char>(std::tolower(static_cast<unsigned char>(f[p])));
    const char c1 = static_cast<char>(std::tolower(static_cast<unsigned char>(f[p + 1])));
    const char c2 = static_cast<char>(std::tolower(static_cast<unsigned char>(f[p + 2])));
    return (c0 == '.' && c1 == 'g' && c2 == 'z');
}

static string normalizeInputPath(string fileS) {
    return normalizePathToken(std::move(fileS));
}

std::istream& safeGetline(std::istream& is, std::string& t) {
    t.clear();
    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if (sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            if (t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

int parseInt(const char** p1) {
    const char* p = *p1;
    while (*p == ' ') p++;

    int acc = 0;
    while (*p >= '0' && *p <= '9')
        acc = acc * 10 + *p++ - '0';

    *p1 = p;
    return acc;
}

std::vector<std::string> header_string_split(const std::string& str, const std::string& sep) {
    std::vector<std::string> tokens;
    tokens.reserve(13);
    size_t start = 0;
    size_t pos = 0;
    while ((pos = str.find_first_of(sep, start)) != std::string::npos) {
        tokens.push_back(str.substr(start, pos - start));
        start = pos + 1;
    }
    if (start < str.length()) {
        tokens.push_back(str.substr(start));
    } else if (start == str.length()) {
        tokens.push_back("0");
    }
    return tokens;
}

void remove_paired_info(string &s, short RP) {
    if (s.empty())  { return; }
    size_t f1 = s.find(" ");
    if (f1 != string::npos) {
        s.erase(f1);
        f1 = string::npos;
    }
    if (s.size() < 2) { return; }
    size_t  tarPos = s.size() - 2;
    f1 = s.rfind("/1");
    if (f1 == string::npos ||  f1 != tarPos) { f1 = s.rfind("/2"); }
    if (f1 == string::npos || f1 != tarPos) { f1 = s.rfind(".1"); }
    if (f1 == string::npos || f1 != tarPos) { f1 = s.rfind(".2"); }
    if (f1 != string::npos && f1 == tarPos) {
        s.erase(f1);
    }
}

std::string header_stem(const std::string& header) {
    if (header.empty()) { return ""; }
    const size_t slash = header.find('/');
    if (slash != std::string::npos) {
        return header.substr(0, slash);
    }

    std::array<size_t, 12> delim_pos;
    size_t start = 0;
    size_t idx = 0;
    while (idx < 11) {
        size_t p = header.find(": ", start);
        if (p == std::string::npos) break;
        delim_pos[idx++] = p;
        start = p + 2;
    }

    auto token_at = [&](size_t i) -> std::string_view {
        if (i == 0) {
            size_t end = (idx > 0) ? delim_pos[0] : header.size();
            return std::string_view(header.data(), end);
        }
        if (i < idx) {
            size_t s = delim_pos[i - 1] + 2;
            size_t e = delim_pos[i];
            return std::string_view(header.data() + s, e - s);
        }
        if (idx > 0 && i == idx) {
            size_t s = delim_pos[idx - 1] + 2;
            size_t e = header.size();
            return std::string_view(header.data() + s, e - s);
        }
        return std::string_view();
    };

    size_t tokens_n = (idx == 0) ? 1 : (idx + 1);

    if (tokens_n == 11) {
        std::string out;
        std::string_view t0 = token_at(0);
        std::string_view t3 = token_at(3);
        std::string_view t4 = token_at(4);
        std::string_view t5 = token_at(5);
        std::string_view t6 = token_at(6);
        out.reserve(t0.size() + t3.size() + t4.size() + t5.size() + t6.size() + 4);
        out.append(t0.data(), t0.size()); out.push_back(':');
        out.append(t3.data(), t3.size()); out.push_back(':');
        out.append(t4.data(), t4.size()); out.push_back(':');
        out.append(t5.data(), t5.size()); out.push_back(':');
        out.append(t6.data(), t6.size());
        return out;
    } else if (tokens_n == 6 || tokens_n == 7) {
        std::string out;
        std::string_view t0 = token_at(0);
        std::string_view t2 = token_at(2);
        std::string_view t3 = token_at(3);
        std::string_view t4 = token_at(4);
        std::string_view t5 = token_at(5);
        out.reserve(t0.size() + t2.size() + t3.size() + t4.size() + t5.size() + 4);
        out.append(t0.data(), t0.size()); out.push_back(':');
        out.append(t2.data(), t2.size()); out.push_back(':');
        out.append(t3.data(), t3.size()); out.push_back(':');
        out.append(t4.data(), t4.size()); out.push_back(':');
        out.append(t5.data(), t5.size());
        return out;
    } else {
        throw std::runtime_error(std::string("fastq header format error: ") + header);
    }
}

void reverseTS(string & Seq) {
    int qs = (int)Seq.length() - 1;
    string S2 = Seq.c_str();
    for (int i = qs; i >= 0; i--) {
        Seq[i] = DNA_trans[(unsigned char)S2[qs - i]];
    }
}
string reverseTS2(const string & Seq) {
    int qs = (int)Seq.length() - 1;
    string S2 = Seq.c_str();
    for (int i = qs; i >= 0; i--) {
        S2[i] = DNA_trans[(unsigned char)Seq[qs - i]];
    }
    return S2;
}

bool any_lowered(const string& is) {
    for (unsigned int i = 0; i < is.length(); i++) {
        char c = is[i];
        if (std::islower(static_cast<unsigned char>(c))) { return true; }
    }
    return false;
}

std::ptrdiff_t len_common_prefix_base(char const a[], char const b[])
{
    if (std::strlen(b) > std::strlen(a)) {
        return std::distance(a, std::mismatch(a, a + std::strlen(a), b).first);
    }
    return std::distance(b, std::mismatch(b, b + std::strlen(b), a).first);
}

string applyFileIT(string x, int it,const string xtr){

    size_t pos = x.find_last_of(".");
    if (pos != string::npos && isGZfile(x)) {
        pos = x.find_last_of(".", pos-1);
    }
    if (it == 0) {
        if (pos == string::npos) {
            return x + xtr;
        }
        return x.substr(0, pos) + xtr + x.substr(pos);
    }
    std::ostringstream ss;
    ss << it;
    if (pos == string::npos) {
        return x + "." + ss.str() + xtr ;
    }
    return x.substr(0,pos) + "."+ss.str() + xtr + x.substr(pos);
}

bool fileExists(const std::string& name, int i, bool extiffail) {
    if (name == "") { return false; }
    std::unique_ptr<FILE, decltype(&fclose)> file(fopen(name.c_str(), "r"), &fclose);
    if (file) {
        return true;
    } else {
        if (extiffail) {
            std::cerr << "ERROR: Could not find file " << name << std::endl;
            if (i >= 0) {
                std::cerr << "on mapping file line " << i << std::endl;
            }
            exit(92);
        }
        return false;
    }
}

string detectSeqFmt(const string inF) {

    vector<string> tfasP = splitByCommas(inF, ';');
    for (size_t i = 0; i < tfasP.size(); i++) {
        vector<string> tfas = splitByCommas(tfasP[i]);
        string fileS = normalizeInputPath(tfas[0]);
        std::unique_ptr<std::istream> fnax;
        string file_type = "test file";
        string tmp("");
        string ret = "";
        if (fileS != "") {
            string fLower = fileS;
            std::transform(fLower.begin(), fLower.end(), fLower.begin(), [](unsigned char c) { return (char)std::tolower(c); });

            // Prefer extension-based detection for gzip files to avoid stream backend issues during auto-detect.
            if (isGZfile(fileS)) {
                string base = fLower.substr(0, fLower.size() - 3);
                if (base.size() >= 6 && (base.compare(base.size() - 6, 6, ".fastq") == 0)) {
                    return "-i_fastq";
                }
                if (base.size() >= 3 && (base.compare(base.size() - 3, 3, ".fq") == 0)) {
                    return "-i_fastq";
                }
                if (base.size() >= 6 && (base.compare(base.size() - 6, 6, ".fasta") == 0)) {
                    return "-i_fna";
                }
                if (base.size() >= 4 && (base.compare(base.size() - 4, 4, ".fna") == 0 || base.compare(base.size() - 3, 3, ".fa") == 0)) {
                    return "-i_fna";
                }
            }

#ifdef _gzipread
            if (isGZfile(fileS)) {
                fnax.reset(new zstr::ifstream(fileS.c_str()));
            }
            else {
                fnax.reset(new std::ifstream(fileS.c_str(), std::ios::in));
            }
#else
            fnax.reset(new std::ifstream(fileS.c_str(), std::ios::in));
#endif

            if (!*(fnax)) {
                std::cerr << "\nCouldn't open " << file_type << " file \"" << fileS << "\"! Skipping..\n";
                continue;
            }
            while (safeGetline(*fnax, tmp)) {
                if (tmp.empty()) { continue; }
                if (tmp[0] == '>') {
                    ret = "-i_fna"; break;
                }
                else if (tmp[0] == '@') {
                    ret = "-i_fastq"; break;
                }
                else if (tmp.length() == 0) {
                    ;
                }
                else {
					bool isGZ = isGZfile(fileS);
                    std::cerr << "\nCouldn't open " << file_type << " file \"" << fileS << "\"! Skipping..\n";
                    exit(888);
                }
            }
        }
        if (ret == "") {
            std::cerr << "Empty input file detected:\n" << fileS << std::endl;
        } else {
            return ret;
        }

    }
    return "empty";
}
