

#include "FastxReader.h"
#include <cctype>
#include <fcntl.h>
#include <cstdio>
//#include <err.h>
//#include <sysexits.h>
#include <iostream>
#include <ctype.h>

using std::string;

inline void StripString(string &str) {
    while (isspace(str.back()))
        str.pop_back();
}

string& FastxRecord::to_string() {
    str_representation.assign(header);
    switch (format) {
        case FORMAT_FASTQ:
            str_representation.append("\n");
            str_representation.append(sequence);
            str_representation.append("\n+\n");
            str_representation.append(quality);
            str_representation.append("\n");
            break;
        default:
            str_representation.append("\n");
            str_representation.append(sequence);
            str_representation.append("\n");
            break;
    }
    return str_representation;
}

BufferedFastxReader::BufferedFastxReader() {
    file_format_ = FORMAT_AUTO_DETECT;
    str_buffer_.reserve(8192);
    block_buffer_ = new char[8192];
    block_buffer_size_ = 8192;
}


BufferedFastxReader::~BufferedFastxReader() {
    delete[] block_buffer_;
}

bool BufferedFastxReader::LoadBlock(std::istream &ifs, size_t block_size) {
    str_stream_.clear();
    str_stream_.str("");
    if (block_buffer_size_ < block_size) {
        delete[] block_buffer_;
        block_buffer_ = new char[block_size];
        block_buffer_size_ = block_size;
    }
    ifs.read(block_buffer_, block_size);
    if (! ifs && ifs.gcount() <= 0)
        return false;
    
    if (file_format_ == FORMAT_AUTO_DETECT) {
        switch (block_buffer_[0]) {
            case '@' : file_format_ = FORMAT_FASTQ; break;
            case '>' : file_format_ = FORMAT_FASTA; break;
            default:
				std::cerr << "sequence reader - unrecognized file format"; exit(55);
        }
    }
    str_buffer_.assign(block_buffer_, ifs.gcount());
    str_stream_ << str_buffer_;
    if (getline(ifs, str_buffer_))
        str_stream_ << str_buffer_ << "\n";
    if (file_format_ == FORMAT_FASTQ) {
        while (getline(ifs, str_buffer_)) {
            str_stream_ << str_buffer_ << "\n";
            if (str_buffer_[0] == '@')
                break;
        }
        int lines_to_read = 0;
        if (getline(ifs, str_buffer_)) {
            str_stream_ << str_buffer_ << "\n";
            lines_to_read = str_buffer_[0] == '@' ? 3 : 2;
            while (lines_to_read-- > 0 && getline(ifs, str_buffer_))
                str_stream_ << str_buffer_ << "\n";
        }
    }
    else {
        while (ifs) {
            if (ifs.peek() == '>')
                break;
            if (getline(ifs, str_buffer_))
                str_stream_ << str_buffer_ << "\n";
        }
    }
    return true;
}

bool BufferedFastxReader::LoadBatch(std::istream &ifs, size_t record_count) {
    str_stream_.clear();
    str_stream_.str("");
    auto valid = false;
    if (file_format_ == FORMAT_AUTO_DETECT) {
        if (!ifs)
            return false;
        switch (ifs.peek()) {
            case '@' :
                file_format_ = FORMAT_FASTQ;
                break;
            case '>' :
                file_format_ = FORMAT_FASTA;
                break;
            case EOF :
                return false;
            default:
				std::cerr << "sequence reader - unrecognized file format"; exit(55);
		}
        valid = true;
    }
    
    size_t line_count = 0;
    while (record_count > 0 && ifs) {
        if (getline(ifs, str_buffer_))
            line_count++;
        valid = true;
        if (file_format_ == FORMAT_FASTQ) {
            if (line_count % 4 == 0)
                record_count--;
        } else {
            if (ifs.peek() == '>')
                record_count--;
        }
        str_stream_ << str_buffer_ << "\n";
    }
    
    return valid;
}

bool BufferedFastxReader::NextSequence(FastxRecord &seq) {
    return BufferedFastxReader::ReadNextSequence
            (str_stream_, seq, str_buffer_, file_format_);
}

bool BufferedFastxReader::ReadNextSequence(std::istream &is, FastxRecord &record,
                                           std::string &str_buffer, FileFormat file_format) {
    if (!getline(is, str_buffer))
        return false;
    StripString(str_buffer);
    if (file_format == FORMAT_AUTO_DETECT) {
        switch (str_buffer[0]) {
            case '@' :
                file_format = FORMAT_FASTQ;
                break;
            case '>' :
                if (!isdigit(is.peek())) {
                    file_format = FORMAT_FASTA;
                } else {
                    file_format = FORMAT_FASTA_QUAL;
                }
                break;
            default:
				std::cerr << "sequence reader - unrecognized file format"; exit(55);
		}
    }
    record.format = file_format;
    if (record.format == FORMAT_FASTQ) {
        if (str_buffer.empty()) // Allow empty line to end file
            return false;
		if (str_buffer[0] != '@') {
			std::cerr << "malformed FASTQ file (exp. '@', saw \""<< str_buffer<< "), aborting"; exit(55);

		}
    } else if (record.format == FORMAT_FASTA || record.format == FORMAT_FASTA_QUAL) {
		if (str_buffer[0] != '>') {
			std::cerr << "malformed FASTQ file (exp. '>', saw \"" << str_buffer << "), aborting"; exit(55);

		}
  }
	else {
		std::cerr << "illegal sequence format encountered in parsing"; exit(55);

	}
    record.header.assign(str_buffer);
    auto first_whitespace_ch = str_buffer.find_first_of(" \t\r", 1);
    auto substr_len = first_whitespace_ch;
    if (substr_len != std::string::npos)
        substr_len--;
    if (str_buffer.size() > 1)
        record.id.assign(str_buffer, 1, substr_len);
    else
        return false;
    
    
    if (record.format == FORMAT_FASTQ) {
        if (!getline(is, str_buffer))
            return false;
        StripString(str_buffer);
        record.sequence.assign(str_buffer);
        if (!getline(is, str_buffer))  //  + line, discard
            return false;
        if (!getline(is, str_buffer))
            return false;
        StripString(str_buffer);
        record.quality.assign(str_buffer);
    } else if (record.format == FORMAT_FASTA) {
        record.quality.assign("");
        record.sequence.assign("");
        while (is && is.peek() != '>') {
            if (!getline(is, str_buffer))
                return !record.sequence.empty();
            StripString(str_buffer);
            record.sequence.append(str_buffer);
        }
    } else if (record.format == FORMAT_FASTA_QUAL) {
        record.quality.assign("");
        record.sequence.assign("");
        while (is && is.peek() != '>') {
            if (!getline(is, str_buffer))
                return !record.sequence.empty();
            StripString(str_buffer);
            record.sequence.append(str_buffer);
            record.sequence.append(" ");
        }
        StripString(record.sequence);
    }
    return true;
}

