#pragma once
#include "../DNAconsts.h"

#ifdef _isa1gzip


#include <streambuf>
#include <istream>
#include <iostream>
//#include <strstream>


#include <isa-l/igzip_lib.h>

class GzipIfstream : public std::istream, private std::streambuf {
public:
	GzipIfstream(const std::string& filename, size_t bufferSize );
	~GzipIfstream();
	bool good() const {
		return (file.good() || decompressedPos < decompressedSize);
	}

	GzipIfstream& read(char* s, std::streamsize n);

	bool eof() const {
		return file.eof() && decompressedPos >= decompressedSize;
	}

protected:
	// Override underflow to manage decompression and buffering
	std::streambuf::int_type underflow() override {
		if (decompressedPos >= decompressedSize) {
			if (!decompressNextChunk()) {
				return std::streambuf::traits_type::eof();  // No more data
			}
		}

		// Cast `decompressedData.data()` from `unsigned char*` to `char*`
		setg(reinterpret_cast<char*>(decompressedData.data()),
			reinterpret_cast<char*>(decompressedData.data()) + decompressedPos,
			reinterpret_cast<char*>(decompressedData.data()) + decompressedSize);

		return std::streambuf::traits_type::to_int_type(*gptr());
	}
	//int underflow()override;
private:
	bool decompressNextChunk();

	std::ifstream file;
	size_t bufferSize;
	std::vector<uint8_t> compressedData;
	std::vector<uint8_t> decompressedData;
	struct inflate_state state;
	struct isal_gzip_header gzipHeader;
	size_t decompressedPos;
	size_t decompressedSize;
};




class GzipOfstream : public std::ostream, private std::streambuf {
public:
    GzipOfstream(const std::string& filename, size_t bufferSize );
    ~GzipOfstream();

protected:
    // Called when the buffer is full, compress and write data
    std::streambuf::int_type overflow(std::streambuf::int_type ch) override {
        if (ch != std::streambuf::traits_type::eof()) {
            *pptr() = ch;  // Add character to buffer
            pbump(1);      // Move put pointer forward
        }

        // Compress the current buffer contents
        return compressAndWrite() ? ch : std::streambuf::traits_type::eof();
    }

    // Flushes any remaining data in the buffer
    int sync() override {
        return compressAndWrite() ? 0 : -1;
    }

private:
	bool compressAndWrite();

    std::ofstream file;
    size_t bufferSize;
    std::vector<uint8_t> uncompressedData;
    std::vector<uint8_t> compressedData;
    struct isal_zstream state;
    struct isal_gzip_header gzipHeader;
};






#endif
