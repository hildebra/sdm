
#include "GZipStr.h"

#ifdef _isa1gzip

GzipIfstream::GzipIfstream(const std::string& filename, size_t bufferSize = 8192)
	: std::istream(this) , bufferSize(bufferSize), decompressedPos(0), decompressedSize(0) {
	cerr << "Reading GzipIfstream::" << filename << endl;
	
	file.open(filename, ios::binary);
	if (!file.is_open()) {
		throw std::runtime_error("Error opening file: " + filename);
	}


	// Initialize ISA-L decompression structures
	isal_inflate_init(&state);
	isal_gzip_header_init(&gzipHeader);
	state.crc_flag = ISAL_GZIP_NO_HDR_VER;

	// Allocate buffers
	compressedData.resize(bufferSize);
	decompressedData.resize(bufferSize * 2); // Adjust size if needed

	state.next_out = decompressedData.data();
	state.avail_out = decompressedData.size();
	cerr<<"GzipIfstream created\n";
}


GzipIfstream& GzipIfstream::read(char* s, std::streamsize n) {
	
	cerr << "GZread" << n;
	std::streamsize bytesRead = 0;

	while (bytesRead < n) {
		if (decompressedPos >= decompressedSize) {
			if (!decompressNextChunk()) {
				break;
			}
		}

		std::streamsize toCopy = std::min(streamsize(n - bytesRead), streamsize(decompressedSize - decompressedPos));
		std::memcpy(s + bytesRead, decompressedData.data() + decompressedPos, toCopy);
		decompressedPos += toCopy;
		bytesRead += toCopy;
	}
	cerr<<"GZRdone";

	return *this;
}
/*int GzipIfstream::underflow() {
	if (decompressedPos >= decompressedSize) {
		if (!decompressNextChunk()) {
			return EOF;  // No more data
		}
	}

	setg(decompressedData.data(),
		decompressedData.data() + decompressedPos,
		decompressedData.data() + decompressedSize);
	return int(*gptr());
}
*/

GzipIfstream::~GzipIfstream() {
	if (file.is_open()) { file.close(); }
}
bool GzipIfstream::decompressNextChunk() {
	if (file.eof()) {
		return false;
	}
	cerr << "Decompr";
	file.read(reinterpret_cast<char*>(compressedData.data()), bufferSize);
	size_t bytesRead = file.gcount();
	if (bytesRead == 0) {
		return false;
	}

	state.next_in = compressedData.data();
	state.avail_in = bytesRead;
	state.next_out = decompressedData.data();
	state.avail_out = decompressedData.size();

	int ret = isal_inflate(&state);
	if (ret != ISAL_DECOMP_OK && ret != ISAL_END_INPUT) {
		throw std::runtime_error("Decompression failed with error code: " + std::to_string(ret));
	}

	decompressedSize = state.next_out - decompressedData.data();
	decompressedPos = 0;

	return decompressedSize > 0;
}











//GzipOfstream----------------
GzipOfstream::GzipOfstream(const std::string& filename, size_t bufferSize = 20 * 1024 * 1024)
	: std::ostream(this), bufferSize(bufferSize) {

	file.open(filename, std::ios::binary);
	if (!file.is_open()) {
		throw std::runtime_error("Error opening file: " + filename);
	}

	// Initialize ISA-L compression structures
	isal_deflate_init(&state);
	isal_write_gzip_header(&state, &gzipHeader);

	// Allocate buffers
	uncompressedData.resize(bufferSize);
	compressedData.resize(bufferSize * 2);  // Adjust as needed

	// Set up the stream buffer pointers
	setp(reinterpret_cast<char*>(uncompressedData.data()),
		reinterpret_cast<char*>(uncompressedData.data() + bufferSize));
}


GzipOfstream::~GzipOfstream() {
	sync();  // Ensure all remaining data is flushed
	if (file.is_open()) {
		file.close();
	}
}

bool GzipOfstream::compressAndWrite() {
	size_t dataSize = pptr() - pbase();
	state.next_in = uncompressedData.data();
	state.avail_in = dataSize;

	// Compress until all data is processed
	while (state.avail_in > 0) {
		state.next_out = compressedData.data();
		state.avail_out = compressedData.size();
		int ret = isal_deflate(&state);
		if (ret != COMP_OK ) {
			throw std::runtime_error("Compression failed with error code: " + std::to_string(ret));
		}

		// Write compressed data to the file
		size_t compressedSize = compressedData.size() - state.avail_out;
		file.write(reinterpret_cast<const char*>(compressedData.data()), compressedSize);
	}

	// Reset the put pointer to the beginning of the buffer
	setp(reinterpret_cast<char*>(uncompressedData.data()),
		reinterpret_cast<char*>(uncompressedData.data() + bufferSize));

	return file.good();
}





#endif