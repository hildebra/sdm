
#pragma once

#include <sstream>

enum FileFormat {
    FORMAT_AUTO_DETECT,
    FORMAT_FASTA,
    FORMAT_FASTQ,
    FORMAT_FASTA_QUAL
};

/*enum FastqVersion {
    VERSION_AUTO_DETECT,
    SANGER,
    SOLEXA,
    ILLUMINA_1_3,
    ILLUMINA_1_5,
};*/

struct FastxRecord {
    FileFormat format;
    std::string header = "";  // header line, including @/>, but not newline
    std::string id = "";      // from first char. after @/> up to first whitespace
    std::string sequence = "";
    std::string quality = "";   // only meaningful for FASTQ seqs
    
    std::string &to_string();

private:
    std::string str_representation;
};

class BufferedFastxReader {
private:
    std::stringstream str_stream_;
    std::string str_buffer_ = "";
    FileFormat file_format_;
    char* block_buffer_;
    size_t block_buffer_size_;
    bool strip_space_ = true;
    
public:
    BufferedFastxReader();
    
    ~BufferedFastxReader();
    
    BufferedFastxReader(const BufferedFastxReader &rhs) = delete;
    
    BufferedFastxReader &operator=(const BufferedFastxReader &rhs) = delete;
    
    bool LoadBatch(std::istream &ifs, size_t record_count);
    
    bool LoadBlock(std::istream &ifs, size_t block_size);
    
    bool NextSequence(FastxRecord &record);
    
    static bool ReadNextSequence(std::istream &is, FastxRecord &record,
                                 std::string &str_buffer_ptr, FileFormat format = FORMAT_AUTO_DETECT);
    
    void auto_detect() { file_format_ = FORMAT_AUTO_DETECT; }
    
    FileFormat file_format() { return file_format_; }
};
