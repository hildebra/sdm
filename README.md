# sdm
Simple demultiplexing
This project started as a simple sequence demultiplexer, because I needed a fast and memory efficient program for this. Since then (5 years ago) it has grown in scope and complexity and is by now part of most of my pipelines, as the initial filtering step, that is able to handle a variety of standard sequence inputs (fasta / fastq all versions / gz of these), detecting inputs and even handling corrupt files (till the point where they are corrupted).
Thus this tool is a basic bioinformatic program, written in C++11 for fast and efficient functionality. It's main purposes are:
- parsing between DNA sequence formats (fasta / fastq)
- demultiplexing multi sample sequence files (based on DNA sequence or header keywords)
- removing sequencing / PCR primer from reads
- filtering DNA sequences based on a multitude of quality parameters
- fast extraction of a subset of DNA sequences, based on header names (also works for paired end reads)
- combining several inputs into one (filtered) sequence file
