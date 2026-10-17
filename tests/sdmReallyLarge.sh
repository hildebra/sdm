#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -o /hpc-home/hildebra/dev/Cpp/sdm/tests/sdmReallyLarge.sh.otxt
#SBATCH -e /hpc-home/hildebra/dev/Cpp/sdm/tests/sdmReallyLarge.sh.etxt
#SBATCH --mem=61440
#SBATCH --export=ALL
#SBATCH -p "nbi-medium,ei-medium,qib-medium,nbi-long,ei-long,qib-long"
#SBATCH --exclude=t512n12
#SBATCH -J 45_SDM701
echo $HOSTNAME;
set -e
ulimit -c 0;


#sdm2.13: 8h50m
OUT=/ei/projects/6/6209c9f1-a78a-4f8e-aa49-16176b0cc246/scratch//sdmTest1234/

mkdir -p $OUT/
rm -f $OUT/*
/hpc-home/hildebra/dev/Cpp/sdm/tests/sdm -i "/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/cloudpool/data/raw/Public/PRJNA763692//SRR17827816_1.fastq.gz,/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/cloudpool/data/raw/Public/PRJNA763692//SRR17827816_2.fastq.gz" -o_fastq $OUT/filtered.1.fq.gz,$OUT/filtered.2.fq.gz -options /hpc-home/hildebra/dev/Cpp/sdm/tests/sdmo_150_.txt -paired 2  -log $OUT/filter.log  -ignore_IO_errors 1 -i_qual_offset auto -binomialFilterBothPairs 1 -threads 4  -illuminaClip 1  -5PR1cut 0 -5PR2cut 0 

