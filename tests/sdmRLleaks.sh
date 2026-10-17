#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -o /hpc-home/hildebra/dev/Cpp/sdm/tests/sdmRLleaks.sh.otxt
#SBATCH -e /hpc-home/hildebra/dev/Cpp/sdm/tests/sdmRLleaks.sh.etxt
#SBATCH --mem=61440
#SBATCH --export=ALL
#SBATCH -p "nbi-medium,ei-medium,qib-medium,nbi-long,ei-long,qib-long"
#SBATCH --exclude=t512n12
#SBATCH -J 45_SDM701
echo $HOSTNAME;
set -e
ulimit -c 0;


#sdm2.13: 8h50m
OUT=/ei/projects/6/6209c9f1-a78a-4f8e-aa49-16176b0cc246/scratch//sdmTest1234D/

mkdir -p $OUT/
rm -f $OUT//filterDone.stone  $OUT/filtered.singl.fq.gz

#--leak-check=yes --tool=drd --read-var-info=yes
valgrind --leak-check=yes /hpc-home/hildebra/dev/Cpp/sdm/tests/sdmD -i "/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/cloudpool/data/raw/Public/PRJNA763692//SRR17827816_1.sm.fastq,/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/cloudpool/data/raw/Public/PRJNA763692//SRR17827816_2.sm.fastq" -o_fastq $OUT/filtered.1.fq.gz,$OUT/filtered.2.fq.gz -options sdmo_150_.txt -paired 2  -log /hpc-home/hildebra/grp/data/projects/postDrama/Mong.C83//LOGandSUB//sdm/filter.log  -ignore_IO_errors 1 -i_qual_offset auto -binomialFilterBothPairs 1 -threads 4  -illuminaClip 1  -5PR1cut 0 -5PR2cut 0 
