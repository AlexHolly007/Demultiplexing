#!/bin/bash

#SBATCH --job-name="demultiplex"
#SBATCH --output='demultiu.out'
#SBATCH --account='bgmp'
#SBATCH --partition='bgmp'

cd /home/alho/bgmp/alho/bioinfo/Bi622/Demultiplexing/Assignment-the-third/

mamba activate base

# /usr/bin/time -v python3 /home/alho/bgmp/alho/bioinfo/Bi622/Demultiplexing/Assignment-the-third/A3_the_third.py -R1 ../TEST-input_FASTQ/Test_R1.fq.gz -R2 ../TEST-input_FASTQ/Test_R2.fq.gz -R3 ../TEST-input_FASTQ/Test_R3.fq.gz -R4 ../TEST-input_FASTQ/Test_R4.fq.gz -i ../TEST-input_FASTQ/indexes.txt -o results -os stat_output.txt
/usr/bin/time -v python3 /home/alho/bgmp/alho/bioinfo/Bi622/Demultiplexing/Assignment-the-third/A3_the_third.py -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -R3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -R4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i /projects/bgmp/shared/2017_sequencing/indexes.txt -o results -os stat_output.txt
