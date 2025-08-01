#!/bin/bash

#SBATCH --job-name="histo_distribution"
#SBATCH --output='histogram.out'
#SBATCH --account='bgmp'
#SBATCH --partition='bgmp'

cd /home/alho/bgmp/alho/bioinfo/Bi622/Demultiplexing/Assignment-the-first/

mamba activate base

/usr/bin/time -v python3 /home/alho/bgmp/alho/bioinfo/Bi622/Demultiplexing/Assignment-the-first/NT_distribution.py
