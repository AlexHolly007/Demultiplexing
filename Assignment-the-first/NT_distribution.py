import gzip
import bioinfo
import numpy as np
import math as m



def populate_list(file: str) -> tuple:
    """Update with your own docstring"""
    my_list = np.zeros(101, dtype=float)
    record_count = 0
    with gzip.open(file, 'r') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                record_count += 1
                for j, letter in enumerate(line.strip()):
                    my_list[j] += bioinfo.convert_phred(str(letter))
    
    phred_means = my_list / record_count
    var = np.zeros(101, dtype=float)

    with gzip.open(file, 'r') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                for j, letter in enumerate(line.strip()):
                    var[j] += (bioinfo.convert_phred(str(letter)) - phred_means)**2

    var = var / record_count

    stdev = m.sqrt(var)
    
    return (phred_means, var, stdev, record_count)


def main():
    # for file in ["/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz",
    # "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz",
    # "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz",
    # "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"]:
    for file in './test.fastq.gz':
        phred_means, var, record_cnt = populate_list(file)

        