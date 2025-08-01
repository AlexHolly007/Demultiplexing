#!/usr/bin/env python3
import gzip
import bioinfo
import numpy as np
import math as m
import matplotlib.pyplot as plt

TESTING = False

def populate_list(file: str) -> tuple:
    """Update with your own docstring"""
    my_list = None
    record_count = 0
    with gzip.open(file, 'rb') as fh:
        first_record = [fh.readline() for _ in range(4)]
        seq_len = len(first_record[1])-1
        print("uh oh")
        if TESTING:
            assert(seq_len == 8), "seq_length not correct"
        my_list = np.zeros(seq_len, dtype=float)
        record_count += 1

        for i, letter in enumerate(first_record[3].strip()):
                    my_list[i] += (ord(chr(letter)) - 33)
        
        while True:
            head = fh.readline()
            if not head:
                break
            _ = fh.readline()
            _ = fh.readline()
            qual= fh.readline()
            record_count += 1
            
            if record_count == 1000:
                print("1000")
            if record_count % 10000000 == 0:
                    print(f"{file[file.find("_R")+1:file.find("_R")+3]} lines p1: {record_count}")
                    
            for j, letter in enumerate(qual.strip()):
                my_list[j] += (ord(chr(letter)) - 33)

    phred_means = my_list / record_count

    if TESTING:
        assert(phred_means[1] == 30.75), "Mean test passed"
        assert(record_count == 4)

    var = np.zeros(seq_len, dtype=float)

    with gzip.open(file, 'rb') as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                if (i/4) % 10000000 == 0:
                    print(f"{file[file.find("_R")+1:file.find("_R")+3]} lines p2: {record_count}")
                for j, letter in enumerate(line.strip()):
                    var[j] += ((ord(chr(letter)) - 33) - phred_means[j])**2

    var = var / record_count

    stdev = np.sqrt(var)

    if TESTING:
        assert(var[1] == 256.6875)
        assert(stdev[1] == 16.021469970012117)
    

    plt.errorbar(range(seq_len), phred_means, stdev, ecolor='r')
    plt.ylim(bottom=0)
    plt.title(f"{file[file.find("_R")+1:file.find("_R")+3]}Mean error plt")
    plt.xlabel("Score position")
    plt.ylabel("Qual score")
    plt.savefig(f"{file[file.find("_R")+1:file.find("_R")+3]}.png")


    return phred_means, var, stdev, record_count

def main():
    if TESTING:
        phred_means, var, stdev, record_cnt = populate_list('../TEST-input_FASTQ/Test_R2.fq.gz')


    else:
        for file in ["/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz",
        "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz",
        "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz",
        "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"]:
            phred_means, var, stdev, record_cnt = populate_list(file)

if __name__ == "__main__":
    main()
        
