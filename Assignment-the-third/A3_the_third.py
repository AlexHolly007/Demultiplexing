#!/usr/bin/env python3
import bioinfo
import argparse
import gzip

C_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

#Reverse compliment of a DNA sequence that is sent in
def reverse_compliment(sequence: str) -> str:
    return str([C_DICT[x] for x in sequence])[::-1]

#Fills a file dictionary with the input index file. 
# Must be in format where 4th column is index name, and 5th is index sequence
def fill_file_dict(index_file):
    file_dict = {}

    return file_dict
    


def main(R1, R2, R3, R4 ,idx_file):
    #This function opens all files that are going to be written to using the index file input
    #creates a dictionary for all file handles. file stream for index A5 would be file_dict['A5_R2']. 
    # 'UK_R2' for unknown R2.   'UM_R1' for unmatched R1
    file_dict = fill_file_dict(idx_file)

    with gzip.open(R2, 'rt') as r2:
        head, seq, extra, qual = [r2.readline() for i in range(4)]



def get_args():
    parser = argparse.ArgumentParser(description='args for Demultiplexing')
    parser.add_argument("-R1", "--read1", required=True, type=str, help="the read1 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R2", "--read2", required=True, type=str, help="the read2 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R3", "--read3", required=True, type=str, help="the read3 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R4", "--read4", required=True, type=str, help="the read4 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument('-c','--cutoff', required=False, default=0, type=str, help="The minimum of what the mean quality score of an index needs to be for it to be not thrown to unknown")
    parser.add_argument('-i','--index_file', required=True, type=str, help="The index file is what gives names and barcodes to each index that can be found")

    return parser.parse_args()
    
#./A3_the_third.py -R1 ../TEST-input_FASTQ/Test_R1.fq.gz -R2 ../TEST-input_FASTQ/Test_R2.fq.gz -R3 ../TEST-input_FASTQ/Test_R3.fq.gz -R4 ../TEST-input_FASTQ/Test_R4.fq.gz -i ../TEST-input_FASTQ/indexes.txt

if __name__ == "__main__":
    args = get_args()
    main(args.read1, args.read2, args.read3, args.read4, args.index_file)