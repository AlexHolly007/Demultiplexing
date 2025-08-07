#!/usr/bin/env python3
import bioinfo
import argparse
import gzip
import os

#maps each nucleotide to its compliment
C_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
#Maps each barcode to its matching group name
BARCODE_MAP = {}
#hold all the open files for each index
FILE_DICT = {}

#holds each found barcode pair. Both matching and index hopped
MATCH_DICT: dict[str, int] = {}


#Reverse compliment of a DNA sequence that is sent in
def reverse_compliment(sequence: str) -> str:
    return ''.join([C_DICT[x] for x in sequence][::-1])

#Fills a file dictionary with the input index file. 
# Must be in format where 4th column is index name, and 5th is index sequence
def fill_file_dict(index_file: str, output_dir: str):

    os.makedirs(output_dir, exist_ok=True)

    with open(index_file) as fi:
        for line in fi:

            #create two files for each index
            index = line.split()[3]
            for i in range(2):
                ft = open(f'{output_dir}/{index}_R{i+1}.fq', 'w')

                FILE_DICT[f'{index}_R{i+1}'] = [ft, 0]

            #Fill the barcode mapping for this index
            barcode = line.split()[4]
            BARCODE_MAP[barcode] = index
    
    for i in range(2):
        ft = open(f'{output_dir}/unknown_R{i+1}.fq', 'w')
        FILE_DICT[f'UK_R{i+1}'] = ft

        ft2 = open(f'{output_dir}/hopped_R{i+1}.fq', 'w')
        FILE_DICT[f'Hop_R{i+1}'] = ft2

    #Fill all possible pairs of indexes, including matching and hopped indexes
    for barcode_key in BARCODE_MAP:
        for barcode_key2 in BARCODE_MAP:
            MATCH_DICT[f'{BARCODE_MAP[barcode_key]}_{BARCODE_MAP[barcode_key2]}'] = 0





#This function will write the the R1 and R2 files for a pair of indexes.
# The file_key parameter, will point to which 2 files that is, so weither the indexes match/ dont match should be already known
def write_reads(file_key, head1, head4, seq1, seq4, qual1, qual4, barcode2, barcode3):
    #R1
    file = FILE_DICT[f'{file_key}_R1']
    file.write(f'{head1}_{barcode2}-{barcode3}\n')
    file.write(f'{seq1}\n')
    file.write(f'+\n{qual1}\n')

    #R2
    file = FILE_DICT[f'{file_key}_R2']
    file.write(f'{head4}_{barcode2}-{barcode3}\n')
    file.write(f'{seq4}\n')
    file.write(f'+\n{qual4}\n')
            
    


def main(R1, R2, R3, R4, cutoff, idx_file, output_dir, output_file):
    #This function opens all files that are going to be written to using the index file input
    #creates a dictionary for all file handles. file stream for index A5 R2 would be file_dict['A5_R2']. 
    # 'UK_R2' for unknown R2.   'Hop_R1' for hopped reads

    #Fills the global variable
    fill_file_dict(idx_file, output_dir)

    with gzip.open(R2, 'rt') as r2, gzip.open(R3, 'rt') as r3, gzip.open(R1, 'rt') as r1, gzip.open(R4, 'rt') as r4:

        properly_matched, index_hopped, unknown_indexes = 0,0,0

        while True:
            #grab records
            head1, seq1, _, qual1 = [r1.readline().strip() for _ in range(4)]
            _, barcode2, _, qual2 = [r2.readline().strip() for _ in range(4)]
            _, barcode3, _, qual3 = [r3.readline().strip() for _ in range(4)]
            head4, seq4, _, qual4 = [r4.readline().strip() for _ in range(4)]

            if not head4:
                break

            #get rc of index3 to see if it matches
            rc_barcode3 = reverse_compliment(barcode3)

            #get median quality of barcode reads
            qual_bar2 = sum((ord(x) - 33) for x in qual2) / len(qual2) 
            qual_bar3 = sum((ord(x) - 33) for x in qual3) / len(qual3) 
            
            #if either doesnt exist, or either quality is too low
            if rc_barcode3 not in BARCODE_MAP or barcode2 not in BARCODE_MAP or qual_bar2 <= cutoff or qual_bar3 <= cutoff:
                #writes to the unknown file
                write_reads('UK', head1, head4, seq1, seq4, qual1, qual4, barcode2, rc_barcode3)
                #statistics
                unknown_indexes += 1
            
            elif barcode2 == rc_barcode3:
                #write to the index file
                write_reads(BARCODE_MAP[barcode2], head1, head4, seq1, seq4, qual1, qual4, barcode2, rc_barcode3)
                #add statistics ....
                MATCH_DICT[f'{BARCODE_MAP[barcode2]}_{BARCODE_MAP[rc_barcode3]}'] += 1
                properly_matched += 1
            
            else:
                #Both barcodes are valid, but they dont match
                write_reads('Hop', head1, head4, seq1, seq4, qual1, qual4, barcode2, rc_barcode3)
                #statistics
                MATCH_DICT[f'{BARCODE_MAP[barcode2]}_{BARCODE_MAP[rc_barcode3]}'] += 1
                index_hopped += 1

    #Print the statistics
    with open(output_file,'w') as fo:
        fo.write('\t-----------------------STATISTIC REPORT FROM DEMULTIPLEXING-----------------\n\n')
        fo.write(f"Quality score cutoff: {cutoff}\n")
        fo.write(f'Total Reads: {unknown_indexes + properly_matched + index_hopped}\n\n')
        fo.write(f'Matched Reads Found: {properly_matched} - ({(properly_matched/(unknown_indexes+properly_matched+index_hopped))*100}%)\n')
        fo.write(f'Index Hopped Reads: {index_hopped} - ({(index_hopped/(unknown_indexes+properly_matched+index_hopped))*100}%)\n')
        fo.write(f'Unknown index Reads: {unknown_indexes} - ({(unknown_indexes/(unknown_indexes+properly_matched+index_hopped))*100}%)\n\n\n')

        #Statistics
        fo.write('\tIndividual Index pair percentages & total values:\n\n')
        for index_pair_key in sorted(MATCH_DICT, key=lambda k: MATCH_DICT[k], reverse=True):
            fo.write(f'{index_pair_key}: {MATCH_DICT[index_pair_key]} - ({MATCH_DICT[index_pair_key]/(index_hopped+properly_matched)}%)\n')



#grab the arguments used
def get_args():
    parser = argparse.ArgumentParser(description='args for Demultiplexing')
    parser.add_argument("-R1", "--read1", required=True, type=str, help="the read1 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R2", "--read2", required=True, type=str, help="the read2 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R3", "--read3", required=True, type=str, help="the read3 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R4", "--read4", required=True, type=str, help="the read4 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument('-c','--cutoff', required=False, default=0, type=str, help="The minimum of what the mean quality score of an index needs to be for it to be not thrown to unknown")
    parser.add_argument('-i','--index_file', required=True, type=str, help="The index file is what gives names and barcodes to each index that can be found")
    parser.add_argument('-o', '--output_dir', required=True, type=str, help='The name of a directory to store the output index files')
    parser.add_argument('-os', '--statistic_output', required=True, type=str, help='The text output file for the statistics')

    return parser.parse_args()
    
#testrun command
#./A3_the_third.py -R1 ../TEST-input_FASTQ/Test_R1.fq.gz -R2 ../TEST-input_FASTQ/Test_R2.fq.gz -R3 ../TEST-input_FASTQ/Test_R3.fq.gz -R4 ../TEST-input_FASTQ/Test_R4.fq.gz -i ../TEST-input_FASTQ/indexes.txt -o results -os stat_output.txt

if __name__ == "__main__":
    args = get_args()
    main(args.read1, args.read2, args.read3, args.read4, args.cutoff, args.index_file, args.output_dir, args.statistic_output)