import bioinfo
import argparse

C_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

def reverse_compliment(sequence: str) -> str:

    return str([C_DICT[x] for x in sequence])[::-1]











def get_args():
    parser = argparse.ArgumentParser(description='args for Demultiplexing')
    parser.add_argument("-R1", "read1", required=True, type=str, help="the read1 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R2", "read2", required=True, type=str, help="the read2 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R3", "read3", required=True, type=str, help="the read3 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument("-R4", "read4", required=True, type=str, help="the read4 file output from sequencing. THIS SHOULD BE A GZ FILE.")
    parser.add_argument('-c','--cutoff', required=False, type=str, help="The minimum of what the mean quality score of an index needs to be for it to be not thrown to unknown")
    parser.add_argument('-i','--index_file', required=True, type=str, help="The index file is what gives names and barcodes to each index that can be found")

    return parser.parse_args()
    
#

if __name__ == "__main__":
    args = get_args()