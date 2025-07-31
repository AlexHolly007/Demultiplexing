import bioinfo
import argparse











def get_args():
    parser = argparse.ArgumentParser(description='args for Demultiplexing')
    parser.add_argument("-R1", "read1", required=True, type=str, help="the read1 file output from sequencing. THIS SHOULD BE A GZ FILE.")

if __name__ == "__main__":
    args = get_args()