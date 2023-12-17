import argparse
import json

'''
Input python sampling_helper.py -fseq -s 
'''

def find_header_index(fasta_file):
    header_index = {}
    with open(fasta_file, 'r') as f:
        while True:
            position = f.tell()
            line = f.readline() 
            if not line:  
                break
            if line.startswith('>'):  
                header = line.strip().lstrip('>').split()[0]
                header_index[header] = position
    return header_index  

def main():
    parser = argparse.ArgumentParser(
    description='Parse a sequence to uncover locations of inverted repeats.')
    parser.add_argument('-fseq', action="store", dest="fseq", type=str,
    required=True)
    parser.add_argument('-s', action="store", dest="s", type=str,
    required=True)

    args = parser.parse_args()
    fasta_seq = args.fseq
    output_file = args.s
    header_index = find_header_index(fasta_seq)

    with open(output_file, 'w') as f:
        json.dump(header_index, f)

if __name__ == "__main__":
    main()