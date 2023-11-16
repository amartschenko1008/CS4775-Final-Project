import argparse
import numpy as np
import random

''' Reads in the sequences from the motif files.

Arguments:
    filename: which filename to read sequences from
Returns:
    output: list of sequences in motif file
'''
def read_fasta(filename):
    # TODO: update as needed to match structure needed for our test files (put test files inside src folder)
    with open(filename, "r") as f:
        output = []
        s = ""
        for l in f.readlines():
            if l.strip()[0] == ">":
                # skip the line that begins with ">"
                if s == "": continue
                output.append(s)
                s = ""
            # keep appending to the current sequence when the line doesn't begin
            # with ">"
            else: s += l.strip()
        output.append(s)
        return output

# TODO: implement necessary functions and helpers

def main():
  # TODO: parse test fasta files and print function results
  print("Hello World!")

if __name__ == '__main__':
    main()