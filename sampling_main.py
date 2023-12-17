import argparse
import random
import pandas as pd
import json

def read_fasta_sites(filename):
  df = pd.read_csv(filename, sep='\t')
  coordinate = pd.DataFrame()
  coordinate = df.iloc[:, :2]
  columns = ['chromosome', 'location(1 base)']
  coordinate.columns = columns
  return coordinate

def extract_rna_editing_seq(fasta_seq, header, editing_point, added_window, header_index):
  start_position = header_index.get(header)
  sequence = ''
  with open(fasta_seq, 'r') as f:
    f.seek(start_position)
    f.readline()

    for line in f:
      if line.startswith('>'):
        break
      sequence += line.strip()

  if editing_point - added_window < 1 or editing_point + added_window > len(sequence):
    return ''
  else: 
    return sequence[editing_point - added_window - 1 : editing_point + added_window - 1]

def extract_random_seq(fasta_seq, header, added_window, header_index):
  start_position = header_index.get(header)
  sequence = ''
  with open(fasta_seq, 'r') as f:
    f.seek(start_position)
    f.readline()

    for line in f:
      if line.startswith('>'):
        break
      sequence += line.strip()
    seq_len = len(sequence)

    if seq_len < 2 * added_window:
      return ''
    
    position = random.randint(0, seq_len - 1)
    while position - added_window < 1 or position + added_window > seq_len:
      position = random.randint(0, seq_len - 1)
    return sequence[position - added_window - 1 : position + added_window - 1]

def filter(seq, threshold = 0.1):
    n_count = seq.upper().count('N')
    return (n_count / len(seq)) < threshold

def extract_seqs(fasta_seq, n_samples, added_window, coordinate, header_index):
  random_samples = []
  edit_rich_samples = []

  num_seqs = len(coordinate)
  already_drawn = set()
  already_random = set()
  sample_count = 0
  
  while sample_count < n_samples:

    # random sampling for editing sites
    editing_site = random.randint(0, num_seqs - 1)
    while editing_site in already_drawn:
      editing_site = random.randint(0, num_seqs - 1)
    already_drawn.add(editing_site)

    # extract RNA editing neighboring sequence
    header = str(coordinate['chromosome'][editing_site])
    editing_point = int(coordinate['location(1 base)'][editing_site])
    seq_tmp_1 = extract_rna_editing_seq(fasta_seq, header, editing_point, added_window, header_index)
    if seq_tmp_1 == '':
      continue
    elif not filter(seq_tmp_1):
      continue

    # random sampling for non-editing sites
    rnd_position = random.randint(0, num_seqs - 1)
    while rnd_position in already_random:
      rnd_position = random.randint(0, num_seqs - 1)
    already_random.add(rnd_position)

    # extract random neighboring sequence
    header = str(coordinate['chromosome'][rnd_position])
    seq_tmp_2 = extract_random_seq(fasta_seq, header, added_window, header_index)
    if seq_tmp_2 == '':
      continue
    elif not filter(seq_tmp_2):
      continue

    edit_rich_samples.append(seq_tmp_1)
    random_samples.append(seq_tmp_2)
    sample_count += 1

    print('processing:', sample_count + 1, '/', n_samples)

  return edit_rich_samples, random_samples
    

def main():
    # DESCRIPTION OF EACH ARGUMENT + SAMPLE SCRIPT COMMAND:
    #
    # fsites:     The file of the rna editing sites 
    #             (make sure to convert to fasta form. You can just rename it to
    #             '<file name>.fasta' and it will autoconvert (took < 1 minute for me))
    #
    # fseq:       The file containing the sequence (also in fasta form)
    # 
    # w:          The added window for inverted repeat search. For example, say we have
    #             detected an RNA editing site as locations 900-1000, and we set
    #             the parameter w = 50, then we will search for inverted repeats inside
    #             the sequence located in locations 850-1050. If w=0 in this case, we
    #             only search for inverted repeats strictly inside the subsequence 
    #             located at 900-1000
    #
    #
    # nsamples:   The number of samples we wish to draw from random parts of the 
    #             sequence and rna editing sites. For example, if nsamples = 2, 
    #             we would return 4 subsequences total: 2 random subsequences and
    #             2 non-random subsequences where the length of the first random 
    #             sample is the same legnth as the first non-random sample, and the 
    #             second random sample is the same length as the second non-random
    #             sample
    #
    # out:        output file
    #
    #
    # Example:  
    #         python3 sampling.py -fsites testFile.txt -fseq test_data/rand100000.fasta -w 1 -nsamples 1 -out out.txt 
              # python sampling_updated.py -fsites /mnt/d/filtered_4_output.txt -fseq /mnt/d/GCF_genomic.fna -w 100 -nsamples 100 -s sample.txt -r random.txt -k index.json
    parser = argparse.ArgumentParser(
    description='Parse a sequence to uncover locations of inverted repeats.')
    parser.add_argument('-fsites', action="store", dest="fsites", type=str, required=True)
    parser.add_argument('-fseq', action="store", dest="fseq", type=str,
    required=True)
    parser.add_argument('-w', action="store", dest="w", type=int,
    required=True)
    parser.add_argument('-nsamples', action="store", dest="nsamples", type=int,
    required=True)
    parser.add_argument('-s', action="store", dest="s", type=str,
    required=True)
    parser.add_argument('-r', action="store", dest="r", type=str,
    required=True)
    parser.add_argument('-k', action="store", dest="k", type=str, required=True)

    args = parser.parse_args()
    fasta_edit_sites = args.fsites
    fasta_seq = args.fseq
    n_samples = args.nsamples
    output_file_1 = args.s
    output_file_2 = args.r
    added_window = args.w
    index_file = args.k

    coordinate = read_fasta_sites(fasta_edit_sites)
    # load json file from helper script
    with open(index_file, 'r') as f:
      header_index = json.load(f)
    non_random, random = extract_seqs(fasta_seq, n_samples, added_window, coordinate, header_index)

    with open(output_file_1, "w") as f:
        for i in range(len(non_random)):
          f.write(f">seq{i}\n")
          f.write(f"{non_random[i]}\n")
    
    with open(output_file_2, "w") as f:
        for i in range(len(random)):
          f.write(f">seq{i}\n")
          f.write(f"{random[i]}\n")

if __name__ == "__main__":
    main()