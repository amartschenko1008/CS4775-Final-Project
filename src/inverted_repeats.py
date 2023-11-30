# # import argparse

# # def read_fasta(filename):
# #     sequences = []
# #     with open(filename, "r") as file:
# #         current_sequence = ""
# #         for line in file:
# #             line = line.strip()
# #             if line.startswith(">"):
# #                 if current_sequence:
# #                     sequences.append(current_sequence)
# #                     current_sequence = ""
# #             else:
# #                 current_sequence += line
# #         if current_sequence:
# #             sequences.append(current_sequence)
# #     return sequences

# # def generate_suffix_array_and_lcp(sequence):
# #     n = len(sequence)
# #     suffix_array = sorted(range(n), key=lambda i: sequence[i:])
# #     lcp = [0] * n

# #     for i in range(1, n):
# #         j = suffix_array[i - 1]
# #         while i < n and j < n and sequence[i] == sequence[j]:
# #             lcp[i] = max(lcp[i], n - j)
# #             i += 1
# #             j += 1

# #     return suffix_array, lcp

# # def find_mismatch_positions(sequence, center_position, gap_size):
# #     positions = []

# #     for gap_size in range(1, max_gap_size + 1):  # Iterate over gap sizes
# #         left_pos = center_position - gap_size
# #         right_pos = center_position + gap_size

# #         if (
# #             left_pos >= 0
# #             and right_pos < len(sequence)
# #             and sequence[left_pos] != sequence[right_pos]
# #         ):
# #             positions.append((left_pos, right_pos))

# #     return positions

# # def kangaroo_method(sequence, SA, LCP, left_pos, right_pos):
# #     common_extension = 0

# #     while (
# #         left_pos + common_extension < len(sequence)
# #         and right_pos + common_extension < len(sequence)
# #         and sequence[SA[left_pos] + common_extension] == sequence[SA[right_pos] + common_extension]
# #     ):
# #         common_extension += 1

# #     return common_extension

# # def find_inverted_repeats(sequence, SA, LCP, mismatch_positions, center_position, max_gap_size, max_size_range):
# #     inverted_repeats = []

# #     for mismatch_pair in mismatch_positions:
# #         left_pos, right_pos = mismatch_pair

# #         for gap_size in range(1, max_gap_size + 1):
# #             # Calculate the length of the potential inverted repeat
# #             length = 2 * gap_size + 1

# #             # Check if the length is within the specified range
# #             if length <= max_size_range:
# #                 # Find the common extension using the kangaroo method
# #                 common_extension = kangaroo_method(sequence, SA, LCP, left_pos, right_pos)

# #                 # Extend the inverted repeat to include the common extension
# #                 left_extension = center_position - gap_size
# #                 right_extension = center_position + gap_size + common_extension

# #                 # Check if the common extension allows further extension
# #                 while (
# #                     left_extension >= 0
# #                     and right_extension < len(sequence)
# #                     and sequence[left_extension] == sequence[right_extension - common_extension]
# #                 ):
# #                     left_extension -= 1
# #                     right_extension += 1

# #                 inverted_repeat = {
# #                     "start": left_extension + 1,
# #                     "end": right_extension - common_extension - 1,
# #                     "length": length + common_extension,
# #                 }
# #                 inverted_repeats.append(inverted_repeat)

# #     return inverted_repeats

# # def output_result(inverted_repeat):
# #     print("Inverted Repeat:", inverted_repeat)

# # def main(filename, max_gap_size, max_size_range):
# #     sequences = read_fasta(filename)

# #     for sequence in sequences:
# #         SA, LCP = generate_suffix_array_and_lcp(sequence)

# #         for center_position in range(len(sequence)):
# #             for gap_size in range(1, max_gap_size + 1):
# #                 mismatch_positions = find_mismatch_positions(sequence, center_position, gap_size)
# #                 inverted_repeats = find_inverted_repeats(sequence, SA, LCP, mismatch_positions, center_position, gap_size, max_size_range)

# #                 for inverted_repeat in inverted_repeats:
# #                     output_result(inverted_repeat)

# # if __name__ == '__main__':
# #     fasta_filename = "test.fasta"  # Replace with the actual filename of your FASTA file
# #     max_gap_size = 2
# #     max_size_range = 6
# #     main(fasta_filename, max_gap_size, max_size_range)
# import getopt
# import sys
# from collections import defaultdict
# from bisect import bisect_left, bisect_right
# from array import array

# # Helper function to print arrays
# def print_array(name, arr):
#     print(name + ":", " ".join(map(str, arr)))

# # Helper function to get the digit count of a number
# def get_digit_count(num):
#     return len(str(num))

# # Helper function to check if a file exists
# def file_exists(file_path):
#     try:
#         with open(file_path, 'r'):
#             return True
#     except FileNotFoundError:
#         return False

# # Function to insert into IUPAC map
# def iupac_map_insert(iupac_map, iupac_to_value, char, values):
#     iupac_map[char] = set(values)
#     index = len(iupac_to_value)
#     iupac_to_value[char] = index

# # Function to set value for matrix
# def set_value_for_matrix(matrix, size, i, j, value):
#     matrix[i * size + j] = value

# # Class to represent Match Matrix
# import getopt
# import sys
# from typing import Set, Tuple, List

# # Class to represent the match matrix and related functions
# class MatchMatrix:
#     match_matrix = []
#     IUPAC_map_count = 0
#     IUPAC_to_value = []

#     @staticmethod
#     def match(char1, char2):
#         i = MatchMatrix.IUPAC_to_value[ord(char1)]
#         j = MatchMatrix.IUPAC_to_value[ord(char2)]
#         return MatchMatrix.match_matrix[i][j]

# # Function to parse command line arguments
# def parse_arguments(argv):
#     input_file = "input.fasta"
#     seq_name = "seq0"
#     min_len = 10
#     max_len = 100
#     max_gap = 100
#     mismatches = 0
#     output_file = "IUPACpal.out"

#     try:
#         opts, args = getopt.getopt(argv, "f:s:m:M:g:x:o:")
#     except getopt.GetoptError:
#         print("Invalid arguments")
#         sys.exit(2)

#     for opt, arg in opts:
#         if opt == '-f':
#             input_file = arg
#         elif opt == '-s':
#             seq_name = arg
#         elif opt == '-m':
#             min_len = int(arg)
#         elif opt == '-M':
#             max_len = int(arg)
#         elif opt == '-g':
#             max_gap = int(arg)
#         elif opt == '-x':
#             mismatches = int(arg)
#         elif opt == '-o':
#             output_file = arg

#     return input_file, seq_name, min_len, max_len, max_gap, mismatches, output_file

# # Function to check if a file exists
# def file_exists(file_path):
#     try:
#         with open(file_path):
#             pass
#         return True
#     except FileNotFoundError:
#         return False

# # Function to read input file and extract the specified sequence
# def extract_sequence(input_file, seq_name):
#     found_seq = False
#     name = ""
#     contents = ""

#     with open(input_file, 'r') as file:
#         for line in file:
#             line = line.strip()
#             if not found_seq:
#                 if len(line) > 0 and line[0] == '>':
#                     parts = line[1:].split()
#                     name = parts[0]
#                 if name == seq_name:
#                     found_seq = True
#                 name = ""
#             else:
#                 if line and line[0] not in {' ', '>', ';'}:
#                     contents += line
#                 else:
#                     break

#     return contents

# # Main function
# def main(argv):
#     input_file, seq_name, min_len, max_len, max_gap, mismatches, output_file = parse_arguments(argv)

#     if not file_exists(input_file):
#         print(f"Error: File '{input_file}' not found.")
#         sys.exit(-1)

#     sequence = extract_sequence(input_file, seq_name)

#     if not sequence:
#         print(f"Error: Sequence '{seq_name}' not found in file '{input_file}'.")
#         sys.exit(-1)

#     n = len(sequence)

#     # Verify arguments are valid
#     if min_len < 2 or min_len > 2**31:
#         print("Error: min_len must be between 2 and 2^31.")
#         sys.exit(-1)

#     if max_len < 0 or max_len > 2**31:
#         print("Error: max_len must be between 0 and 2^31.")
#         sys.exit(-1)

#     if max_gap < 0 or max_gap > 2**31:
#         print("Error: max_gap must be between 0 and 2^31.")
#         sys.exit(-1)

#     if mismatches < 0 or mismatches > 2**31:
#         print("Error: mismatches must be between 0 and 2^31.")
#         sys.exit(-1)

#     if min_len >= n or max_len < min_len or max_gap >= n or mismatches >= n or mismatches >= min_len:
#         print("Error: Invalid arguments with respect to sequence length.")
#         sys.exit(-1)

#     # Optionally display user-given options
#     if True:
#         print("\ninput_file:", input_file)
#         print("seq_name:", seq_name)

#         # Optionally print sequence data
#         if False:
#             print("sequence:", sequence)

#         print("min_len:", min_len)
#         print("max_len:", max_len)
#         print("max_gap:", max_gap)
#         print("mismatches:", mismatches)
#         print("output_file:", output_file, "\n")

# def build_match_matrix(iupac_map):
#     iupac_map_count = len(iupac_map)
#     match_matrix = [[False] * iupac_map_count for _ in range(iupac_map_count)]

#     for i, (char1, set1) in enumerate(iupac_map.items()):
#         for j, (char2, set2) in enumerate(iupac_map.items()):
#             match = any(c1 == c2 for c1 in set1 for c2 in set2)
#             match_matrix[i][j] = match

#     return match_matrix

# # Build complement array
# def build_complement_array():
#     complement = {
#         'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'u': 'a',
#         'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
#         'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b',
#         'n': 'n', '*': 'n', '-': 'n'
#     }
#     return complement

# # Build IUPAC_map and IUPAC_to_value
# def build_iupac_map():
#     iupac_map = {
#         'a': {'a'}, 'c': {'c'}, 'g': {'g'}, 't': {'t'},
#         'u': {'t'}, 'r': {'a', 'g'}, 'y': {'c', 't'},
#         's': {'g', 'c'}, 'w': {'a', 't'}, 'k': {'g', 't'},
#         'm': {'a', 'c'}, 'b': {'c', 'g', 't'}, 'd': {'a', 'g', 't'},
#         'h': {'a', 'c', 't'}, 'v': {'a', 'c', 'g'},
#         'n': {'a', 'c', 'g', 't'}, '*': {'a', 'c', 'g', 't'},
#         '-': {'a', 'c', 'g', 't'}, '$': {'$'}, '#': {'#'}
#     }
#     iupac_to_value = {char: i for i, char in enumerate(iupac_map)}

#     return iupac_map, iupac_to_value

# # Optionally print match matrix
# def print_match_matrix(match_matrix, iupac_map):
#     letters = list(iupac_map.keys())

#     print("Match Matrix:")
#     print("  " + " ".join(letters))
    
#     for i, row in enumerate(match_matrix):
#         print(letters[i] + " " + " ".join(map(str, row)))

#     print("\n\n")

# def calculate_suffix_array(s):
#     n = len(s)
#     sa = list(range(n))
#     sa.sort(key=lambda i: s[i:])
#     return sa

# # Placeholder for inverse suffix array calculation
# def calculate_inverse_suffix_array(sa):
#     n = len(sa)
#     inv_sa = [0] * n
#     for i in range(n):
#         inv_sa[sa[i]] = i
#     return inv_sa

# # Placeholder for Kasai's algorithm for LCP calculation
# def calculate_lcp(s, sa, inv_sa):
#     n = len(s)
#     rank = [0] * n
#     k = 0
#     lcp = [0] * n

#     for i in range(n):
#         rank[sa[i]] = i

#     for i in range(n):
#         if rank[i] == n - 1:
#             k = 0
#             continue

#         j = sa[rank[i] + 1]
#         while i + k < n and j + k < n and s[i + k] == s[j + k]:
#             k += 1

#         lcp[rank[i]] = k
#         if k > 0:
#             k -= 1

#     return lcp

# # Placeholder for sparse table RMQ of LCP
# def calculate_rmq_lcp(lcp):
#     n = len(lcp)
#     log_n = (n - 1).bit_length()
#     sparse_table = [[0] * log_n for _ in range(n)]

#     for i in range(n):
#         sparse_table[i][0] = lcp[i]

#     for j in range(1, log_n):
#         for i in range(n - (1 << j) + 1):
#             sparse_table[i][j] = min(sparse_table[i][j - 1], sparse_table[i + (1 << (j - 1))][j - 1])

#     def query(i, j):
#         k = (j - i).bit_length() - 1
#         return min(sparse_table[i][k], sparse_table[j - (1 << k) + 1][k])

#     return query

# def test_calculate_suffix_array():
#     s = "banana"
#     expected_result = [5, 3, 1, 0, 4, 2]
#     assert calculate_suffix_array(s) == expected_result

# def calculate_inverse_suffix_array(sa):
#     n = len(sa)
#     inv_sa = [0] * n

#     for i in range(n):
#         inv_sa[sa[i]] = i

#     return inv_sa

# def test_calculate_inverse_suffix_array():
#     sa = [5, 3, 1, 0, 4, 2]
#     expected_result = [3, 5, 2, 1, 4, 0]
#     assert calculate_inverse_suffix_array(sa) == expected_result

# def test_calculate_lcp():
#     s = "banana"
#     sa = [5, 3, 1, 0, 4, 2]
#     inv_sa = [3, 5, 2, 1, 4, 0]
#     expected_result = [0, 1, 3, 0, 0, 2]
#     assert calculate_lcp(s, sa, inv_sa) == expected_result

# def test_calculate_rmq_lcp():
#     lcp = [0, 1, 3, 0, 0, 2]
#     rmq = calculate_rmq_lcp(lcp)
#     assert rmq(1, 3) == 1
#     assert rmq(2, 5) == 0
#     assert rmq(0, 4) == 0

# if __name__ == "__main__":
#     test_calculate_suffix_array()
#     test_calculate_inverse_suffix_array()
#     test_calculate_lcp()
#     test_calculate_rmq_lcp()
#     main(sys.argv[1:])
# import re

# def read_fasta(file_path):
#     """
#     Read a FASTA file and return the DNA sequence.
#     """
#     with open(file_path, 'r') as file:
#         lines = file.readlines()

#     # Assuming that the sequence is in a single line after the header
#     sequence = ''.join(line.strip() for line in lines[1:])
#     return sequence

# def find_inverted_repeats(dna_sequence, min_length=2, max_length=None, max_mismatches=0):
#     inverted_repeats = []

#     if max_length is None:
#         max_length = len(dna_sequence)

#     for length in range(min_length, max_length + 1):
#         for i in range(len(dna_sequence) - length + 1):
#             subsequence = dna_sequence[i:i + length]
#             reverse_complement = get_reverse_complement(subsequence)
#             mismatches = sum(base1 != base2 for base1, base2 in zip(subsequence, reverse_complement))
    
#             if reverse_complement in dna_sequence[i+length:] and subsequence not in inverted_repeats:
#                 starts = find_all_substrings(dna_sequence, reverse_complementi)
#                 for s in starts:
#                     if i + length - 1 < s: 
#                         inverted_repeats.append((i, i + length - 1, s, s + length - 1))
    
#     return inverted_repeats

# def find_all_substrings(string, substring):
#     # Initialize an empty list to store the 
#     # indices of all occurrences of the substring.
#     indices = []
#     # Set the starting index i to 0.
#     i = 0
#     # Use a while loop to keep searching for 
#     # the substring in the string.
#     while i < len(string):
#         # Use the find() method to find the first 
#         #occurrence of the substring in the string
#         j = string.find(substring, i)
#         # If find() returns -1, it means that there  
#         # are no more occurrences of the substring in 
#         # the string, so break out of the loop.
#         if j == -1:
#             break
#         indices.append(j)
#         i = j + len(substring)
#     # Return the list of indices.
#     return indices

# def get_reverse_complement(sequence):
#     complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
#     reverse_complement = ''.join(complement_dict[base] for base in reversed(sequence))
#     return reverse_complement

# def print_inverted_repeats(inverted_repeats, dna_sequence):
#     print("Palindromes:")
#     for idx, (start_first, end_first, start_sec, end_sec) in enumerate(inverted_repeats, start=1):
#         first = dna_sequence[start_first: end_first + 1]
#         sec = dna_sequence[start_sec: end_sec + 1][::-1]
#         mismatches = "".join("|" if base1 != base2 else " " for base1, base2 in zip(first, sec))

        
#         print(f"{start_first:2d} {first:<2s} {end_first:2d}")
#         print(f"   {mismatches}")
#         print(f"{end_sec:2d} {sec:<2s} {start_sec:2d}\n")

# def main():
#     fasta_file_path = 'test.fasta'
#     dna_sequence = read_fasta(fasta_file_path)
#     min_palindrome_half_length = 3
#     max_palindrome_half_length = 12
#     max_mismatch = 1

#     inverted_repeats = find_inverted_repeats(dna_sequence, min_length=min_palindrome_half_length, max_length=max_palindrome_half_length, max_mismatches=0)
#     print(inverted_repeats)
#     if inverted_repeats:
#         print_inverted_repeats(inverted_repeats, dna_sequence)
#     else:
#         print("No inverted repeats found.")


# if __name__ == "__main__":
#     main()
#------------------------------------------------------------------------------------------
# def find_inverted_repeats(dna_sequence, min_len, max_len, max_gap, mismatches):
#     repeats = []

#     for length in range(min_len, max_len + 1):
#         for start1 in range(len(dna_sequence) - length + 1):
#             end1 = start1 + length
#             for start2 in range(len(dna_sequence) - length + 1):
#                 end2 = start2 + length

#                 # Check for valid gap
#                 if abs(start2 - end1) <= max_gap:
#                     sequence1 = dna_sequence[start1:end1]
#                     sequence2 = reverse_complement(dna_sequence[start2:end2])

#                     # Check for mismatches
#                     mismatch_count = sum(1 for base1, base2 in zip(sequence1, sequence2) if base1 != base2)
#                     if mismatch_count <= mismatches:
#                         repeats.append((start1, end1, start2, end2))

#     return repeats

# def reverse_complement(sequence):
#     complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
#     return ''.join(complement_dict[base] for base in reversed(sequence))

# # Example usage:
# dna_sequence = "ATGCGATCGATCGATCGATCGTAGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
# #dna_sequence = "ATGCAT"
# min_len = 3
# max_len = 7
# max_gap = 2
# mismatches = 0

# result = find_inverted_repeats(dna_sequence, min_len, max_len, max_gap, mismatches)
# print("Inverted repeats found:")
# for repeat in result:
#     print(f"Repeat 1: {repeat[0]}-{repeat[1]}, {dna_sequence[repeat[0]: repeat[1] + 1] if repeat[1] != repeat[2] else dna_sequence[repeat[0]: repeat[1]]}, Repeat 2: {repeat[2]}-{repeat[3]}, {dna_sequence[repeat[2]: repeat[3] + 1]}")
#------------------------------------------------------------------------------------------
import re
def readFASTA(filename): #function to extract sequence from FASTA
    fo = open(filename)
    lines = fo.readlines() #read all lines in file
    seq = '' #empty string to store extracted sequence
    
    for line in lines:
        if line.startswith('>'): #ignore the header line which starts with '>'
            pass
        else:
            line = re.sub('\n','',line) #replace newlines with nothing
            seq += line #update seq variable with lines

    return seq #return extracted sequence from fasta

def readGB(filename): #function to extract sequence from Genbank
    gb_seq = '^\s+\d+\s+(([a-z_]+\s*)+)' #regular expression (re) pattern to match sequence part    #\s+\d+\s+([a-zA-Z\s]+)+
    fo = open(filename)
    lines = fo.readlines() #read all lines in file
    seq = '' #empty string to store extracted sequence
    
    for line in lines:
        sequenceline = re.search(gb_seq,line) #use re pattern to search for sequence in file
        if sequenceline: #if pattern found
            grp1 = sequenceline.group(1) #for the first subgroup
            sequenceline1 = re.sub('[\s]','',grp1) #replace whitespace with nothing
            seq += sequenceline1 #update seq variable with lines

    return seq #return extracted sequence from genbank

def reverseComplement(seq): #function to commpute reverse complement
    Base = {'a':'t','t':'a','g':'c','c':'g', '_':'_'} #dictionary to store base pair for A,G,C,T and spacer region
    ComplementSeq = "" #empty string to store complement sequence

    for i in range(0, len(seq)): #for every base in the whole original sequence length
        pair = seq[i] #new variable to store each base in original sequence
        ComplementSeq = ComplementSeq + Base[pair] #concatenate the complement base pair together using values from dictionary
    reverseComplementSeq = ComplementSeq[::-1] #reverse the ComplementSeq
    
    return reverseComplementSeq #return the reverse complement sequence

def CommonSequence(seq,revCom,minLength): #function to extract common sequence between the original and reverse complement
    seqLength = len(seq) #compute the sequence length
    commonSequence = [] #empty list to store common sequence

    for i in range(seqLength,minLength-1,-1): #loop from the reverse of sequence
    #until min palindrome length
        for k in range(seqLength-i+1): #loop in the length of short sequence
            if (seq[k:i+k] in revCom): #true if base present in reverse complement
                flag = 1
                for m in range(len(commonSequence)): #loop in the length of list
                    if seq[k:i+k] in commonSequence[m]: #if base is already present in list
                        flag = 0 
                        break #break the loop

                if flag == 1: #if base is not already present in list
                    commonSequence.append(seq[k:i+k]) #add base to the list

    if len(commonSequence): #true if list is not empty
        return(commonSequence) #return list that contains common sequences
    else: #false if list is empty
        pass

def AllPalindrome(allMatches): #function to find all palindromes
    allPalindrome = [] #empty list to store all palindromes
    for sequence in allMatches: #for every sequence in all the common sequence
        #check if that particular sequence is equivalent to its reverse complement (means its a palindrome)
        #and if that sequence does not exist in the list already
        if sequence == reverseComplement(sequence) and sequence not in allPalindrome: #true
            allPalindrome.append(sequence) #add that sequence to the list

    return allPalindrome #return all the palindromes in the whole sequence

def NormalPalindrome(allPalindrome): #function to find normal palindromes (without spacer region)
    normalPalindrome = [] #empty list to store normal palindromes
    for sequence in allPalindrome: #for every sequence in all the palindromes
        if '_' not in sequence: #filter out palindromes that doesnt contain '_'
            normalPalindrome.append(sequence) #add that palindrome to the list

    if len(normalPalindrome): #print out all the normal palindromes if available
        normalPalindrome = ', '.join(normalPalindrome) #convert normal palindrome list to string for output
        print("\nNormal palindromes (non-repeating): \n",normalPalindrome,"\n") 
    else:
         print("There are no normal palindromes that can be detected.\n")

def SpacerPalindrome(allPalindrome): #function to find spacer palindromes
    allPalindrome = ' '.join(allPalindrome) #convert list of all palindromes to string for re
    spacerPalindrome = re.findall(r'[agct]+_+[agct]+',allPalindrome) #find all spacer palindromes using re

    if len(spacerPalindrome): #print out all the spacer palindromes if available
        spacerPalindrome = ', '.join(spacerPalindrome) #convert spacer palindrome list to string for output
        print("Reverse-complement non-repeating palindromes with an intervening spacer region: \n",spacerPalindrome,"\n")
    else:
        print("There is no reverse-complement non-repeating palindromes with an intervening spacer region that can be detected.\n")

seq = readFASTA('rand1000.fasta') #call function fileInput()
minLength = 3
revCom = reverseComplement(seq) #call function reverseComplement()
allMatches = CommonSequence(seq,revCom,minLength) #call function CommonSequence()
allPalindrome = AllPalindrome(allMatches) #call function AllPalindrome()
NormalPalindrome(allPalindrome) #call function NormalPalindrome()
SpacerPalindrome(allPalindrome) #call function SpacerPalindrome()