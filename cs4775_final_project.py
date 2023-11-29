#!/usr/bin/env python3


import argparse
import json

'''Computes the score and alignment of two strings.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def NeedlemanWunsch(x, y, s, d):
    a = len(x)
    b = len(y)
    ''' Recurrence matrix, redefine/use as necessary. '''
    m = [[0] * (b + 1) for i in range(0, a + 1)]
    ''' Traceback matrix, redefine/use as necessary. '''
    t = [[0] * (b + 1) for i in range(0, a + 1)]

    for i in range(1, a + 1):
        m[i][0] , t[i][0] = m[i - 1][0] - d, 1
    for j in range(1, b + 1):
        m[0][j] , t[0][j] = m[0][j - 1] - d, 2
    
    for i in range(1, a + 1):
        for j in range(1, b + 1):
            row = m[i - 1][j] - d
            col = m[i][j - 1] - d
            diag = m[i - 1][j - 1] + s[x[i - 1]][y[j - 1]]
            m[i][j] = max(row, col, diag)
            if m[i][j] == row:
                t[i][j] = 1 
            elif m[i][j] == col:
                t[i][j] = 2
            else:
                t[i][j] = 0
    i = len(x)
    j = len(y)
    a_x = []
    a_y = []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and t[i][j] == 0:
            a_x.append(x[i - 1])
            a_y.append(y[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and t[i][j] == 1:
            a_x.append(x[i - 1])
            a_y.append("-")
            i -= 1
        else:
            a_x.append("-")
            a_y.append(y[j-1])
            j -= 1

    a_x = ''.join(a_x[::-1])
    a_y = ''.join(a_y[::-1])
    return m[a][b], (a_x, a_y)


'''Prints two aligned sequences formatted for convenient inspection.
Arguments:
    a_x: the first sequence aligned
    a_y: the second sequence aligned
Outputs:
    Prints aligned sequences (80 characters per line) to console
'''
def print_alignment(a_x, a_y):
    assert len(a_x) == len(a_y), "Sequence alignment lengths must be the same."
    for i in range(1 + (len(a_x) // 80)):
        start = i * 80
        end = (i + 1) * 80
        print(a_x[start:end])
        print(a_y[start:end])
        print()


'''def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for two sequences with a linear gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = NeedlemanWunsch(x, y, s, d)
    print("Alignment:")
    print_alignment(a_x, a_y)
    print("Score: " + str(score))


if __name__ == "__main__":
    main()'''


#Smith-Waterman algorithm for local alignment in python
'''Computes the score and local alignment of two strings.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    d: the gap opening/extension penalty
    s: the match score
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def SmithWaterman(x, y, d, s):
    a = len(x)
    b = len(y)
    m = [[0] * (b + 1) for i in range(0, a + 1)]
    score = 0
    optLoc = (0,0)
    for i in range(1, a + 1):
         for j in range(1, b + 1):
              row = m[i - 1][j] - d
              col = m[i][j - 1] - d
              diag = m[i - 1][j - 1] + (s if x[i-1] == y[j-1] else 0)
              m[i][j] = max(row, col, diag, 0)
              if m[i][j] >= score:
                  score = m[i][j]
                  optLoc = (i,j)
    seq = ''
    i = optLoc[0]
    j = optLoc[1]
    while(i > 0 or j > 0):
        if m[i][j] == 0:
            break
        elif m[i][j] == m[i-1][j] - d:
            i -= 1
            seq += x[i]
        elif m[i][j] == m[i][j-1] - d:
            j -= 1
            seq += y[j]
        else:
            i -= 1
            j -= 1
            seq += x[i]
    return seq[::-1]

def removeDuplicates(fastq):
    fastdict = {}
    for read in fastq:
        if read not in fastdict:
            fastdict[read] = 1
        else:
            fastdict[read] += 1
    return fastdict.keys()    

