from Bio import SeqIO
import numpy as np

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-1):
    n = len(seq1)
    m = len(seq2)

    score_matrix = np.zeros((n + 1, m + 1))
    traceback_matrix = np.full((n + 1, m + 1), '')

    for i in range(1, n + 1):
        score_matrix[i][0] = gap_penalty * i
        traceback_matrix[i][0] = 'U'

    for j in range(1, m + 1):
        score_matrix[0][j] = gap_penalty * j
        traceback_matrix[0][j] = 'L'

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            diag = score_matrix[i - 1][j - 1] + match
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(diag, up, left)

            if score_matrix[i][j] == diag:
                traceback_matrix[i][j] = 'D'  #diagonal
            elif score_matrix[i][j] == up:
                traceback_matrix[i][j] = 'U'  #up
            else:
                traceback_matrix[i][j] = 'L'  #left

    align1, align2 = '', ''
    i, j = n, m

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 'D':
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'U':
            align1 += seq1[i - 1]
            align2 += '-'
            i -= 1
        else:
            align1 += '-'
            align2 += seq2[j - 1]
            j -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]

    #print(score_matrix)
    #print(traceback_matrix)

    return score_matrix[n][m], align1, align2

with open('/home/bumblebee/studia/semestr5/strukturalna/seqs.fa', 'r') as file:
    sequences = list(SeqIO.parse(file, "fasta"))

if len(sequences) < 2:
    print("Plik FASTA musi zawierać dwie sekwencje.")
else:
    seq1 = str(sequences[0].seq)
    seq2 = str(sequences[1].seq)

    score, alignment1, alignment2 = needleman_wunsch(seq1, seq2)

    print(f"\nWynik dopasowania: {score}")
    print(f"Alignment 1: {alignment1}")
    print(f"Alignment 2: {alignment2}")
