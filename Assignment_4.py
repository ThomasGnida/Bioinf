"""
Assignment in Sequence Alignment

Title: A comprehensive implementation, evaluation and comparison of sequence alignment methods

Tasks
A. Implementations — Write modular, documented Python implementations with unit tests.
B. Experiments — Compare accuracy, runtime, memory, and parameter sensitivity.
C. Analysis — Complexity analysis, biological use cases, discussion of strengths/weaknesses.

"""

import numpy as np
import Bio
import pandas as pd


def benchmark_algorithm(sequence_array, algorithm_method):
    runtime = 0
    memory = 0


    return accuracy, runtime, memory, parameter_sensitivity, complexity


def score(a, b, match=1, mismatch=-1):
    """
    Score for matching characters.
    :param a: Character from seq1
    :param b: Character from seq2
    :param match: Score for match
    :param mismatch: Score for mismatch
    :return: Score
    """
    return match if a == b else mismatch

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Global alignment using Needleman-Wunsch algorithm. NW creates a scoring matrix by determining
    the "best" solution for each next step. Once the matrix is calculated the algorithm reconstructs
    the aligned sequence with the highest score.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Penalty for mismatch
    :param gap: Penalty for gap
    :return: Aligned seq1, aligned seq2, score
    """
    # Initialize the scoring matrix
    n = len(seq1) + 1
    m = len(seq2) + 1
    matrix = np.zeros(n, m)

    # Initialize the first row and column with gaps to account for gaps at the start of alignment
    for i in range(n):
        matrix[i][0] = i * gap
    for j in range(m):
        matrix[0][j] = j * gap

    # Fill the scoring matrix
    for i in range(1, n):
        for j in range(1, m):
            match_score = matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
            delete = matrix[i - 1][j] + gap
            insert = matrix[i][j - 1] + gap
            matrix[i][j] = max(match_score, delete, insert)

    #Traverse the created scoring matrix backwards to get optimal alignment (bottom-right to top-left)
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = len(seq1), len(seq2)

    while i > 0 and j > 0:
        current_score = matrix[i][j]
        if current_score == matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == matrix[i - 1][j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Add remaining gaps if necessary
    while i > 0:
        aligned_seq1.append(seq1[i - 1])
        aligned_seq2.append('-')
        i -= 1
    while j > 0:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[j - 1])
        j -= 1

    # Reverse the aligned sequences
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, matrix[-1][-1]

def smith_waterman(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Local alignment using Smith-Waterman algorithm.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap
    :return: Aligned seq1, aligned seq2, score
    """
    n = len(seq1) + 1
    m = len(seq2) + 1
    matrix = np.zeros((n, m))
    max_score = 0
    max_pos = None

    # Fill the scoring matrix
    for i in range(1, n):
        for j in range(1, m):
            match_score = matrix[i - 1][j - 1] + score(seq1[i-1], seq2(j-1))
            delete = matrix[i - 1][j] + gap
            insert = matrix[i][j - 1] + gap
            matrix[i][j] = max(0, match_score, delete, insert)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)
                #If the maximum score is highest in place it marks the optimal local alignment

    # Traceback
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos

    while matrix[i][j] > 0:
        if matrix[i][j] == matrix[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i - 1][j] + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse the aligned sequences
    aligned_seq1 = ''.join(reversed(aligned_seq1))
    aligned_seq2 = ''.join(reversed(aligned_seq2))

return aligned_seq1, aligned_seq2, max_score

def gotoh(seq1, seq2, match=1, mismatch=-1, gap_open=-5, gap_extend=-1):
    return
def banded_dp(seq1, seq2, k, match=1, mismatch=-1, gap=-2):
    return

def blast_seed_extend(seq1, seq2, k=11, match=1, mismatch=-1, gap=-2):
    return

def greedy():
    return

def progressive_alignment():
    return

def iterative_refinement():
    return

def profile_hmm_alignment():
    return