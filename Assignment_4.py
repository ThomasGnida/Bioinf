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
    matrix = np.zeros((n, m))

    # Initialize the first row and column with gap penalties
    # This accounts for gaps at the beginning of the sequences
    for i in range(n):
        matrix[i, 0] = i * gap  # Gaps for seq1
    for j in range(m):
        matrix[0, j] = j * gap  # Gaps for seq2

    # Fill the scoring matrix
    for i in range(1, n):
        for j in range(1, m):
            # Calculate the possible scores for this position
            match_score = matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
            delete = matrix[i - 1, j] + gap  # Gap in seq2
            insert = matrix[i, j - 1] + gap  # Gap in seq1
            matrix[i, j] = max(match_score, delete, insert)

    #Traverse the created scoring matrix backwards to get optimal alignment (bottom-right to top-left)
    aligned_seq1 = []
    aligned_seq2 = []
    i = len(seq1)
    j = len(seq2)

    while i > 0 and j > 0:
        score_current = matrix[i, j]
        score_diag = matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
        if score_current == score_diag:
            # Match or mismatch
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif score_current == matrix[i - 1, j] + gap:
            # Gap in seq2
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            # Gap in seq1
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Add any remaining characters with gaps
    while i > 0:
        aligned_seq1.append(seq1[i - 1])
        aligned_seq2.append('-')
        i -= 1
    while j > 0:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[j - 1])
        j -= 1

    # Since we built the sequences backwards, reverse them
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    aligned_seq1 = ''.join(aligned_seq1)
    aligned_seq2 = ''.join(aligned_seq2)

    return aligned_seq1, aligned_seq2, matrix[n-1, m-1]

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
    # Set up the matrix for local alignment
    # Similar to NW but we allow negative scores to be set to zero
    n = len(seq1) + 1
    m = len(seq2) + 1
    matrix = np.zeros((n, m))
    max_score = 0
    max_pos = None

    # Fill the matrix, resetting to zero if negative
    for i in range(1, n):
        for j in range(1, m):
            # Calculate the possible scores
            match_score = matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
            delete = matrix[i - 1, j] + gap  # Gap in seq2
            insert = matrix[i, j - 1] + gap  # Gap in seq1
            matrix[i, j] = max(0, match_score, delete, insert)

            # Keep track of the highest score position
            if matrix[i, j] > max_score:
                max_score = matrix[i, j]
                max_pos = (i, j)
                # This marks the starting point of the optimal local alignment

    # If no positive score, return empty
    if max_score == 0:
        return "", "", 0

    # Traceback
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos

    while matrix[i, j] > 0:
        score_current = matrix[i, j]
        score_diag = matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
        if score_current == score_diag:
            # Match or mismatch
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif score_current == matrix[i - 1, j] + gap:
            # Gap in seq2
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            # Gap in seq1
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    # Reverse the sequences
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    aligned_seq1 = ''.join(aligned_seq1)
    aligned_seq2 = ''.join(aligned_seq2)

    return aligned_seq1, aligned_seq2, max_score

def gotoh(seq1, seq2, match=1, mismatch=-1, gap_open=-5, gap_extend=-1):
    """
    Global alignment using Gotoh algorithm with affine gap penalties.
    Uses three matrices: M for matches, Ix for insertions (gaps in seq2), Iy for deletions (gaps in seq1).
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Penalty for mismatch
    :param gap_open: Penalty for opening a gap
    :param gap_extend: Penalty for extending a gap
    :return: Aligned seq1, aligned seq2, score
    """
    # Set up the Matrices
    n = len(seq1) + 1
    m = len(seq2) + 1
    M = np.full((n, m), -np.inf)
    Ix = np.full((n, m), -np.inf)
    Iy = np.full((n, m), -np.inf)
    M[0, 0] = 0

    # Initialize first row (gaps in seq1, so Iy)
    for j in range(1, m):
        Iy[0, j] = gap_open + (j - 1) * gap_extend

    # Initialize first column (gaps in seq2, so Ix)
    for i in range(1, n):
        Ix[i, 0] = gap_open + (i - 1) * gap_extend

    # Fill the matrices
    for i in range(1, n):
        for j in range(1, m):
            # Score for match/mismatch
            s = score(seq1[i - 1], seq2[j - 1], match, mismatch)
            # M matrix: from diagonal
            M[i, j] = max(M[i - 1, j - 1] + s, Ix[i - 1, j - 1] + s, Iy[i - 1, j - 1] + s)
            # Ix matrix: insertions (gaps in seq2)
            Ix[i, j] = max(M[i - 1, j] + gap_open, Ix[i - 1, j] + gap_extend)
            # Iy matrix: deletions (gaps in seq1)
            Iy[i, j] = max(M[i, j - 1] + gap_open, Iy[i, j - 1] + gap_extend)

    # The final score is the max at the end
    final_score = max(M[n - 1, m - 1], Ix[n - 1, m - 1], Iy[n - 1, m - 1])

    # Traceback to reconstruct alignment
    aligned_seq1 = []
    aligned_seq2 = []
    i = n - 1
    j = m - 1

    # Determine which matrix to start from
    state = 0
    if final_score == M[i, j]:
        state = 0
    elif final_score == Ix[i, j]:
        state = 1
    else:
        state = 2
    while i > 0 or j > 0:
        if state == 0:  # M
            # Determine the source before moving
            s = score(seq1[i - 1], seq2[j - 1], match, mismatch)
            if M[i, j] == M[i - 1, j - 1] + s:
                next_state = 0
            elif M[i, j] == Ix[i - 1, j - 1] + s:
                next_state = 1
            else:
                next_state = 2
            # Append characters
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
            state = next_state
        elif state == 1:  # Ix
            # Determine source for Ix
            if Ix[i, j] == M[i - 1, j] + gap_open:
                next_state = 0
            else:
                next_state = 1
            # Append
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
            state = next_state
        else:  # Iy
            # Determine source for Iy
            if Iy[i, j] == M[i, j - 1] + gap_open:
                next_state = 0
            else:
                next_state = 2
            # Append
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
            state = next_state

    # Reverse the sequences
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    aligned_seq1 = ''.join(aligned_seq1)
    aligned_seq2 = ''.join(aligned_seq2)

    return aligned_seq1, aligned_seq2, final_score

def banded_dp(seq1, seq2, k, match=1, mismatch=-1, gap=-2):
    """
    Banded dynamic programming for global alignment, restricting to band around diagonal.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param k: Band width
    :param match: Score for match
    :param mismatch: Penalty for mismatch
    :param gap: Penalty for gap
    :return: Aligned seq1, aligned seq2, score
    """
    n = len(seq1) + 1
    m = len(seq2) + 1
    matrix = np.full((n, m), -np.inf)
    # Initialize the band
    for i in range(n):
        for j in range(m):
            if abs(i - j) <= k:
                if i == 0 and j == 0:
                    matrix[i, j] = 0
                elif i == 0:
                    matrix[i, j] = j * gap
                elif j == 0:
                    matrix[i, j] = i * gap
    # Fill the matrix within the band
    for i in range(1, n):
        for j in range(max(1, i - k), min(m, i + k + 1)):
            match_score = matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
            delete = matrix[i - 1, j] + gap
            insert = matrix[i, j - 1] + gap
            matrix[i, j] = max(match_score, delete, insert)
    # Traceback
    aligned_seq1 = []
    aligned_seq2 = []
    i = n - 1
    j = m - 1
    while i > 0 or j > 0:
        if i == 0:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        elif j == 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            score_current = matrix[i, j]
            score_diag = matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match, mismatch)
            if score_current == score_diag:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif score_current == matrix[i - 1, j] + gap:
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            else:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    aligned_seq1 = ''.join(aligned_seq1)
    aligned_seq2 = ''.join(aligned_seq2)
    return aligned_seq1, aligned_seq2, matrix[n-1, m-1]

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