"""
Learning Objectives
• Implement classical dynamic programming alignment algorithms (global & local).
• Implement affine gap penalties and banded/space-efficient variants.
• Implement heuristic pairwise alignment (BLAST-style).
• Implement multiple sequence alignment strategies (progressive, iterative refinement,		profile HMM).
• Benchmark algorithms on datasets with reproducible experiments.
• Write a structured report with figures, analysis, and recommendations

Pairwise Alignment (implement in R or Python):
1. Needleman–Wunsch (global)
2. Smith–Waterman (local)
3. Gotoh (affine gap penalties)
4. Banded/space-efficient alignment (Hirschberg/banded DP)
Heuristic/Approximate:
5. BLAST-style seed-and-extend
6. Greedy/minimizer-based approximate alignment
Multiple Sequence Alignment (MSA):
7. Progressive alignment (ClustalW-style)
8. Iterative refinement (MUSCLE-style)
9. Profile HMM alignment (Viterbi decoding)

Tasks
A. Implementations — Write modular, documented Python implementations with unit tests.
B. Experiments — Compare accuracy, runtime, memory, and parameter sensitivity.
C. Analysis — Complexity analysis, biological use cases, discussion of strengths/weaknesses.

"""

import numpy as np

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
    Global alignment using Needleman-Wunsch algorithm.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap
    :return: Aligned seq1, aligned seq2, score
    """
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m+1, n+1))
    # Initialize
    for i in range(1, m+1):
        dp[i][0] = dp[i-1][0] + gap
    for j in range(1, n+1):
        dp[0][j] = dp[0][j-1] + gap
    # Fill
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch)
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            dp[i][j] = max(diag, up, left)
    # Traceback
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2, dp[m][n]

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
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m+1, n+1))
    max_score = 0
    max_i, max_j = 0, 0
    # Fill
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch)
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            dp[i][j] = max(0, diag, up, left)
            if dp[i][j] > max_score:
                max_score = dp[i][j]
                max_i, max_j = i, j
    # Traceback from max
    align1, align2 = '', ''
    i, j = max_i, max_j
    while dp[i][j] > 0:
        if dp[i][j] == dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2, max_score

def gotoh(seq1, seq2, match=1, mismatch=-1, gap_open=-5, gap_extend=-1):
    """
    Global alignment with affine gap penalties using Gotoh algorithm.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap_open: Penalty for opening a gap
    :param gap_extend: Penalty for extending a gap
    :return: Aligned seq1, aligned seq2, score
    """
    m, n = len(seq1), len(seq2)
    M = np.full((m+1, n+1), -np.inf)
    Ix = np.full((m+1, n+1), -np.inf)
    Iy = np.full((m+1, n+1), -np.inf)
    
    # Initialize
    M[0][0] = 0
    Ix[0][0] = gap_open
    Iy[0][0] = gap_open
    for i in range(1, m+1):
        M[i][0] = -np.inf
        Ix[i][0] = gap_open + (i - 1) * gap_extend
        Iy[i][0] = gap_open + i * gap_extend
    for j in range(1, n+1):
        M[0][j] = -np.inf
        Ix[0][j] = gap_open + j * gap_extend
        Iy[0][j] = gap_open + (j - 1) * gap_extend
    
    # Fill
    for i in range(1, m+1):
        for j in range(1, n+1):
            s = score(seq1[i-1], seq2[j-1], match, mismatch)
            M[i][j] = s + max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1])
            Ix[i][j] = max(gap_open + M[i][j-1], gap_extend + Ix[i][j-1])
            Iy[i][j] = max(gap_open + M[i-1][j], gap_extend + Iy[i-1][j])
    
    # Find max score
    final_score = max(M[m][n], Ix[m][n], Iy[m][n])
    
    # Traceback
    align1, align2 = '', ''
    i, j = m, n
    current_matrix = 'M' if final_score == M[m][n] else ('Ix' if final_score == Ix[m][n] else 'Iy')
    
    while i > 0 or j > 0:
        if current_matrix == 'M':
            if i > 0 and j > 0 and M[i][j] == score(seq1[i-1], seq2[j-1], match, mismatch) + M[i-1][j-1]:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
                current_matrix = 'M'
            elif i > 0 and j > 0 and M[i][j] == score(seq1[i-1], seq2[j-1], match, mismatch) + Ix[i-1][j-1]:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
                current_matrix = 'Ix'
            elif i > 0 and j > 0 and M[i][j] == score(seq1[i-1], seq2[j-1], match, mismatch) + Iy[i-1][j-1]:
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
                current_matrix = 'Iy'
        elif current_matrix == 'Ix':
            if Ix[i][j] == gap_open + M[i][j-1]:
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1
                current_matrix = 'M'
            else:
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1
                current_matrix = 'Ix'
        elif current_matrix == 'Iy':
            if Iy[i][j] == gap_open + M[i-1][j]:
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
                current_matrix = 'M'
            else:
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
                current_matrix = 'Iy'
        if i == 0 and j == 0:
            break
    return align1, align2, final_score

def banded_dp(seq1, seq2, k, match=1, mismatch=-1, gap=-2):
    """
    Banded dynamic programming for approximate alignment.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param k: Band width
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap
    :return: Aligned seq1, aligned seq2, score
    """
    m, n = len(seq1), len(seq2)
    dp = np.full((m+1, n+1), -np.inf)
    # Initialize within band
    dp[0][0] = 0
    for i in range(1, min(m+1, k+1)):
        dp[i][0] = dp[i-1][0] + gap
    for j in range(1, min(n+1, k+1)):
        dp[0][j] = dp[0][j-1] + gap
    # Fill within band
    for i in range(1, m+1):
        for j in range(max(1, i-k), min(n+1, i+k+1)):
            diag = dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch)
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap
            dp[i][j] = max(diag, up, left)
    # Traceback
    align1, align2 = '', ''
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2, dp[m][n]

def blast_seed_extend(seq1, seq2, k=11, match=1, mismatch=-1, gap=-2):
    """
    BLAST-style seed-and-extend heuristic alignment.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param k: Seed length
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap
    :return: Aligned seq1, aligned seq2, score
    """
    seeds = []
    for i in range(len(seq1) - k + 1):
        substr = seq1[i:i+k]
        pos = seq2.find(substr)
        while pos != -1:
            seeds.append((i, pos))
            pos = seq2.find(substr, pos + 1)
    best_align = ('', '', 0)
    for start1, start2 in seeds:
        # Extend left
        left1, left2 = '', ''
        i, j = start1 - 1, start2 - 1
        while i >= 0 and j >= 0 and seq1[i] == seq2[j]:
            left1 = seq1[i] + left1
            left2 = seq2[j] + left2
            i -= 1
            j -= 1
        # Extend right
        right1, right2 = '', ''
        i, j = start1 + k, start2 + k
        while i < len(seq1) and j < len(seq2) and seq1[i] == seq2[j]:
            right1 += seq1[i]
            right2 += seq2[j]
            i += 1
            j += 1
        # Combine
        align1 = left1 + seq1[start1:start1+k] + right1
        align2 = left2 + seq2[start2:start2+k] + right2
        score = len(align1) * match
        if score > best_align[2]:
            best_align = (align1, align2, score)
    return best_align
