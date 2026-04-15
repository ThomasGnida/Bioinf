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
    Banded dynamic programming for global alignment,
    restricting to band around diagonal.
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
    """
    BLAST seed and extend algorithm for local alignment.
    Finds seeds of length k, then extends ungapped.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param k: Seed length
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap (not used in this simple version)
    :return: Aligned seq1, aligned seq2, score
    """
    n = len(seq1)
    m = len(seq2)
    best_score = 0
    best_align1 = ""
    best_align2 = ""
    # Find seeds
    for i in range(n - k + 1):
        for j in range(m - k + 1):
            if seq1[i:i+k] == seq2[j:j+k]:
                # Extend left
                left_i = i
                left_j = j
                while left_i > 0 and left_j > 0 and seq1[left_i-1] == seq2[left_j-1]:
                    left_i -= 1
                    left_j -= 1
                # Extend right
                right_i = i + k
                right_j = j + k
                while right_i < n and right_j < m and seq1[right_i] == seq2[right_j]:
                    right_i += 1
                    right_j += 1
                # Calculate score
                align_len = right_i - left_i
                score = align_len * match
                if score > best_score:
                    best_score = score
                    best_align1 = seq1[left_i:right_i]
                    best_align2 = seq2[left_j:right_j]
    if best_score == 0:
        return "", "", 0
    return best_align1, best_align2, best_score

def greedy(seq1, seq2, match=1, mismatch=-1, gap=-2):
    """
    Greedy alignment algorithm. Processes sequences left to right,
    greedily matching characters when possible.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap
    :return: Aligned seq1, aligned seq2, score
    """
    n = len(seq1)
    m = len(seq2)
    aligned_seq1 = []
    aligned_seq2 = []
    score_val = 0
    i = 0
    j = 0
    # Greedily process both sequences
    while i < n and j < m:
        if seq1[i] == seq2[j]:
            # Match found, align them
            aligned_seq1.append(seq1[i])
            aligned_seq2.append(seq2[j])
            score_val += match
            i += 1
            j += 1
        else:
            # Try to find a match ahead
            # Look for the next match in seq2 that appears in remaining seq1
            found = False
            for ii in range(i, min(i + 5, n)):  # Look ahead up to 5 chars
                for jj in range(j, min(j + 5, m)):  # Look ahead up to 5 chars
                    if seq1[ii] == seq2[jj]:
                        # Found a match ahead, add gaps until we reach it
                        for _ in range(ii - i):
                            aligned_seq1.append(seq1[i])
                            aligned_seq2.append('-')
                            score_val += gap
                            i += 1
                        for _ in range(jj - j):
                            aligned_seq1.append('-')
                            aligned_seq2.append(seq2[j])
                            score_val += gap
                            j += 1
                        found = True
                        break
                if found:
                    break
            if not found:
                # No match found ahead, add mismatch and move on
                aligned_seq1.append(seq1[i])
                aligned_seq2.append(seq2[j])
                score_val += mismatch
                i += 1
                j += 1
    # Add remaining characters as gaps
    while i < n:
        aligned_seq1.append(seq1[i])
        aligned_seq2.append('-')
        score_val += gap
        i += 1
    while j < m:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[j])
        score_val += gap
        j += 1
    aligned_seq1_str = ''.join(aligned_seq1)
    aligned_seq2_str = ''.join(aligned_seq2)
    return aligned_seq1_str, aligned_seq2_str, score_val

def progressive_alignment(sequences, match=1, mismatch=-1, gap=-2):
    """
    Progressive alignment for multiple sequences using a guide tree approach.
    Builds pairwise alignments and combines them based on sequence similarity.
    :param sequences: List of sequences to align
    :param match: Score for match
    :param mismatch: Penalty for mismatch
    :param gap: Penalty for gap
    :return: List of aligned sequences, total alignment score
    """
    if not sequences:
        return [], 0
    if len(sequences) == 1:
        return sequences, 0
    
    # Build pairwise distance matrix to find most similar sequences
    n_seqs = len(sequences)
    distances = {}
    
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            # Use Needleman-Wunsch to get alignment score, convert to distance
            _, _, align_score = needleman_wunsch(sequences[i], sequences[j], match, mismatch, gap)
            distances[(i, j)] = -align_score
            distances[(j, i)] = -align_score
    
    # Build guide tree by repeatedly combining most similar sequences
    active_seqs = list(range(n_seqs))
    aligned_dict = {i: sequences[i] for i in range(n_seqs)}
    next_idx = n_seqs
    
    while len(active_seqs) > 1:
        # Find closest pair
        min_dist = np.inf
        merge_i, merge_j = 0, 1
        for ii in range(len(active_seqs)):
            for jj in range(ii + 1, len(active_seqs)):
                idx_i = active_seqs[ii]
                idx_j = active_seqs[jj]
                dist = distances.get((idx_i, idx_j), distances.get((idx_j, idx_i), np.inf))
                if dist < min_dist:
                    min_dist = dist
                    merge_i = ii
                    merge_j = jj
        
        # Align the two closest sequences/profiles
        idx_i = active_seqs[merge_i]
        idx_j = active_seqs[merge_j]
        seq_i = aligned_dict[idx_i]
        seq_j = aligned_dict[idx_j]
        
        align_i, align_j, _ = needleman_wunsch(seq_i, seq_j, match, mismatch, gap)
        
        # Store the merged alignment under a new index
        aligned_dict[next_idx] = align_i
        aligned_dict[idx_i] = align_i
        aligned_dict[idx_j] = align_j
        
        # Update distances for remaining sequences
        for k in active_seqs:
            if k != idx_i and k != idx_j:
                dist_i_k = distances.get((idx_i, k), distances.get((k, idx_i), 0))
                dist_j_k = distances.get((idx_j, k), distances.get((k, idx_j), 0))
                new_dist = (dist_i_k + dist_j_k) / 2.0
                distances[(next_idx, k)] = new_dist
                distances[(k, next_idx)] = new_dist
        
        # Remove merged sequences and add new node
        active_seqs.remove(idx_i)
        active_seqs.remove(idx_j)
        active_seqs.append(next_idx)
        next_idx += 1
    
    # Reconstruct final alignment from original indices
    result_aligned = [aligned_dict[idx] for idx in range(n_seqs)]
    final_score = sum(needleman_wunsch(result_aligned[i], result_aligned[j], 
                                       match, mismatch, gap)[2] 
                      for i in range(len(result_aligned))
                      for j in range(i + 1, len(result_aligned)))
    
    return result_aligned, final_score

def iterative_refinement(sequences, match=1, mismatch=-1, gap=-2, max_iterations=10):
    """
    Iterative refinement of multiple sequence alignment.
    Removes one sequence at a time and realigns it, keeping improvement if score increases.
    :param sequences: List of sequences to align
    :param match: Score for match
    :param mismatch: Penalty for mismatch
    :param gap: Penalty for gap
    :param max_iterations: Maximum number of refinement iterations
    :return: List of refined aligned sequences, final alignment score
    """
    if len(sequences) <= 1:
        return sequences, 0
    
    # Start with progressive alignment
    current_align, current_score = progressive_alignment(sequences, match, mismatch, gap)
    
    for iteration in range(max_iterations):
        improved = False
        
        # Try removing and realigning each sequence
        for idx in range(len(current_align)):
            # Create alignment without this sequence
            other_seqs = [current_align[i] for i in range(len(current_align)) if i != idx]
            removed_seq = current_align[idx]
            
            # Align the removed sequence to the alignment of others
            # Simple approach: align to first other sequence
            if other_seqs:
                align_removed, align_other, new_pair_score = needleman_wunsch(
                    removed_seq, other_seqs[0], match, mismatch, gap
                )
                
                # Recalculate total score
                new_align = other_seqs[:]
                new_align[0] = align_other
                new_align.insert(idx, align_removed)
                
                # Calculate score of new alignment
                new_score = 0
                for i in range(len(new_align)):
                    for j in range(i + 1, len(new_align)):
                        # Count matches in aligned columns
                        col_score = 0
                        min_len = min(len(new_align[i]), len(new_align[j]))
                        for col in range(min_len):
                            if new_align[i][col] != '-' and new_align[j][col] != '-':
                                col_score += score(new_align[i][col], new_align[j][col], match, mismatch)
                            elif new_align[i][col] == '-' or new_align[j][col] == '-':
                                col_score += gap
                        new_score += col_score
                
                # Accept if improved
                if new_score > current_score:
                    current_align = new_align
                    current_score = new_score
                    improved = True
        
        # If no improvement found, stop
        if not improved:
            break
    
    return current_align, current_score

def profile_hmm_alignment():
    return