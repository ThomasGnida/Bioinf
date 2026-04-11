"""
Assignment in Sequence Alignment

Title: A comprehensive implementation, evaluation and comparison of sequence alignment methods
Submission: Single PDF report + code repository (Git).

Learning Objectives
• Implement classical dynamic programming alignment algorithms (global & local).
• Implement affine gap penalties and banded/space-efficient variants.
• Implement heuristic pairwise alignment (BLAST-style).
• Implement multiple sequence alignment strategies (progressive, iterative refinement,		profile HMM).
• Benchmark algorithms on datasets with reproducible experiments.
• Write a structured report with figures, analysis, and recommendations.

Required Methods

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
import Bio
import pandas as pd


def benchmark_algorithm(sequence_array, algorithm_method):
    accuracy = 0
    runtime = 0
    memory = 0
    parameter_sensitivity = 0
    complexity = 0


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
    Global alignment using Needleman-Wunsch algorithm.
    :param seq1: First sequence
    :param seq2: Second sequence
    :param match: Score for match
    :param mismatch: Score for mismatch
    :param gap: Score for gap
    :return: Aligned seq1, aligned seq2, score
    """
    return

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

    return

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

    return
def banded_dp(seq1, seq2, k, match=1, mismatch=-1, gap=-2):
    """

    """

    return

def blast_seed_extend(seq1, seq2, k=11, match=1, mismatch=-1, gap=-2):
    """

    """
    return

def greedy():
    return

def progressive_alignment():
    return

def iterative_refinement():
    return

def profile_hmm_alignment():
    return