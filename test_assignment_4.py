import unittest
from Assignment_4 import needleman_wunsch, smith_waterman, gotoh, score, banded_dp, blast_seed_extend, greedy, progressive_alignment, iterative_refinement

class TestSequenceAlignment(unittest.TestCase):

    def test_score_function(self):
        # Test match
        self.assertEqual(score('A', 'A'), 1)
        # Test mismatch
        self.assertEqual(score('A', 'T'), -1)
        # Test with custom scores
        self.assertEqual(score('A', 'A', match=2, mismatch=-3), 2)
        self.assertEqual(score('A', 'T', match=2, mismatch=-3), -3)

    def test_needleman_wunsch_perfect_match(self):
        seq1 = "ATCG"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = needleman_wunsch(seq1, seq2)
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 4)  # 4 matches

    def test_needleman_wunsch_with_gaps(self):
        seq1 = "ATCG"
        seq2 = "ACG"
        aligned1, aligned2, score_val = needleman_wunsch(seq1, seq2)
        # Expected: ATCG and A-CG, score = 1 + (-2) + 1 + 1 = 1
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "A-CG")
        self.assertEqual(score_val, 1)

    def test_needleman_wunsch_different_lengths(self):
        seq1 = "AT"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = needleman_wunsch(seq1, seq2)
        self.assertEqual(aligned1, "AT--")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 2 + 2*(-2))  # 2 matches, 2 gaps

    def test_smith_waterman_perfect_match(self):
        seq1 = "ATCG"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = smith_waterman(seq1, seq2)
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 4)

    def test_smith_waterman_local(self):
        seq1 = "ATCGT"
        seq2 = "TCG"
        aligned1, aligned2, score_val = smith_waterman(seq1, seq2)
        # Should align TCG with TCG
        self.assertEqual(aligned1, "TCG")
        self.assertEqual(aligned2, "TCG")
        self.assertEqual(score_val, 3)

    def test_smith_waterman_no_match(self):
        seq1 = "AAA"
        seq2 = "TTT"
        aligned1, aligned2, score_val = smith_waterman(seq1, seq2)
        # No positive score, should return empty or minimal
        self.assertEqual(score_val, 0)
        self.assertEqual(aligned1, "")
        self.assertEqual(aligned2, "")

    def test_gotoh_perfect_match(self):
        seq1 = "ATCG"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = gotoh(seq1, seq2)
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 4)

    def test_gotoh_with_affine_gaps(self):
        seq1 = "ATCG"
        seq2 = "ACG"
        aligned1, aligned2, score_val = gotoh(seq1, seq2, gap_open=-5, gap_extend=-1)
        # With affine, gap of 1 should be -5, score = 1 + 1 + 1 -5 = -2
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "A-CG")
        self.assertEqual(score_val, -2)

    def test_gotoh_different_lengths(self):
        seq1 = "AT"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = gotoh(seq1, seq2, gap_open=-5, gap_extend=-1)
        # Gaps at end: two gaps, but since affine, opening once, extending once: -5 -1 = -6
        self.assertEqual(aligned1, "AT--")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 2 - 6)  # 2 matches, gap penalty

    def test_banded_dp_perfect_match(self):
        seq1 = "ATCG"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = banded_dp(seq1, seq2, k=1)
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 4)

    def test_banded_dp_with_gaps(self):
        seq1 = "ATCG"
        seq2 = "ACG"
        aligned1, aligned2, score_val = banded_dp(seq1, seq2, k=1)
        # Expected: ATCG and A-CG, score = 1 + (-2) + 1 + 1 = 1
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "A-CG")
        self.assertEqual(score_val, 1)

    def test_banded_dp_different_lengths(self):
        seq1 = "AT"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = banded_dp(seq1, seq2, k=2)
        self.assertEqual(aligned1, "AT--")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 2 + 2*(-2))  # 2 matches, 2 gaps

    def test_blast_seed_extend_perfect_match(self):
        seq1 = "ATCGATCGATCG"
        seq2 = "ATCGATCGATCG"
        aligned1, aligned2, score_val = blast_seed_extend(seq1, seq2, k=4)
        self.assertEqual(aligned1, "ATCGATCGATCG")
        self.assertEqual(aligned2, "ATCGATCGATCG")
        self.assertEqual(score_val, 12)

    def test_blast_seed_extend_partial(self):
        seq1 = "ATCGTAAAA"
        seq2 = "ATCGAAAAA"
        aligned1, aligned2, score_val = blast_seed_extend(seq1, seq2, k=4)
        # Should find ATCG, then stops at mismatch
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 4)

    def test_blast_seed_extend_no_match(self):
        seq1 = "AAAA"
        seq2 = "TTTT"
        aligned1, aligned2, score_val = blast_seed_extend(seq1, seq2, k=2)
        self.assertEqual(score_val, 0)
        self.assertEqual(aligned1, "")
        self.assertEqual(aligned2, "")

    def test_greedy_perfect_match(self):
        seq1 = "ATCG"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = greedy(seq1, seq2, 1, -1, -2)
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 4)

    def test_greedy_with_mismatch(self):
        seq1 = "ATCG"
        seq2 = "ACCG"
        aligned1, aligned2, score_val = greedy(seq1, seq2, 1, -1, -2)
        # Simple greedy: no looka head, no gaps unless sequence ends
        # A match, T vs C mismatch, then C match, G match
        self.assertEqual(aligned1, "ATCG")
        self.assertEqual(aligned2, "ACCG")
        self.assertEqual(score_val, 1 - 1 + 1 + 1) # = 2

    def test_greedy_different_lengths(self):
        seq1 = "AT"
        seq2 = "ATCG"
        aligned1, aligned2, score_val = greedy(seq1, seq2, 1, -1, -2)
        # A matches, T matches, then gap for CG
        self.assertEqual(aligned1, "AT--")
        self.assertEqual(aligned2, "ATCG")
        self.assertEqual(score_val, 1 + 1 + (-2) + (-2))  # 2 matches, 2 gaps

    def test_progressive_alignment(self):
        """Test progressive alignment with diverse sequences"""
        seqs = ['ATCGATCG', 'AGCGATCG', 'ATCGGTCG', 'AGCGGTCG']
        result = progressive_alignment(seqs)
        aligned_seqs, score_val = result
        self.assertEqual(len(aligned_seqs), 4)
        self.assertGreater(len(aligned_seqs[0]), 0)
        self.assertGreater(score_val, 0)

    def test_iterative_refinement(self):
        """Test iterative refinement with diverse sequences"""
        seqs = ['ATCGATCG', 'AGCGATCG', 'ATCGGTCG', 'AGCGGTCG']
        result = iterative_refinement(seqs, max_iterations=5)
        aligned_seqs, score_val = result
        self.assertEqual(len(aligned_seqs), 4)
        self.assertGreater(len(aligned_seqs[0]), 0)
        self.assertGreater(score_val, 0)


if __name__ == '__main__':
    unittest.main()
