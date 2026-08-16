"""
Unit tests for the Coherence Invariant on 4x4 open and torus lattices.
Tests predictions P-A (window bounded by S), P-B (boundary role concentration),
and witness side-separation.
"""
import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiments.coherence_window_4x4 import run_sampling, ROLE4


class TestCoherenceWindow4x4(unittest.TestCase):

    def test_window_bound_and_witness_separation_open(self):
        """Verify that on 4x4 open lattice, max BFS depth <= S, and witness sides do not mix."""
        res = run_sampling(ta=2, tp=1, periodic=False, n_trials=50, seed=42)
        max_d = max(k for k in res['depth_hist'].keys() if k is not None)
        self.assertLessEqual(max_d, res['S'], f"Max depth {max_d} exceeded S={res['S']}")
        self.assertEqual(res['ceil_v_witnessed'], 0, "Ceiling holds must never be witnessed on v-side")
        self.assertEqual(res['h1_u_witnessed'], 0, "Age-1 holds must never be witnessed on u-side")
        if res['ceil_roles']['corner'] > 0 and res['ceil_roles']['edge'] > 0 and res['ceil_roles']['centre'] > 0:
            r_corner = res['ceil_roles']['corner'] / 4.0
            r_edge = res['ceil_roles']['edge'] / 8.0
            r_centre = res['ceil_roles']['centre'] / 4.0
            self.assertGreaterEqual(r_corner, r_edge, "Corner per-capita rate should be >= edge rate")

    def test_window_bound_and_witness_separation_torus(self):
        """Verify that on 4x4 periodic torus, max BFS depth <= S, and witness sides do not mix."""
        res = run_sampling(ta=2, tp=1, periodic=True, n_trials=50, seed=42)
        max_d = max(k for k in res['depth_hist'].keys() if k is not None)
        self.assertLessEqual(max_d, res['S'], f"Max depth {max_d} exceeded S={res['S']}")
        self.assertEqual(res['ceil_v_witnessed'], 0, "Ceiling holds must never be witnessed on v-side")
        self.assertEqual(res['h1_u_witnessed'], 0, "Age-1 holds must never be witnessed on u-side")


if __name__ == "__main__":
    unittest.main()
