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

    def test_saved_npz_archive_structure_and_swap_size_fields(self):
        """Verify that the committed NPZ archive contains the required swap-size distribution fields for P-C."""
        npz_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "result", "topology", "coherence_window_4x4.npz")
        self.assertTrue(os.path.exists(npz_path), f"Missing artifact {npz_path}")
        d = np.load(npz_path)

        # Check required keys for P-C comparison
        required_keys = [
            "open_ceil_sizes_keys", "open_ceil_sizes_vals",
            "torus_ceil_sizes_keys", "torus_ceil_sizes_vals",
            "open_h1_sizes_keys", "open_h1_sizes_vals",
            "torus_h1_sizes_keys", "torus_h1_sizes_vals"
        ]
        for rk in required_keys:
            self.assertIn(rk, d.files, f"Missing key {rk} in {npz_path}")

        # Compute swap-size statistics from saved archive
        open_keys, open_vals = d["open_ceil_sizes_keys"], d["open_ceil_sizes_vals"]
        torus_keys, torus_vals = d["torus_ceil_sizes_keys"], d["torus_ceil_sizes_vals"]

        open_mean = np.sum(open_keys * open_vals) / np.sum(open_vals)
        torus_mean = np.sum(torus_keys * torus_vals) / np.sum(torus_vals)

        p_single_open = open_vals[0] / np.sum(open_vals)  # key 1 is index 0
        p_single_torus = torus_vals[0] / np.sum(torus_vals)  # key 1 is index 0

        # Verify exact reproducibility from archive
        self.assertEqual(int(d["open_ceil_u"].item()), 115)
        self.assertEqual(int(d["open_ceil_total"].item()), 115)
        self.assertEqual(int(d["open_ceil_v"].item()), 0)
        self.assertEqual(int(d["open_h1_v"].item()), 246)
        self.assertEqual(int(d["open_h1_total"].item()), 246)
        self.assertEqual(int(d["open_h1_u"].item()), 0)

        self.assertEqual(int(d["torus_ceil_u"].item()), 22)
        self.assertEqual(int(d["torus_ceil_total"].item()), 22)
        self.assertEqual(int(d["torus_ceil_v"].item()), 0)

        self.assertEqual(int(open_keys[-1]), 9, "Uncapped search must reveal max swap size 9 tail")
        self.assertAlmostEqual(open_mean, 2.2696, places=3)
        self.assertAlmostEqual(torus_mean, 1.7727, places=3)
        self.assertAlmostEqual(p_single_open, 0.3739, places=3)
        self.assertAlmostEqual(p_single_torus, 0.5000, places=3)


if __name__ == "__main__":
    unittest.main()
