"""
Unit tests for Phase 2 closed-loop plasticity benchmarks.
"""

import unittest
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiments.closed_loop_plasticity import (
    evaluate_e1_seed,
    evaluate_e5_seed,
    fit_linear_decoder
)


class TestClosedLoopBenchmarks(unittest.TestCase):

    def test_fit_linear_decoder(self):
        """Test linear decoder fitting on synthetic data."""
        rng = np.random.default_rng(42)
        X_train = rng.standard_normal((100, 10))
        y_train = (X_train[:, 0] > 0).astype(int)
        
        X_test = rng.standard_normal((50, 10))
        y_test = (X_test[:, 0] > 0).astype(int)

        acc = fit_linear_decoder(X_train, y_train, X_test, y_test, seed=42)
        self.assertGreater(acc, 0.70)

    def test_e1_seed_execution(self):
        """Test E1 single-seed evaluation returns valid accuracy and RIR bounds."""
        f_acc, t_acc, rir, rewards = evaluate_e1_seed(seed=0, axes=("tau", "theta", "w"))
        self.assertGreaterEqual(f_acc, 0.0)
        self.assertLessEqual(f_acc, 1.0)
        self.assertGreaterEqual(t_acc, 0.0)
        self.assertLessEqual(t_acc, 1.0)
        self.assertGreaterEqual(rir, 0.0)
        self.assertEqual(len(rewards), 300)

    def test_e5_seed_execution(self):
        """Test E5 executive switch single-seed evaluation returns valid bounds."""
        f_acc, t_acc, rir, rewards = evaluate_e5_seed(seed=0, axes=("tau", "theta", "w"))
        self.assertGreaterEqual(f_acc, 0.0)
        self.assertLessEqual(f_acc, 1.0)
        self.assertGreaterEqual(t_acc, 0.0)
        self.assertLessEqual(t_acc, 1.0)
        self.assertGreaterEqual(rir, 0.0)
        self.assertEqual(len(rewards), 300)


if __name__ == "__main__":
    unittest.main()
