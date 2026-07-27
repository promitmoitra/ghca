"""
Unit tests for Phase 3 sequential learning anti-forgetting protocol.
"""

import unittest
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from experiments.sequential_closed_loop import evaluate_sequential_seed


class TestSequentialClosedLoop(unittest.TestCase):

    def test_evaluate_sequential_seed(self):
        """Test single seed sequential learning execution and dictionary return structure."""
        res = evaluate_sequential_seed(seed=0, axes=("tau", "theta", "w"))
        self.assertIn("acc_A1", res)
        self.assertIn("acc_B", res)
        self.assertIn("acc_A2", res)
        self.assertIn("retention_pct", res)
        
        self.assertGreaterEqual(res["acc_A1"], 0.0)
        self.assertLessEqual(res["acc_A1"], 1.0)
        self.assertGreaterEqual(res["acc_A2"], 0.0)
        self.assertLessEqual(res["acc_A2"], 1.0)
        self.assertGreaterEqual(res["retention_pct"], 0.0)


if __name__ == "__main__":
    unittest.main()
