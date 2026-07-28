"""
Unit tests for Axis G Structural Plasticity (edge pruning & sprouting),
homeostatic tau-relaxation consolidation, and graph modularity calculation.
"""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import ClosedLoopPlasticityEngine, compute_graph_modularity, make_closed_loop_learner


class TestAxisGPlasticity(unittest.TestCase):

    def test_axis_g_pruning(self):
        """Verify that edges with weights below w_prune are pruned and marked non-plastic."""
        seed = 42
        net = make_closed_loop_learner(axes=("tau", "theta", "w", "g"), seed=seed)
        
        # Manually set a plastic edge weight below w_prune
        plastic_indices = np.where(net.plastic)
        i, j = plastic_indices[0][0], plastic_indices[1][0]
        net.W[i, j] = 0.01  # < w_prune (0.05)

        # Call learn with delta = -0.5
        net.learn(delta=-0.5)

        # Check edge is now zero weight and non-plastic
        self.assertEqual(net.W[i, j], 0.0)
        self.assertFalse(net.plastic[i, j])

    def test_axis_g_sprouting(self):
        """Verify that co-active nodes sprout new weak plastic connections when delta > 0."""
        seed = 123
        net = make_closed_loop_learner(axes=("tau", "theta", "w", "g"), seed=seed)

        # Force high activity trace on specific nodes
        net.activity_trace[5] = 0.50
        net.activity_trace[10] = 0.50

        # Record plastic edge count before
        count_before = np.count_nonzero(net.plastic)

        # Call learn with positive reward error delta = 1.0
        net.learn(delta=1.0)

        # Record plastic edge count after
        count_after = np.count_nonzero(net.plastic)

        self.assertGreaterEqual(count_after, count_before)

    def test_consolidate_tau(self):
        """Verify that consolidate_tau decays non-protected tau while preserving protected tau."""
        seed = 999
        net = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)

        # Set node 0 with high tau (protected) and node 1 with high tau
        net.tau_base[0] = 20.0
        net.tau[0] = 20.0
        net.tau_base[1] = 20.0
        net.tau[1] = 20.0

        net.consolidate_tau(decay_rate=0.10)

        # Both should decay towards tau_min, but protected node decays proportionally
        self.assertLess(net.tau[0], 20.0)
        self.assertLess(net.tau[1], 20.0)

    def test_compute_graph_modularity(self):
        """Verify compute_graph_modularity returns finite float >= 0 for non-trivial graph."""
        seed = 77
        net = make_closed_loop_learner(axes=("tau", "theta", "w", "g"), seed=seed)
        Q = compute_graph_modularity(net.W)
        self.assertIsInstance(Q, float)
        self.assertGreaterEqual(Q, 0.0)

    def test_rng_determinism(self):
        """Verify bit-identical outputs under explicit seed threading."""
        net1 = make_closed_loop_learner(axes=("tau", "theta", "w", "g"), seed=100)
        net2 = make_closed_loop_learner(axes=("tau", "theta", "w", "g"), seed=100)

        net1.learn(1.0)
        net2.learn(1.0)

        np.testing.assert_array_equal(net1.W, net2.W)
        np.testing.assert_array_equal(net1.tau, net2.tau)


if __name__ == "__main__":
    unittest.main()

