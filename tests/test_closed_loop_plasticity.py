"""
Unit tests for ClosedLoopPlasticityEngine and multi-axis substrate plasticity rules.
"""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import (
    ClosedLoopPlasticityEngine,
    compute_readout_independence_ratio,
    make_closed_loop_learner
)


class TestClosedLoopPlasticity(unittest.TestCase):

    def test_seed_reproducibility(self):
        """Verify that two learners initialized with identical seeds yield bit-identical states."""
        seed = 42
        learner1 = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)
        learner2 = make_closed_loop_learner(axes=("tau", "theta", "w"), seed=seed)

        rng1 = np.random.default_rng(seed + 10)
        rng2 = np.random.default_rng(seed + 10)

        for _ in range(50):
            drive = rng1.random(learner1.N) < 0.05
            _ = rng2.random(learner2.N) < 0.05  # keep RNG aligned

            learner1.step_learn(drive)
            learner2.step_learn(drive)

            delta1 = rng1.standard_normal()
            delta2 = rng2.standard_normal()

            learner1.learn(delta1)
            learner2.learn(delta2)

        np.testing.assert_array_equal(learner1.phi, learner2.phi)
        np.testing.assert_array_almost_equal(learner1.tau, learner2.tau)
        np.testing.assert_array_almost_equal(learner1.theta, learner2.theta)
        np.testing.assert_array_almost_equal(learner1.W, learner2.W)

    def test_axes_isolation(self):
        """Verify that enabling single axes modifies only the corresponding parameters."""
        seed = 123
        
        # Tau only
        learner_tau = make_closed_loop_learner(axes=("tau",), seed=seed)
        w_orig = learner_tau.W.copy()
        theta_orig = learner_tau.theta.copy()
        
        learner_tau.step_learn(drive=np.ones(learner_tau.N, dtype=bool))
        learner_tau.learn(delta=1.0)
        
        np.testing.assert_array_equal(learner_tau.W, w_orig)
        np.testing.assert_array_equal(learner_tau.theta, theta_orig)

        # W only
        learner_w = make_closed_loop_learner(axes=("w",), seed=seed)
        tau_orig = learner_w.tau.copy()
        theta_orig_w = learner_w.theta.copy()
        
        learner_w.step_learn(drive=np.ones(learner_w.N, dtype=bool))
        learner_w.learn(delta=1.0)
        
        np.testing.assert_array_equal(learner_w.tau, tau_orig)
        np.testing.assert_array_equal(learner_w.theta, theta_orig_w)

    def test_parameter_bounds(self):
        """Verify that parameter updates remain strictly within specified min/max bounds."""
        learner = make_closed_loop_learner(
            axes=("tau", "theta", "w"),
            seed=999
        )

        # Force extreme positive and negative updates
        for delta in [100.0, -100.0, 500.0, -500.0]:
            learner.step_learn(drive=np.ones(learner.N, dtype=bool))
            learner.learn(delta)

            assert np.all(learner.tau >= learner.tau_min)
            assert np.all(learner.tau <= learner.tau_max)
            assert np.all(learner.theta >= learner.theta_bounds[0])
            assert np.all(learner.theta <= learner.theta_bounds[1])
            assert np.all(learner.W >= 0.0)
            assert np.all(learner.W <= learner.w_max)

    def test_readout_independence_ratio(self):
        """Verify Readout Independence Ratio (RIR) calculation logic."""
        rir_perfect = compute_readout_independence_ratio(0.90, 0.90)
        self.assertAlmostEqual(rir_perfect, 1.0, places=4)

        rir_partial = compute_readout_independence_ratio(0.75, 0.90)
        self.assertAlmostEqual(rir_partial, 0.8333, places=3)

        rir_zero = compute_readout_independence_ratio(0.0, 0.90)
        self.assertAlmostEqual(rir_zero, 0.0, places=4)


if __name__ == "__main__":
    unittest.main()
