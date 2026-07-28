"""
Unit tests for chance-corrected retention metrics and motor drive activation.
"""

import unittest
import numpy as np
from experiments.closed_loop_structural import (
    make_scaled_closed_loop_learner,
    run_trial,
    task_mapping,
    compute_chance_corrected_retention,
)


class TestRetentionMetrics(unittest.TestCase):
    def test_chance_corrected_retention_perfect(self):
        r = compute_chance_corrected_retention(0.90, 0.90, chance=0.50)
        self.assertAlmostEqual(r, 1.0, places=5)

    def test_chance_corrected_retention_chance_level(self):
        r = compute_chance_corrected_retention(0.50, 0.90, chance=0.50)
        self.assertAlmostEqual(r, 0.0, places=5)

    def test_chance_corrected_retention_negative(self):
        r = compute_chance_corrected_retention(0.30, 0.90, chance=0.50)
        self.assertAlmostEqual(r, -0.5, places=5)

    def test_chance_corrected_retention_unlearned_init(self):
        r = compute_chance_corrected_retention(0.50, 0.45, chance=0.50)
        self.assertEqual(r, 0.0)

    def test_task1_initial_acquisition_above_target(self):
        """Verify that fixed run_trial achieves Task 1 accuracy > 0.70."""
        scores = []
        for seed in range(3):
            engine, N_actual, beta_1 = make_scaled_closed_loop_learner(n_h=100, seed=seed)
            rng = np.random.default_rng(seed)
            rewards = []
            for t in range(200):
                x = int(rng.integers(2))
                target = task_mapping(0, x)
                r = run_trial(engine, x, target, rng, train=True)
                rewards.append(r)
            scores.append(np.mean(rewards[-50:]))
        
        mean_acc = np.mean(scores)
        self.assertGreater(mean_acc, 0.70, f"Task 1 accuracy {mean_acc:.3f} was <= 0.70 target")


if __name__ == "__main__":
    unittest.main()
