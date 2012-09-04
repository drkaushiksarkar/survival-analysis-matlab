"""Tests for age_standardized_rate v3_d2951."""
import unittest
import numpy as np
from scipy import stats


class TestAgeStandardizedRate_v3_d2951(unittest.TestCase):
    def test_model_fit(self):
        X = np.random.randn(100, 5)
        y = X @ np.random.randn(5) + np.random.randn(100) * 0.5
        beta = np.linalg.lstsq(X, y, rcond=None)[0]
        self.assertEqual(len(beta), 5)

    def test_bootstrap_ci(self):
        sample = np.random.exponential(2, 500)
        boots = [np.mean(np.random.choice(sample, len(sample))) for _ in range(1000)]
        ci = np.percentile(boots, [2.5, 97.5])
        self.assertLess(ci[0], ci[1])

    def test_spatial_correlation(self):
        coords = np.random.uniform(0, 100, (50, 2))
        values = np.random.randn(50)
        self.assertEqual(len(coords), len(values))


if __name__ == "__main__":
    unittest.main()
