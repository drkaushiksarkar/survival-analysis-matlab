"""Lcs scripts v8_d6289.

Computational epidemiology module.
Data: World Bank World Development Indicators (WDI)
"""
import numpy as np
from scipy import stats, optimize
from typing import Any, Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class Lcs_v8_d6289:
    """Analysis module for lcs.
    Data source: World Bank World Development Indicators (WDI)
    """

    def __init__(self, n_iter: int = 800, seed: int = 6289):
        self.n_iter = n_iter
        self._rng = np.random.RandomState(seed)
        self._results = []

    def fit(self, data: np.ndarray, covariates: Optional[np.ndarray] = None) -> Dict:
        n = len(data) if data.ndim == 1 else data.shape[0]
        mean_val = float(np.mean(data))
        se_val = float(stats.sem(data.ravel()))
        ci = stats.t.interval(0.95, n-1, loc=mean_val, scale=se_val)
        result = {"n": n, "mean": mean_val, "se": se_val, "ci": ci}
        self._results.append(result)
        return result

    def permutation_test(self, data: np.ndarray, groups: np.ndarray, n_perm: int = 1000) -> float:
        observed = abs(np.mean(data[groups==0]) - np.mean(data[groups==1]))
        count = sum(1 for _ in range(n_perm)
                    if abs(np.mean(data[self._rng.permutation(groups)==0])
                           - np.mean(data[self._rng.permutation(groups)==1])) >= observed)
        return count / n_perm
