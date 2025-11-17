"""
Utility functions for NE ratio calculations with bootstrap support.
"""

import re
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class BootstrapCalculator:
    """Calculate bootstrap statistics for Nₑ ratios."""
    
    @staticmethod
    def bootstrap_ci(data: np.ndarray, statistic_func: callable,
                    n_bootstraps: int = 1000, confidence: float = 0.95,
                    random_state: Optional[int] = None) -> Tuple[float, float, float]:
        """
        Calculate bootstrap confidence interval.
        
        Args:
            data: Input data array
            statistic_func: Function to calculate statistic
            n_bootstraps: Number of bootstrap samples
            confidence: Confidence level
            random_state: Random seed
            
        Returns:
            Tuple of (mean, lower_ci, upper_ci)
        """
        if random_state is not None:
            np.random.seed(random_state)
            
        n = len(data)
        bootstrap_stats = []
        
        for _ in range(n_bootstraps):
            # Sample with replacement
            bootstrap_sample = np.random.choice(data, size=n, replace=True)
            stat = statistic_func(bootstrap_sample)
            bootstrap_stats.append(stat)
        
        bootstrap_stats = np.array(bootstrap_stats)
        mean_stat = np.mean(bootstrap_stats)
        
        # Percentile method for CI
        alpha = (1 - confidence) / 2
        lower = np.percentile(bootstrap_stats, alpha * 100)
        upper = np.percentile(bootstrap_stats, (1 - alpha) * 100)
        
        return mean_stat, lower, upper
    
    @staticmethod
    def bootstrap_ne_ratio(pi_comp: np.ndarray, pi_A: np.ndarray,
                          D_comp: float, D_A: float,
                          n_bootstraps: int = 1000) -> Tuple[float, float, float]:
        """
        Bootstrap Nₑ ratio for a specific compartment.
        """
        def calculate_ratio(pi_comp_sample, pi_A_sample):
            pi_comp_mean = np.mean(pi_comp_sample)
            pi_A_mean = np.mean(pi_A_sample)
            return (pi_comp_mean * D_A) / (pi_A_mean * D_comp)
        
        # Generate bootstrap samples
        bootstrap_ratios = []
        n_comp = len(pi_comp)
        n_A = len(pi_A)
        
        for _ in range(n_bootstraps):
            pi_comp_boot = np.random.choice(pi_comp, size=n_comp, replace=True)
            pi_A_boot = np.random.choice(pi_A, size=n_A, replace=True)
            ratio = calculate_ratio(pi_comp_boot, pi_A_boot)
            bootstrap_ratios.append(ratio)
        
        bootstrap_ratios = np.array(bootstrap_ratios)
        mean_ratio = np.mean(bootstrap_ratios)
        lower_ci = np.percentile(bootstrap_ratios, 2.5)
        upper_ci = np.percentile(bootstrap_ratios, 97.5)
        
        return mean_ratio, lower_ci, upper_ci