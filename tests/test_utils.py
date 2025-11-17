"""
Tests for utility functions.
"""

import pytest
import numpy as np
import pandas as pd
from ne_calculator.utils import StatisticsCalculator, BootstrapCalculator


class TestStatisticsCalculator:
    """Test StatisticsCalculator class."""
    
    def test_jackknife_variance_basic(self):
        """Test basic jackknife variance calculation."""
        # Simple test case
        estimates = [1.0, 2.0, 3.0, 4.0, 5.0]
        mean, se = StatisticsCalculator.jackknife_variance(estimates)
        
        assert abs(mean - 3.0) < 0.001
        assert se > 0
    
    def test_jackknife_variance_single_value(self):
        """Test jackknife with single value."""
        estimates = [5.0]
        mean, se = StatisticsCalculator.jackknife_variance(estimates)
        
        assert mean == 5.0
        assert se == 0.0
    
    def test_jackknife_variance_two_values(self):
        """Test jackknife with two values."""
        estimates = [2.0, 4.0]
        mean, se = StatisticsCalculator.jackknife_variance(estimates)
        
        assert mean == 3.0
        assert se > 0


class TestBootstrapCalculator:
    """Test BootstrapCalculator class."""
    
    def setup_method(self):
        """Setup test data."""
        self.bootstrap_calc = BootstrapCalculator()
        np.random.seed(42)  # For reproducible tests
    
    def test_bootstrap_ci_basic(self):
        """Test basic bootstrap confidence interval calculation."""
        data = np.random.normal(0, 1, 100)
        
        def mean_statistic(x):
            return np.mean(x)
        
        mean, lower, upper = self.bootstrap_calc.bootstrap_ci(
            data, mean_statistic, n_bootstraps=100
        )
        
        assert lower < mean < upper
        assert abs(mean - np.mean(data)) < 0.5  # Should be close to true mean
    
    def test_bootstrap_ne_ratio(self):
        """Test bootstrap NE ratio calculation."""
        # Create test data
        pi_comp = np.random.normal(0.001, 0.0001, 50)
        pi_A = np.random.normal(0.0015, 0.0001, 50)
        D_comp = 0.012
        D_A = 0.015
        
        mean_ratio, lower_ci, upper_ci = self.bootstrap_calc.bootstrap_ne_ratio(
            pi_comp, pi_A, D_comp, D_A, n_bootstraps=100
        )
        
        # Check that results are reasonable
        assert lower_ci < mean_ratio < upper_ci
        assert mean_ratio > 0
    
    def test_bootstrap_ci_with_random_state(self):
        """Test that random state produces reproducible results."""
        data = np.array([1, 2, 3, 4, 5])
        
        def stat_func(x):
            return np.mean(x)
        
        # First call with random state
        mean1, lower1, upper1 = self.bootstrap_calc.bootstrap_ci(
            data, stat_func, n_bootstraps=50, random_state=123
        )
        
        # Second call with same random state
        mean2, lower2, upper2 = self.bootstrap_calc.bootstrap_ci(
            data, stat_func, n_bootstraps=50, random_state=123
        )
        
        # Should be identical
        assert mean1 == mean2
        assert lower1 == lower2
        assert upper1 == upper2


def test_data_validation_integration():
    """Integration test for data validation with real data structures."""
    from ne_calculator.utils import DataValidator
    
    validator = DataValidator()
    
    # Test with realistic data structure
    valid_diversity = pd.DataFrame({
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'] * 10,
        'compartment': ['A', 'X', 'Y', 'MT'] * 10,
        'start': list(range(0, 2000000, 50000)),
        'end': list(range(50000, 2050000, 50000)),
        'pi': np.random.uniform(0.0001, 0.003, 40),
        'species': ['test_species'] * 40
    })
    
    valid_divergence = pd.DataFrame({
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'],
        'compartment': ['A', 'X', 'Y', 'MT'],
        'D': [0.015, 0.012, 0.018, 0.025],
        'species': ['test_species'] * 4
    })
    
    # Should not raise exceptions
    assert validator.validate_diversity_data(valid_diversity) is True
    assert validator.validate_divergence_data(valid_divergence) is True


# Run tests with: pytest tests/ -v