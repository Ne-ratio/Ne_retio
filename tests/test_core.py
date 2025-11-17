"""
Tests for core NE ratio calculation functionality.
"""

import pytest
import pandas as pd
import numpy as np
from ne_calculator.core import NECalculator
from ne_calculator.utils import DataValidator


class TestNECalculator:
    """Test NECalculator class."""
    
    def setup_method(self):
        """Setup test data."""
        self.calculator = NECalculator()
        
        # Create test diversity data
        self.diversity_data = pd.DataFrame({
            'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'] * 5,
            'compartment': ['A', 'X', 'Y', 'MT'] * 5,
            'pi': [0.0012, 0.0009, 0.0003, 0.0020] * 5,
            'species': ['human'] * 20
        })
        
        # Create test divergence data
        self.divergence_data = pd.DataFrame({
            'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'],
            'compartment': ['A', 'X', 'Y', 'MT'],
            'D': [0.015, 0.012, 0.018, 0.025],
            'species': ['human'] * 4
        })
    
    def test_compute_ne_ratios_basic(self):
        """Test basic NE ratio calculation."""
        results = self.calculator.compute_ne_ratios(
            self.diversity_data, self.divergence_data
        )
        
        # Check structure
        expected_columns = ['species', 'compartment', 'ne_ratio', 'pi_comp', 
                          'pi_A', 'D_comp', 'D_A', 'n_windows_comp', 'n_windows_A']
        assert all(col in results.columns for col in expected_columns)
        
        # Check values
        assert len(results) == 3  # X, Y, MT compartments
        assert all(results['ne_ratio'] > 0)
        
        # Specific check for X chromosome ratio
        x_ratio = results[results['compartment'] == 'X']['ne_ratio'].iloc[0]
        expected_x_ratio = (0.0009 * 0.015) / (0.0012 * 0.012)
        assert abs(x_ratio - expected_x_ratio) < 0.001
    
    def test_compute_ne_ratios_multiple_species(self):
        """Test NE ratio calculation with multiple species."""
        # Add chimpanzee data
        chimp_diversity = self.diversity_data.copy()
        chimp_diversity['species'] = 'chimp'
        chimp_diversity['pi'] = [0.0015, 0.0011, 0.0001, 0.0028] * 5
        
        chimp_divergence = self.divergence_data.copy()
        chimp_divergence['species'] = 'chimp'
        chimp_divergence['D'] = [0.015, 0.012, 0.020, 0.028]
        
        combined_diversity = pd.concat([self.diversity_data, chimp_diversity])
        combined_divergence = pd.concat([self.divergence_data, chimp_divergence])
        
        results = self.calculator.compute_ne_ratios(
            combined_diversity, combined_divergence
        )
        
        # Check both species are present
        assert set(results['species']) == {'human', 'chimp'}
        assert len(results) == 6  # 3 compartments × 2 species
        
        # Check chimp has lower Y chromosome ratio (expected for polygynous species)
        human_y = results[(results['species'] == 'human') & 
                         (results['compartment'] == 'Y')]['ne_ratio'].iloc[0]
        chimp_y = results[(results['species'] == 'chimp') & 
                         (results['compartment'] == 'Y')]['ne_ratio'].iloc[0]
        
        # Chimps should have lower Y chromosome effective population size
        assert chimp_y < human_y
    
    def test_bootstrap_confidence_intervals(self):
        """Test bootstrap confidence interval calculation."""
        results = self.calculator.bootstrap_confidence_intervals(
            self.diversity_data, self.divergence_data,
            n_bootstraps=100,  # Use fewer for faster tests
            confidence=0.95
        )
        
        # Check structure
        expected_columns = ['species', 'compartment', 'ne_ratio_mean', 
                          'ne_ratio_lower', 'ne_ratio_upper', 'bootstrap_se', 'n_bootstraps']
        assert all(col in results.columns for col in expected_columns)
        
        # Check confidence intervals are reasonable
        for _, row in results.iterrows():
            assert row['ne_ratio_lower'] < row['ne_ratio_mean'] < row['ne_ratio_upper']
            assert row['bootstrap_se'] > 0
    
    def test_insufficient_data(self):
        """Test behavior with insufficient data for bootstrap."""
        # Create minimal data that shouldn't support bootstrap
        minimal_diversity = self.diversity_data.head(5)  # Only 5 windows
        minimal_divergence = self.divergence_data
        
        results = self.calculator.bootstrap_confidence_intervals(
            minimal_diversity, minimal_divergence,
            n_bootstraps=10
        )
        
        # Should return empty DataFrame or handle gracefully
        assert isinstance(results, pd.DataFrame)


class TestDataValidator:
    """Test DataValidator class."""
    
    def setup_method(self):
        """Setup test data."""
        self.validator = DataValidator()
        
        self.valid_diversity = pd.DataFrame({
            'chromosome': ['chr1', 'chrX', 'chrY'],
            'compartment': ['A', 'X', 'Y'],
            'start': [0, 0, 0],
            'end': [50000, 50000, 50000],
            'pi': [0.001, 0.002, 0.003]
        })
        
        self.valid_divergence = pd.DataFrame({
            'chromosome': ['chr1', 'chrX', 'chrY'],
            'compartment': ['A', 'X', 'Y'],
            'D': [0.01, 0.02, 0.03]
        })
    
    def test_validate_diversity_data_valid(self):
        """Test validation of valid diversity data."""
        assert self.validator.validate_diversity_data(self.valid_diversity) is True
    
    def test_validate_diversity_data_missing_columns(self):
        """Test validation with missing columns."""
        invalid_data = self.valid_diversity.drop('pi', axis=1)
        
        with pytest.raises(ValueError, match="Missing required columns"):
            self.validator.validate_diversity_data(invalid_data)
    
    def test_validate_diversity_data_negative_pi(self):
        """Test validation with negative π values."""
        invalid_data = self.valid_diversity.copy()
        invalid_data.loc[0, 'pi'] = -0.001
        
        with pytest.raises(ValueError, match="Negative π values"):
            self.validator.validate_diversity_data(invalid_data)
    
    def test_validate_diversity_data_invalid_compartment(self):
        """Test validation with invalid compartment values."""
        invalid_data = self.valid_diversity.copy()
        invalid_data.loc[0, 'compartment'] = 'Z'  # Invalid compartment
        
        with pytest.raises(ValueError, match="Invalid compartments"):
            self.validator.validate_diversity_data(invalid_data)
    
    def test_validate_divergence_data_valid(self):
        """Test validation of valid divergence data."""
        assert self.validator.validate_divergence_data(self.valid_divergence) is True
    
    def test_validate_divergence_data_negative_divergence(self):
        """Test validation with negative divergence values."""
        invalid_data = self.valid_divergence.copy()
        invalid_data.loc[0, 'D'] = -0.01
        
        with pytest.raises(ValueError, match="Negative divergence values"):
            self.validator.validate_divergence_data(invalid_data)


def test_integration_with_test_data():
    """Integration test with simulated test data."""
    calculator = NECalculator()
    
    # Load test data
    diversity_df = pd.read_csv('test_data/diversity_simulated.csv')
    divergence_df = pd.read_csv('test_data/divergence_simulated.csv')
    
    # Calculate ratios
    results = calculator.compute_ne_ratios(diversity_df, divergence_df)
    
    # Basic checks
    assert len(results) > 0
    assert all(results['ne_ratio'] > 0)
    assert set(results['species']) == {'human', 'chimp', 'bonobo'}
    assert set(results['compartment']) == {'X', 'Y', 'MT'}
    
    # Check that bootstrap works with sufficient data
    bootstrap_results = calculator.bootstrap_confidence_intervals(
        diversity_df, divergence_df, n_bootstraps=50  # Small number for speed
    )
    assert len(bootstrap_results) > 0