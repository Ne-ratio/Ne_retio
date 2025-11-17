"""
Pytest configuration and fixtures for NE ratio calculator tests.
"""

import pytest
import pandas as pd
import numpy as np
from ne_calculator import NECalculator, DiversityCalculator


@pytest.fixture
def sample_diversity_data():
    """Fixture providing sample diversity data for tests."""
    return pd.DataFrame({
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'] * 5,
        'compartment': ['A', 'X', 'Y', 'MT'] * 5,
        'start': list(range(0, 1000000, 50000)),
        'end': list(range(50000, 1050000, 50000)),
        'pi': np.random.uniform(0.0005, 0.0025, 20),
        'n_sites': np.random.randint(20, 50, 20),
        'species': ['test_species'] * 20
    })


@pytest.fixture
def sample_divergence_data():
    """Fixture providing sample divergence data for tests."""
    return pd.DataFrame({
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'],
        'compartment': ['A', 'X', 'Y', 'MT'],
        'D': [0.015, 0.012, 0.018, 0.025],
        'species': ['test_species'] * 4
    })


@pytest.fixture
def ne_calculator():
    """Fixture providing NECalculator instance."""
    return NECalculator()


@pytest.fixture
def diversity_calculator():
    """Fixture providing DiversityCalculator instance."""
    return DiversityCalculator()


@pytest.fixture
def multi_species_data():
    """Fixture providing multi-species test data."""
    species = ['human', 'chimp', 'bonobo']
    data_frames = []
    
    for species_name in species:
        # Diversity data
        div_df = pd.DataFrame({
            'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'] * 5,
            'compartment': ['A', 'X', 'Y', 'MT'] * 5,
            'pi': np.random.uniform(0.0005, 0.0025, 20),
            'species': [species_name] * 20
        })
        data_frames.append(div_df)
    
    diversity_data = pd.concat(data_frames, ignore_index=True)
    
    # Divergence data
    divergence_data = pd.DataFrame({
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'] * 3,
        'compartment': ['A', 'X', 'Y', 'MT'] * 3,
        'D': np.random.uniform(0.010, 0.020, 12),
        'species': [sp for sp in species for _ in range(4)]
    })
    
    return diversity_data, divergence_data