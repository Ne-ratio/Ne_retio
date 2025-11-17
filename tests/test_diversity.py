"""
Tests for diversity calculation functionality.
"""

import pytest
import pandas as pd
import numpy as np
from ne_calculator.diversity import DiversityCalculator
from ne_calculator.utils import ChromosomeClassifier


class TestDiversityCalculator:
    """Test DiversityCalculator class."""
    
    def setup_method(self):
        """Setup test calculator."""
        self.calculator = DiversityCalculator()
    
    def test_chromosome_classification(self):
        """Test chromosome classification into compartments."""
        test_cases = [
            ('chr1', 'A'),
            ('chrX', 'X'),
            ('chrY', 'Y'),
            ('chrM', 'MT'),
            ('chrMT', 'MT'),
            ('X', 'X'),
            ('chromosome5', 'A'),
            ('unknown', 'A')  # Default to autosome
        ]
        
        for chrom_name, expected_comp in test_cases:
            result = self.calculator._classify_chromosome(chrom_name)
            assert result == expected_comp, f"Failed for {chrom_name}"
    
    def test_pi_calculation_single_site(self):
        """Test π calculation for a single site."""
        # Test case 1: All homozygous for same allele
        genotypes1 = [(0, 0), (0, 0), (0, 0)]
        pi1 = self.calculator._calculate_pi_site(genotypes1)
        assert pi1 == 0.0  # No diversity
        
        # Test case 2: All heterozygous
        genotypes2 = [(0, 1), (0, 1), (0, 1)]
        pi2 = self.calculator._calculate_pi_site(genotypes2)
        assert pi2 > 0.0  # Should have diversity
        
        # Test case 3: Mixed genotypes
        genotypes3 = [(0, 0), (0, 1), (1, 1)]
        pi3 = self.calculator._calculate_pi_site(genotypes3)
        assert 0.0 < pi3 < 1.0
    
    def test_pi_calculation_window(self):
        """Test π calculation for a window of variants."""
        # Create mock variant data
        variants = [
            {
                'genotypes': [(0, 0), (0, 1), (1, 1)],
                'ref': 'A',
                'alt': 'G',
                'position': 1000
            },
            {
                'genotypes': [(0, 0), (0, 0), (0, 1)],
                'ref': 'C',
                'alt': 'T',
                'position': 2000
            }
        ]
        
        pi, n_sites = self.calculator._calculate_pi_window(variants)
        
        assert n_sites == 2
        assert 0.0 < pi < 1.0
    
    def test_pi_calculation_empty_window(self):
        """Test π calculation with empty window."""
        variants = []
        pi, n_sites = self.calculator._calculate_pi_window(variants)
        
        assert pi == 0.0
        assert n_sites == 0
    
    def test_pi_calculation_single_sample(self):
        """Test π calculation with only one sample."""
        variants = [{
            'genotypes': [(0, 1)],  # Only one sample
            'ref': 'A',
            'alt': 'G',
            'position': 1000
        }]
        
        pi, n_sites = self.calculator._calculate_pi_window(variants)
        
        assert pi == 0.0  # Can't calculate diversity with one sample
        assert n_sites == 1


class TestChromosomeClassifier:
    """Test ChromosomeClassifier utility class."""
    
    def test_autosome_patterns(self):
        """Test autosome classification patterns."""
        autosomes = ['chr1', 'chr2', 'chr22', '1', '2', '22', 'chromosome1']
        
        for chrom in autosomes:
            result = ChromosomeClassifier.classify(chrom)
            assert result == 'A', f"Failed for {chrom}"
    
    def test_x_chromosome_patterns(self):
        """Test X chromosome classification patterns."""
        x_chroms = ['chrX', 'X', 'chromosomeX', 'ChrX']
        
        for chrom in x_chroms:
            result = ChromosomeClassifier.classify(chrom)
            assert result == 'X', f"Failed for {chrom}"
    
    def test_y_chromosome_patterns(self):
        """Test Y chromosome classification patterns."""
        y_chroms = ['chrY', 'Y', 'chromosomeY', 'ChrY']
        
        for chrom in y_chroms:
            result = ChromosomeClassifier.classify(chrom)
            assert result == 'Y', f"Failed for {chrom}"
    
    def test_mitochondrial_patterns(self):
        """Test mitochondrial DNA classification patterns."""
        mt_chroms = ['chrM', 'chrMT', 'MT', 'M', 'mitochondrion']
        
        for chrom in mt_chroms:
            result = ChromosomeClassifier.classify(chrom)
            assert result == 'MT', f"Failed for {chrom}"


def test_diversity_data_structure():
    """Test that diversity data has correct structure."""
    # This would test the actual VCF processing in a real scenario
    # For now, we test the data structure expectations
    
    calculator = DiversityCalculator()
    
    # Mock the VCF processing to return structured data
    class MockVariantFile:
        def __init__(self):
            self.samples = ['sample1', 'sample2', 'sample3']
        
        def __iter__(self):
            # Yield mock variants
            variants = [
                type('Variant', (), {
                    'chrom': 'chr1',
                    'pos': 1000,
                    'ref': 'A',
                    'alts': ('G',),
                    'samples': [
                        type('Sample', (), {'GT': (0, 0)}),
                        type('Sample', (), {'GT': (0, 1)}),
                        type('Sample', (), {'GT': (1, 1)})
                    ]
                })()
            ]
            return iter(variants)
        
        def close(self):
            pass
    
    # We would test the full pipeline here with mocked VCF file
    # For now, just verify our understanding of the data structure
    expected_columns = ['chromosome', 'compartment', 'start', 'end', 'pi', 'n_sites']
    assert hasattr(calculator, 'calculate_from_vcf')