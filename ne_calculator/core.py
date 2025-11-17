"""
Core Nₑ ratio calculation logic with bootstrap support.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
import logging
from .utils import StatisticsCalculator

logger = logging.getLogger(__name__)


class NECalculator:
    """
    Main class for calculating effective population size ratios with bootstrap support.
    """
    
    def __init__(self, ploidy_config: Optional[Dict] = None):
        """
        Initialize Nₑ calculator.
        
        Args:
            ploidy_config: Dictionary specifying ploidy for each compartment
        """
        self.ploidy_config = ploidy_config or {
            'A': 2,  # Autosomes
            'X': 1.5,  # X chromosome (haploid in males, diploid in females)
            'Y': 1,   # Y chromosome
            'MT': 1   # Mitochondrial DNA
        }
    
    def compute_ne_ratios(self, diversity_df: pd.DataFrame, 
                         divergence_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute Nₑ ratios using diversity and divergence data.
        """
        # Merge data on chromosome/compartment
        merged_df = pd.merge(
            diversity_df, divergence_df, 
            on=['chromosome', 'compartment', 'species'], 
            suffixes=('_pi', '_D')
        )
        
        results = []
        
        for species in merged_df['species'].unique():
            species_data = merged_df[merged_df['species'] == species]
            
            for comp in ['X', 'Y', 'MT']:
                comp_data = species_data[species_data['compartment'] == comp]
                autosome_data = species_data[species_data['compartment'] == 'A']
                
                if len(comp_data) > 0 and len(autosome_data) > 0:
                    pi_A = autosome_data['pi'].mean()
                    D_A = autosome_data['D'].mean()
                    pi_comp = comp_data['pi'].mean()
                    D_comp = comp_data['D'].mean()
                    
                    ne_ratio = (pi_comp * D_A) / (pi_A * D_comp)
                    
                    results.append({
                        'species': species,
                        'compartment': comp,
                        'ne_ratio': ne_ratio,
                        'pi_comp': pi_comp,
                        'pi_A': pi_A,
                        'D_comp': D_comp,
                        'D_A': D_A,
                        'n_windows_comp': len(comp_data),
                        'n_windows_A': len(autosome_data)
                    })
        
        return pd.DataFrame(results)
    
    def bootstrap_confidence_intervals(self, diversity_df: pd.DataFrame,
                                     divergence_df: pd.DataFrame,
                                     n_bootstraps: int = 1000,
                                     confidence: float = 0.95) -> pd.DataFrame:
        """
        Calculate bootstrap confidence intervals for Nₑ ratios.
        
        Args:
            diversity_df: Diversity data with multiple windows per chromosome
            divergence_df: Divergence data
            n_bootstraps: Number of bootstrap iterations
            confidence: Confidence level (0.95 for 95% CI)
            
        Returns:
            DataFrame with bootstrap confidence intervals
        """
        # Merge data
        merged_df = pd.merge(
            diversity_df, divergence_df,
            on=['chromosome', 'compartment', 'species'],
            suffixes=('_pi', '_D')
        )
        
        bootstrap_results = []
        
        for species in merged_df['species'].unique():
            species_data = merged_df[merged_df['species'] == species]
            
            for comp in ['X', 'Y', 'MT']:
                comp_data = species_data[species_data['compartment'] == comp]
                autosome_data = species_data[species_data['compartment'] == 'A']
                
                if len(comp_data) >= 10 and len(autosome_data) >= 10:  # Minimum for bootstrap
                    bootstrap_estimates = []
                    
                    for _ in range(n_bootstraps):
                        # Resample with replacement
                        comp_bootstrap = comp_data.sample(n=len(comp_data), replace=True)
                        autosome_bootstrap = autosome_data.sample(n=len(autosome_data), replace=True)
                        
                        pi_A = autosome_bootstrap['pi'].mean()
                        D_A = autosome_bootstrap['D'].mean()
                        pi_comp = comp_bootstrap['pi'].mean()
                        D_comp = comp_bootstrap['D'].mean()
                        
                        ne_ratio = (pi_comp * D_A) / (pi_A * D_comp)
                        bootstrap_estimates.append(ne_ratio)
                    
                    # Calculate confidence intervals
                    bootstrap_estimates = np.array(bootstrap_estimates)
                    lower = np.percentile(bootstrap_estimates, (1 - confidence) * 50)
                    upper = np.percentile(bootstrap_estimates, (1 + confidence) * 50)
                    mean_estimate = np.mean(bootstrap_estimates)
                    
                    bootstrap_results.append({
                        'species': species,
                        'compartment': comp,
                        'ne_ratio_mean': mean_estimate,
                        'ne_ratio_lower': lower,
                        'ne_ratio_upper': upper,
                        'bootstrap_se': np.std(bootstrap_estimates),
                        'n_bootstraps': n_bootstraps
                    })
        
        return pd.DataFrame(bootstrap_results)
    
    def jackknife_confidence_intervals(self, diversity_df: pd.DataFrame,
                                     divergence_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate jackknife confidence intervals for Nₑ ratios.
        """
        # Implementation using the StatisticsCalculator
        pass