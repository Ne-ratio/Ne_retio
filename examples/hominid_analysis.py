"""
Example analysis of hominid data (humans, chimpanzees, bonobos).
"""

import pandas as pd
import numpy as np
from ne_calculator import NECalculator

def analyze_hominid_data():
    """Analyze hominid genomic data."""
    print("Hominid Nₑ Ratio Analysis")
    
    # Example data structure for hominids
    species_data = {
        'human': {
            'diversity': {'A': 0.0012, 'X': 0.0009, 'Y': 0.0003, 'MT': 0.002},
            'divergence': {'A': 0.015, 'X': 0.012, 'Y': 0.018, 'MT': 0.025}
        },
        'chimp': {
            'diversity': {'A': 0.0015, 'X': 0.0011, 'Y': 0.0001, 'MT': 0.0028},
            'divergence': {'A': 0.015, 'X': 0.012, 'Y': 0.020, 'MT': 0.028}
        },
        'bonobo': {
            'diversity': {'A': 0.0014, 'X': 0.0010, 'Y': 0.0002, 'MT': 0.0025},
            'divergence': {'A': 0.015, 'X': 0.012, 'Y': 0.019, 'MT': 0.026}
        }
    }
    
    calculator = NECalculator()
    results = []
    
    for species, data in species_data.items():
        # Create DataFrames
        diversity_list = []
        divergence_list = []
        
        for comp in ['A', 'X', 'Y', 'MT']:
            diversity_list.append({
                'chromosome': f'chr{comp}',
                'compartment': comp,
                'pi': data['diversity'][comp]
            })
            
            divergence_list.append({
                'chromosome': f'chr{comp}',
                'compartment': comp,
                'D': data['divergence'][comp]
            })
        
        diversity_df = pd.DataFrame(diversity_list)
        divergence_df = pd.DataFrame(divergence_list)
        
        # Calculate ratios
        ne_ratios = calculator.compute_ne_ratios(diversity_df, divergence_df)
        ne_ratios['species'] = species
        
        results.append(ne_ratios)
    
    # Combine results
    all_results = pd.concat(results, ignore_index=True)
    
    # Pivot for comparison
    pivot_df = all_results.pivot_table(
        index='species', 
        columns='compartment', 
        values='ne_ratio'
    )
    
    print("\nNₑ Ratios by Species:")
    print(pivot_df)
    
    return all_results, pivot_df

if __name__ == "__main__":
    results, pivot = analyze_hominid_data()