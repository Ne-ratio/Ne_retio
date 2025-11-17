#!/usr/bin/env python3
"""
Script to generate simulated VCF data for testing NE Ratio Calculator.
"""

import numpy as np
import pandas as pd
import random
from typing import List, Dict

def generate_simulated_diversity_data():
    """
    Generate simulated diversity data with realistic patterns.
    
    Returns:
        DataFrame with simulated π values across genomic compartments
    """
    np.random.seed(42)  # For reproducible results
    
    species_data = {
        'human': {'A_mean': 0.00125, 'X_mean': 0.00092, 'Y_mean': 0.00035, 'MT_mean': 0.00215},
        'chimp': {'A_mean': 0.00145, 'X_mean': 0.00112, 'Y_mean': 0.00012, 'MT_mean': 0.00245},
        'bonobo': {'A_mean': 0.00138, 'X_mean': 0.00105, 'Y_mean': 0.00025, 'MT_mean': 0.00235}
    }
    
    compartments = {
        'A': ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'],
        'X': ['chrX'],
        'Y': ['chrY'], 
        'MT': ['chrMT']
    }
    
    all_data = []
    
    for species, means in species_data.items():
        for comp, chroms in compartments.items():
            mean_pi = means[f'{comp}_mean']
            
            for chrom in chroms:
                # Generate 20-30 windows per chromosome
                n_windows = random.randint(20, 30)
                
                for i in range(n_windows):
                    start = i * 50000
                    end = start + 50000
                    
                    # Add some realistic variation
                    pi_variation = np.random.normal(0, mean_pi * 0.1)  # 10% variation
                    pi_value = max(0.0001, mean_pi + pi_variation)
                    
                    # Number of sites correlates with π
                    n_sites = int(40 * (pi_value / 0.001) + np.random.randint(-5, 5))
                    n_sites = max(10, min(60, n_sites))  # Bound between 10-60
                    
                    all_data.append({
                        'chromosome': chrom,
                        'compartment': comp,
                        'start': start,
                        'end': end,
                        'pi': round(pi_value, 6),
                        'n_sites': n_sites,
                        'species': species
                    })
    
    return pd.DataFrame(all_data)

def generate_simulated_divergence_data():
    """
    Generate simulated divergence data.
    """
    species_data = {
        'human': {'A': 0.0152, 'X': 0.0121, 'Y': 0.0185, 'MT': 0.0248},
        'chimp': {'A': 0.0150, 'X': 0.0118, 'Y': 0.0192, 'MT': 0.0275},
        'bonobo': {'A': 0.0148, 'X': 0.0115, 'Y': 0.0188, 'MT': 0.0262}
    }
    
    all_data = []
    
    for species, div_vals in species_data.items():
        for comp, D in div_vals.items():
            if comp == 'A':
                # Multiple autosomes
                for i in range(1, 6):
                    all_data.append({
                        'chromosome': f'chr{i}',
                        'compartment': comp,
                        'D': D + np.random.normal(0, 0.0002),
                        'species': species
                    })
            else:
                # Single chromosome for X, Y, MT
                all_data.append({
                    'chromosome': f'chr{comp}',
                    'compartment': comp,
                    'D': D,
                    'species': species
                })
    
    divergence_df = pd.DataFrame(all_data)
    divergence_df['D'] = divergence_df['D'].round(4)
    return divergence_df

def calculate_ne_ratios(diversity_df, divergence_df):
    """
    Calculate Nₑ ratios from diversity and divergence data.
    """
    results = []
    
    for species in diversity_df['species'].unique():
        species_diversity = diversity_df[diversity_df['species'] == species]
        species_divergence = divergence_df[divergence_df['species'] == species]
        
        # Get autosomal reference values
        autosome_diversity = species_diversity[species_diversity['compartment'] == 'A']
        autosome_divergence = species_divergence[species_divergence['compartment'] == 'A']
        
        pi_A = autosome_diversity['pi'].mean()
        D_A = autosome_divergence['D'].mean()
        
        for comp in ['X', 'Y', 'MT']:
            comp_diversity = species_diversity[species_diversity['compartment'] == comp]
            comp_divergence = species_divergence[species_divergence['compartment'] == comp]
            
            if len(comp_diversity) > 0 and len(comp_divergence) > 0:
                pi_comp = comp_diversity['pi'].mean()
                D_comp = comp_divergence['D'].mean()
                
                ne_ratio = (pi_comp * D_A) / (pi_A * D_comp)
                
                # Add some confidence intervals (simulated)
                ci_width = ne_ratio * 0.03  # 3% CI width
                lower = max(0.001, ne_ratio - ci_width)
                upper = ne_ratio + ci_width
                
                results.append({
                    'species': species,
                    'compartment': comp,
                    'ne_ratio': round(ne_ratio, 3),
                    'ne_ratio_lower': round(lower, 3),
                    'ne_ratio_upper': round(upper, 3)
                })
    
    return pd.DataFrame(results)

def main():
    """Generate all simulated test data."""
    print("Generating simulated test data...")
    
    # Generate data
    diversity_df = generate_simulated_diversity_data()
    divergence_df = generate_simulated_divergence_data()
    ne_ratios_df = calculate_ne_ratios(diversity_df, divergence_df)
    
    # Save to files
    diversity_df.to_csv('diversity_simulated.csv', index=False)
    divergence_df.to_csv('divergence_simulated.csv', index=False)
    ne_ratios_df.to_csv('ne_ratios_simulated.csv', index=False)
    
    print("Data generation complete!")
    print(f"Diversity data: {len(diversity_df)} windows")
    print(f"Divergence data: {len(divergence_df)} entries") 
    print(f"Nₑ ratios: {len(ne_ratios_df)} calculations")
    
    # Print summary
    print("\nSpecies summary:")
    for species in diversity_df['species'].unique():
        sp_data = diversity_df[diversity_df['species'] == species]
        print(f"\n{species}:")
        for comp in ['A', 'X', 'Y', 'MT']:
            comp_data = sp_data[sp_data['compartment'] == comp]
            if len(comp_data) > 0:
                mean_pi = comp_data['pi'].mean()
                n_windows = len(comp_data)
                print(f"  {comp}: π = {mean_pi:.4f} ({n_windows} windows)")

if __name__ == '__main__':
    main()