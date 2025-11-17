"""
Example usage of NE Ratio Calculator.
"""

import pandas as pd
import matplotlib.pyplot as plt
from ne_calculator import NECalculator, DiversityCalculator

def basic_usage():
    """Basic usage example."""
    print("NE Ratio Calculator - Basic Usage Example")
    
    # Initialize calculator
    calculator = NECalculator()
    
    # Example data (in practice, load from VCF files)
    diversity_data = {
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'],
        'compartment': ['A', 'X', 'Y', 'MT'],
        'pi': [0.0012, 0.0009, 0.0001, 0.0025]
    }
    
    divergence_data = {
        'chromosome': ['chr1', 'chrX', 'chrY', 'chrMT'],
        'compartment': ['A', 'X', 'Y', 'MT'],
        'D': [0.015, 0.012, 0.020, 0.030]
    }
    
    diversity_df = pd.DataFrame(diversity_data)
    divergence_df = pd.DataFrame(divergence_data)
    
    # Calculate Nₑ ratios
    results = calculator.compute_ne_ratios(diversity_df, divergence_df)
    print("\nNₑ Ratio Results:")
    print(results)
    
    return results

def plot_results(results):
    """Plot Nₑ ratio results."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    compartments = results['compartment']
    ratios = results['ne_ratio']
    
    colors = {'A': 'blue', 'X': 'red', 'Y': 'green', 'MT': 'orange'}
    
    for i, (comp, ratio) in enumerate(zip(compartments, ratios)):
        ax.bar(i, ratio, color=colors.get(comp, 'gray'), label=comp)
    
    ax.set_xlabel('Genomic Compartment')
    ax.set_ylabel('Nₑ Ratio (relative to autosomes)')
    ax.set_title('Effective Population Size Ratios')
    ax.set_xticks(range(len(compartments)))
    ax.set_xticklabels(compartments)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('ne_ratios_example.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    results = basic_usage()
    plot_results(results)