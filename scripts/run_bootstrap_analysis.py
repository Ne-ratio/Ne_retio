#!/usr/bin/env python3
"""
Script for running bootstrap analysis on NE ratios with confidence intervals.
"""

import argparse
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent.parent))

from ne_calculator import NECalculator

def setup_logging(verbose: bool):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def plot_bootstrap_results(results_df, output_dir):
    """Plot bootstrap results with confidence intervals."""
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Nₑ ratios by species and compartment
    sns.barplot(data=results_df, x='compartment', y='ne_ratio_mean', 
                hue='species', ax=axes[0])
    axes[0].set_title('Nₑ Ratios by Genomic Compartment')
    axes[0].set_ylabel('Nₑ Ratio (relative to autosomes)')
    axes[0].legend(title='Species')
    
    # Plot 2: Confidence intervals
    for i, (_, row) in enumerate(results_df.iterrows()):
        x_pos = i % 3  # 3 compartments
        species_offset = {'human': -0.2, 'chimp': 0, 'bonobo': 0.2}
        offset = species_offset.get(row['species'], 0)
        
        axes[1].errorbar(
            x=x_pos + offset,
            y=row['ne_ratio_mean'],
            yerr=[[row['ne_ratio_mean'] - row['ne_ratio_lower']], 
                  [row['ne_ratio_upper'] - row['ne_ratio_mean']]],
            fmt='o',
            capsize=5,
            label=row['species'] if i % 3 == 0 else ""
        )
    
    axes[1].set_xticks([0, 1, 2])
    axes[1].set_xticklabels(['X', 'Y', 'MT'])
    axes[1].set_title('Nₑ Ratios with 95% Confidence Intervals')
    axes[1].set_ylabel('Nₑ Ratio')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = output_dir / 'bootstrap_results.png'
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.show()
    
    return plot_file

def main():
    parser = argparse.ArgumentParser(
        description='Run bootstrap analysis for NE ratios'
    )
    
    parser.add_argument('--diversity', required=True, 
                       help='Diversity CSV file')
    parser.add_argument('--divergence', required=True,
                       help='Divergence CSV file')
    parser.add_argument('--output-dir', default='bootstrap_results',
                       help='Output directory')
    parser.add_argument('--n-bootstraps', type=int, default=1000,
                       help='Number of bootstrap iterations')
    parser.add_argument('--confidence', type=float, default=0.95,
                       help='Confidence level')
    parser.add_argument('--plot', action='store_true',
                       help='Generate plots')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Load data
        logger.info("Loading data...")
        diversity_df = pd.read_csv(args.diversity)
        divergence_df = pd.read_csv(args.divergence)
        
        # Run bootstrap analysis
        logger.info(f"Running bootstrap analysis ({args.n_bootstraps} iterations)...")
        calculator = NECalculator()
        
        bootstrap_results = calculator.bootstrap_confidence_intervals(
            diversity_df, 
            divergence_df,
            n_bootstraps=args.n_bootstraps,
            confidence=args.confidence
        )
        
        # Save results
        results_file = output_dir / 'bootstrap_ne_ratios.csv'
        bootstrap_results.to_csv(results_file, index=False)
        logger.info(f"Bootstrap results saved to {results_file}")
        
        # Generate summary statistics
        summary_file = output_dir / 'bootstrap_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("Bootstrap Analysis Summary\n")
            f.write("=========================\n\n")
            f.write(f"Number of bootstrap iterations: {args.n_bootstraps}\n")
            f.write(f"Confidence level: {args.confidence}\n")
            f.write(f"Total calculations: {len(bootstrap_results)}\n\n")
            
            f.write("Nₑ Ratio Estimates with Confidence Intervals:\n")
            for _, row in bootstrap_results.iterrows():
                f.write(f"\n{row['species']} - {row['compartment']}:\n")
                f.write(f"  Mean: {row['ne_ratio_mean']:.4f}\n")
                f.write(f"  95% CI: [{row['ne_ratio_lower']:.4f}, {row['ne_ratio_upper']:.4f}]\n")
                f.write(f"  Bootstrap SE: {row['bootstrap_se']:.4f}\n")
        
        logger.info(f"Summary saved to {summary_file}")
        
        # Generate plots if requested
        if args.plot:
            logger.info("Generating plots...")
            plot_file = plot_bootstrap_results(bootstrap_results, output_dir)
            logger.info(f"Plot saved to {plot_file}")
        
        logger.info("Bootstrap analysis complete!")
        
    except Exception as e:
        logger.error(f"Bootstrap analysis failed: {e}")
        raise

if __name__ == '__main__':
    main()