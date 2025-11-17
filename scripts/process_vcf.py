#!/usr/bin/env python3
"""
VCF Processing Script for NE Ratio Calculator

This script processes VCF files to calculate nucleotide diversity (π)
in non-overlapping windows across genomic compartments.
"""

import argparse
import logging
import sys
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import pysam
from tqdm import tqdm

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from ne_calculator.diversity import DiversityCalculator
from ne_calculator.utils import ChromosomeClassifier, DataValidator


class VCFProcessor:
    """
    Process VCF files to calculate nucleotide diversity.
    """
    
    def __init__(self, window_size: int = 50000, min_sites: int = 10):
        self.window_size = window_size
        self.min_sites = min_sites
        self.diversity_calc = DiversityCalculator()
        self.logger = logging.getLogger(__name__)
    
    def parse_vcf_metadata(self, vcf_file: str) -> Dict:
        """
        Parse VCF metadata to get sample information and basic statistics.
        
        Args:
            vcf_file: Path to VCF file
            
        Returns:
            Dictionary with VCF metadata
        """
        try:
            vcf = pysam.VariantFile(vcf_file)
            metadata = {
                'samples': list(vcf.header.samples),
                'n_samples': len(vcf.header.samples),
                'contigs': list(vcf.header.contigs),
                'file_format': str(vcf.header.version)
            }
            vcf.close()
            return metadata
        except Exception as e:
            self.logger.error(f"Error reading VCF metadata: {e}")
            raise
    
    def calculate_diversity_from_vcf(self, vcf_file: str, 
                                   samples: Optional[List[str]] = None,
                                   chromosomes: Optional[List[str]] = None,
                                   max_variants: Optional[int] = None) -> pd.DataFrame:
        """
        Calculate nucleotide diversity from VCF file.
        
        Args:
            vcf_file: Path to VCF file
            samples: List of samples to include (None = all)
            chromosomes: List of chromosomes to process (None = all)
            max_variants: Maximum number of variants to process (for testing)
            
        Returns:
            DataFrame with diversity calculations
        """
        self.logger.info(f"Processing VCF file: {vcf_file}")
        
        try:
            # Open VCF file
            vcf = pysam.VariantFile(vcf_file)
            
            # Filter samples if provided
            if samples:
                available_samples = set(vcf.header.samples)
                selected_samples = [s for s in samples if s in available_samples]
                if len(selected_samples) != len(samples):
                    missing = set(samples) - available_samples
                    self.logger.warning(f"Missing samples: {missing}")
                vcf.subset_samples(selected_samples)
                self.logger.info(f"Processing {len(selected_samples)} samples")
            else:
                self.logger.info(f"Processing all {len(vcf.header.samples)} samples")
            
            # Initialize data structures
            chromosome_data = defaultdict(list)
            variant_count = 0
            
            # Process variants
            self.logger.info("Reading variants...")
            for variant in tqdm(vcf):
                if max_variants and variant_count >= max_variants:
                    self.logger.info(f"Reached maximum variant limit: {max_variants}")
                    break
                
                chrom = variant.chrom
                
                # Filter chromosomes if specified
                if chromosomes and chrom not in chromosomes:
                    continue
                
                pos = variant.pos
                
                # Skip multi-allelic sites for simplicity
                if len(variant.alts) > 1:
                    continue
                
                # Skip indels
                if any(len(alt) != len(variant.ref) for alt in variant.alts):
                    continue
                
                # Extract genotype data
                genotypes = []
                for sample in variant.samples:
                    gt = variant.samples[sample].get('GT', None)
                    if gt is not None and None not in gt:  # Skip missing genotypes
                        genotypes.append(gt)
                
                if len(genotypes) >= 2:  # Need at least 2 samples for diversity
                    chromosome_data[chrom].append({
                        'position': pos,
                        'genotypes': genotypes,
                        'ref': variant.ref,
                        'alt': variant.alts[0] if variant.alts else '.'
                    })
                    variant_count += 1
            
            vcf.close()
            self.logger.info(f"Processed {variant_count} variants across {len(chromosome_data)} chromosomes")
            
            # Calculate diversity in windows
            self.logger.info("Calculating nucleotide diversity in windows...")
            results = self._calculate_diversity_windows(chromosome_data)
            
            return results
            
        except Exception as e:
            self.logger.error(f"Error processing VCF file: {e}")
            if 'vcf' in locals():
                vcf.close()
            raise
    
    def _calculate_diversity_windows(self, chromosome_data: Dict) -> pd.DataFrame:
        """
        Calculate diversity in non-overlapping windows.
        
        Args:
            chromosome_data: Dictionary of chromosome -> variant data
            
        Returns:
            DataFrame with window-based diversity calculations
        """
        results = []
        
        for chrom, variants in chromosome_data.items():
            if not variants:
                continue
                
            comp = ChromosomeClassifier.classify(chrom)
            variants.sort(key=lambda x: x['position'])
            
            positions = [v['position'] for v in variants]
            min_pos, max_pos = min(positions), max(positions)
            
            self.logger.debug(f"Processing {chrom} ({comp}): {len(variants)} variants, "
                            f"positions {min_pos}-{max_pos}")
            
            n_windows = 0
            n_windows_with_data = 0
            
            for window_start in range(min_pos, max_pos, self.window_size):
                window_end = window_start + self.window_size
                window_variants = [
                    v for v in variants 
                    if window_start <= v['position'] < window_end
                ]
                
                n_windows += 1
                
                if len(window_variants) >= self.min_sites:
                    pi = self.diversity_calc._calculate_pi_window(window_variants)
                    
                    results.append({
                        'chromosome': chrom,
                        'compartment': comp,
                        'start': window_start,
                        'end': window_end,
                        'pi': pi,
                        'n_sites': len(window_variants),
                        'n_samples': len(window_variants[0]['genotypes']) if window_variants else 0
                    })
                    n_windows_with_data += 1
            
            self.logger.info(f"Chromosome {chrom}: {n_windows_with_data}/{n_windows} "
                           f"windows passed filters (≥{self.min_sites} sites)")
        
        return pd.DataFrame(results)
    
    def filter_low_quality_windows(self, diversity_df: pd.DataFrame, 
                                 min_sites: int = 10,
                                 max_missing: float = 0.5) -> pd.DataFrame:
        """
        Filter low-quality diversity windows.
        
        Args:
            diversity_df: Diversity DataFrame
            min_sites: Minimum number of sites per window
            max_missing: Maximum allowed proportion of missing data
            
        Returns:
            Filtered DataFrame
        """
        initial_count = len(diversity_df)
        
        # Filter by minimum sites
        filtered_df = diversity_df[diversity_df['n_sites'] >= min_sites]
        
        # Additional filtering could be added here (missing data, etc.)
        
        filtered_count = len(filtered_df)
        removed_count = initial_count - filtered_count
        
        self.logger.info(f"Filtered {removed_count} windows "
                       f"({filtered_count} remaining, {removed_count/initial_count*100:.1f}% removed)")
        
        return filtered_df
    
    def generate_summary_statistics(self, diversity_df: pd.DataFrame) -> Dict:
        """
        Generate summary statistics for diversity data.
        
        Args:
            diversity_df: Diversity DataFrame
            
        Returns:
            Dictionary with summary statistics
        """
        if diversity_df.empty:
            return {}
        
        summary = {
            'total_windows': len(diversity_df),
            'total_sites': diversity_df['n_sites'].sum(),
            'mean_pi': diversity_df['pi'].mean(),
            'std_pi': diversity_df['pi'].std(),
            'median_pi': diversity_df['pi'].median()
        }
        
        # Statistics by compartment
        compartment_stats = {}
        for comp in diversity_df['compartment'].unique():
            comp_data = diversity_df[diversity_df['compartment'] == comp]
            compartment_stats[comp] = {
                'n_windows': len(comp_data),
                'mean_pi': comp_data['pi'].mean(),
                'std_pi': comp_data['pi'].std(),
                'n_sites': comp_data['n_sites'].sum()
            }
        
        summary['compartments'] = compartment_stats
        
        return summary


def setup_logging(verbose: bool = False, log_file: Optional[str] = None):
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def main():
    """Main command line interface."""
    parser = argparse.ArgumentParser(
        description='Process VCF files to calculate nucleotide diversity',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with default parameters
  python process_vcf.py --vcf input.vcf.gz --output diversity.csv
  
  # Process specific chromosomes with custom window size
  python process_vcf.py --vcf input.vcf.gz --output diversity.csv \\
    --chromosomes chr1 chr2 chrX --window-size 100000
  
  # Process only specific samples
  python process_vcf.py --vcf input.vcf.gz --output diversity.csv \\
    --samples sample1 sample2 sample3
  
  # Generate detailed summary and plots
  python process_vcf.py --vcf input.vcf.gz --output diversity.csv \\
    --summary --plot --verbose
        """
    )
    
    # Required arguments
    parser.add_argument('--vcf', required=True, help='Input VCF file (can be gzipped)')
    parser.add_argument('--output', required=True, help='Output CSV file')
    
    # Processing options
    parser.add_argument('--window-size', type=int, default=50000,
                       help='Window size in base pairs (default: 50000)')
    parser.add_argument('--min-sites', type=int, default=10,
                       help='Minimum number of sites per window (default: 10)')
    parser.add_argument('--samples', nargs='+', 
                       help='Specific samples to include (default: all)')
    parser.add_argument('--chromosomes', nargs='+',
                       help='Specific chromosomes to process (default: all)')
    parser.add_argument('--max-variants', type=int,
                       help='Maximum number of variants to process (for testing)')
    
    # Output options
    parser.add_argument('--summary', action='store_true',
                       help='Generate summary statistics')
    parser.add_argument('--plot', action='store_true',
                       help='Generate diversity distribution plots')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')
    parser.add_argument('--log-file', help='Log file path')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose, args.log_file)
    logger = logging.getLogger(__name__)
    
    try:
        # Initialize processor
        processor = VCFProcessor(
            window_size=args.window_size,
            min_sites=args.min_sites
        )
        
        # Parse VCF metadata
        metadata = processor.parse_vcf_metadata(args.vcf)
        logger.info(f"VCF metadata: {metadata['n_samples']} samples, "
                   f"{len(metadata['contigs'])} contigs")
        
        # Calculate diversity
        diversity_df = processor.calculate_diversity_from_vcf(
            vcf_file=args.vcf,
            samples=args.samples,
            chromosomes=args.chromosomes,
            max_variants=args.max_variants
        )
        
        if diversity_df.empty:
            logger.warning("No diversity data generated. Check input parameters.")
            sys.exit(1)
        
        # Filter low-quality windows
        filtered_df = processor.filter_low_quality_windows(diversity_df)
        
        # Save results
        filtered_df.to_csv(args.output, index=False)
        logger.info(f"Diversity results saved to {args.output}")
        logger.info(f"Total windows: {len(filtered_df)}")
        
        # Generate summary statistics
        if args.summary:
            summary = processor.generate_summary_statistics(filtered_df)
            
            summary_file = Path(args.output).with_suffix('.summary.txt')
            with open(summary_file, 'w') as f:
                f.write("Diversity Analysis Summary\n")
                f.write("=========================\n\n")
                f.write(f"Input VCF: {args.vcf}\n")
                f.write(f"Window size: {args.window_size} bp\n")
                f.write(f"Minimum sites per window: {args.min_sites}\n")
                f.write(f"Total windows: {summary['total_windows']}\n")
                f.write(f"Total sites: {summary['total_sites']}\n")
                f.write(f"Mean π: {summary['mean_pi']:.6f}\n")
                f.write(f"Standard deviation π: {summary['std_pi']:.6f}\n\n")
                
                f.write("By Genomic Compartment:\n")
                for comp, stats in summary['compartments'].items():
                    f.write(f"  {comp}:\n")
                    f.write(f"    Windows: {stats['n_windows']}\n")
                    f.write(f"    Mean π: {stats['mean_pi']:.6f}\n")
                    f.write(f"    Sites: {stats['n_sites']}\n")
            
            logger.info(f"Summary saved to {summary_file}")
        
        # Generate plots
        if args.plot:
            try:
                plot_diversity_distribution(filtered_df, args.output)
                logger.info("Diversity plots generated")
            except Exception as e:
                logger.warning(f"Could not generate plots: {e}")
        
        logger.info("VCF processing completed successfully!")
        
    except Exception as e:
        logger.error(f"VCF processing failed: {e}")
        sys.exit(1)


def plot_diversity_distribution(diversity_df: pd.DataFrame, output_base: str):
    """
    Generate diversity distribution plots.
    
    Args:
        diversity_df: Diversity DataFrame
        output_base: Base path for output files
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    output_path = Path(output_base)
    
    # Set style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Plot 1: Diversity by compartment
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Boxplot of π by compartment
    sns.boxplot(data=diversity_df, x='compartment', y='pi', ax=axes[0,0])
    axes[0,0].set_title('Nucleotide Diversity by Genomic Compartment')
    axes[0,0].set_ylabel('π')
    axes[0,0].set_xlabel('Compartment')
    
    # Histogram of π values
    diversity_df['pi'].hist(bins=50, ax=axes[0,1])
    axes[0,1].set_title('Distribution of π Values')
    axes[0,1].set_xlabel('π')
    axes[0,1].set_ylabel('Frequency')
    
    # Sites vs diversity
    axes[1,0].scatter(diversity_df['n_sites'], diversity_df['pi'], alpha=0.6)
    axes[1,0].set_xlabel('Number of Sites per Window')
    axes[1,0].set_ylabel('π')
    axes[1,0].set_title('Sites vs Diversity')
    
    # Mean diversity by chromosome
    chrom_stats = diversity_df.groupby('chromosome')['pi'].agg(['mean', 'count'])
    chrom_stats = chrom_stats[chrom_stats['count'] >= 5]  # Only chromosomes with enough data
    if len(chrom_stats) > 0:
        chrom_stats['mean'].plot(kind='bar', ax=axes[1,1])
        axes[1,1].set_title('Mean π by Chromosome')
        axes[1,1].set_ylabel('Mean π')
        axes[1,1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plot_file = output_path.with_suffix('.plots.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create compartment-specific plots if we have multiple compartments
    compartments = diversity_df['compartment'].unique()
    if len(compartments) > 1:
        fig, axes = plt.subplots(1, len(compartments), figsize=(4*len(compartments), 4))
        if len(compartments) == 1:
            axes = [axes]
        
        for i, comp in enumerate(compartments):
            comp_data = diversity_df[diversity_df['compartment'] == comp]
            comp_data['pi'].hist(bins=30, ax=axes[i])
            axes[i].set_title(f'{comp} Compartment')
            axes[i].set_xlabel('π')
            axes[i].set_ylabel('Frequency')
        
        plt.tight_layout()
        comp_plot_file = output_path.with_suffix('.compartments.png')
        plt.savefig(comp_plot_file, dpi=300, bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    main()