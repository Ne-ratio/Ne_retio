"""
Nucleotide diversity calculation from VCF files.
"""

import numpy as np
import pandas as pd
import pysam
from typing import List, Dict, Optional, Tuple
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


class DiversityCalculator:
    """Calculate nucleotide diversity (π) from VCF files."""
    
    def __init__(self):
        self.compartments = ['A', 'X', 'Y', 'MT']
    
    def calculate_from_vcf(self, vcf_file: str, window_size: int = 50000,
                          samples: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Calculate π from VCF file in non-overlapping windows.
        
        Args:
            vcf_file: Path to VCF file
            window_size: Window size in bp
            samples: List of samples to include
            
        Returns:
            DataFrame with diversity calculations
        """
        try:
            vcf = pysam.VariantFile(vcf_file)
        except Exception as e:
            logger.error(f"Error opening VCF file: {e}")
            raise
        
        # Filter samples if provided
        if samples:
            vcf.subset_samples(samples)
        
        results = []
        chromosome_data = {}
        
        # First pass: collect variant data by chromosome
        logger.info("Reading VCF file...")
        for variant in tqdm(vcf):
            chrom = variant.chrom
            pos = variant.pos
            
            if chrom not in chromosome_data:
                chromosome_data[chrom] = []
            
            # Extract genotype data
            genotypes = []
            for sample in variant.samples:
                gt = variant.samples[sample]['GT']
                if None not in gt:  # Skip missing genotypes
                    genotypes.append(gt)
            
            if len(genotypes) > 0:
                chromosome_data[chrom].append({
                    'position': pos,
                    'genotypes': genotypes,
                    'ref': variant.ref,
                    'alt': variant.alts[0] if variant.alts else '.'
                })
        
        # Second pass: calculate π in windows
        logger.info("Calculating nucleotide diversity...")
        for chrom, variants in chromosome_data.items():
            if not variants:
                continue
                
            comp = self._classify_chromosome(chrom)
            variants.sort(key=lambda x: x['position'])
            
            positions = [v['position'] for v in variants]
            min_pos, max_pos = min(positions), max(positions)
            
            for window_start in range(min_pos, max_pos, window_size):
                window_end = window_start + window_size
                window_variants = [
                    v for v in variants 
                    if window_start <= v['position'] < window_end
                ]
                
                if len(window_variants) >= 10:  # Minimum sites threshold
                    pi = self._calculate_pi_window(window_variants)
                    
                    results.append({
                        'chromosome': chrom,
                        'compartment': comp,
                        'start': window_start,
                        'end': window_end,
                        'pi': pi,
                        'n_sites': len(window_variants),
                        'n_samples': len(window_variants[0]['genotypes']) if window_variants else 0
                    })
        
        vcf.close()
        return pd.DataFrame(results)
    
    def _classify_chromosome(self, chrom_name: str) -> str:
        """Classify chromosome into genomic compartment."""
        chrom_lower = chrom_name.lower()
        
        if 'x' in chrom_lower:
            return 'X'
        elif 'y' in chrom_lower:
            return 'Y'
        elif 'm' in chrom_lower or 'mt' in chrom_lower:
            return 'MT'
        else:
            return 'A'
    
    def _calculate_pi_window(self, variants: List[Dict]) -> float:
        """Calculate π for a window of variants."""
        if not variants:
            return 0.0
        
        n_samples = len(variants[0]['genotypes'])
        if n_samples < 2:
            return 0.0
        
        total_pi = 0.0
        total_sites = 0
        
        for variant in variants:
            genotypes = variant['genotypes']
            site_pi = self._calculate_pi_site(genotypes)
            total_pi += site_pi
            total_sites += 1
        
        return total_pi / total_sites if total_sites > 0 else 0.0
    
    def _calculate_pi_site(self, genotypes: List[Tuple]) -> float:
        """Calculate π for a single site."""
        n_samples = len(genotypes)
        if n_samples < 2:
            return 0.0
        
        # Count allele frequencies
        allele_counts = {}
        total_alleles = 0
        
        for gt in genotypes:
            for allele in gt:
                if allele is not None:
                    allele_counts[allele] = allele_counts.get(allele, 0) + 1
                    total_alleles += 1
        
        if total_alleles == 0:
            return 0.0
        
        # Calculate expected heterozygosity
        pi = 1.0
        for count in allele_counts.values():
            freq = count / total_alleles
            pi -= freq ** 2
        
        return pi