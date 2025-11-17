# Estimating Effective Population Size Ratios Across Genomic Compartments

This repository contains scripts and data examples for the paper:

> **A phylogenetically informed framework for estimating effective population size ratios across genomic compartments**  
> Liangyu Lu, Wennan Dai, Yuhao Pan, Zheng Yan  
> *State Key Laboratory of Grassland Agro-Ecosystems & College of Ecology, Lanzhou University*

## Overview

We present a phylogenetic framework to estimate relative effective population size (*Nₑ*) ratios across genomic compartments (autosomes, X, Y, mtDNA) using interspecies divergence to calibrate mutation rates.

Key features:
- Uses diversity (π) and divergence (D) to cancel out mutation rate and time
- Provides unbiased *Nₑ* ratio estimates
- Applicable to non-model organisms

## File Structure

- `scripts/`: Python scripts for calculating diversity, divergence, and *Nₑ* ratios
- `data/`: Example input/output files
- `figures/`: Example output plot

## Usage

1. Install Dependencies

```bash
pip install -r requirements.txt
2. Compute Nucleotide Diversity (π)
Use compute_diversity.py to calculate π from VCF files:
bash
python scripts/compute_diversity.py --vcf <input.vcf> --windows 50000 --out diversity.csv
3. Compute Divergence (D)
Use compute_divergence.py to calculate divergence from an outgroup:
bash
python scripts/compute_divergence.py --alignment <alignment.fa> --outgroup macaque --out divergence.csv
4. Estimate Nₑ Ratios
Use Ne_ratio_calculator.py to compute ratios:
bash
python scripts/Ne_ratio_calculator.py --diversity diversity.csv --divergence divergence.csv --out Ne_ratios.csv
Example Data
See data/diversity_example.csv and data/divergence_example.csv for expected input formats.
Output
The final output (Ne_ratios_output.csv) includes:
Species
NₑX:NₑA
NₑY:NₑA
NₑMT:NₑA

# Citation
Lu, L., Dai, W., Pan, Y., & Yan, Z. (in preparation). A Divergence-Calibrated Framework for Robust Cross-Species Comparison of Effective Population Size Ratios. Manuscript in preparation for Methods in Ecology and Evolution.

(The manuscript and associated tool are available at: https://github.com/Ne-ratio/Ne_ratio)

# Contact
For questions, please contact Zheng Yan (zheng.yan@lzu.edu.cn).
