# Ne Ratio Calculator

This repository contains a computational framework for estimating relative effective population size (Ne) ratios across different genomic compartments using data from VCF files. The method uses interspecies divergence to calibrate for mutation rate variation.

## Introduction

The framework estimates Ne ratios between autosomes, X chromosome, Y chromosome, and mitochondrial DNA (mtDNA). The primary advantage of the method is its use of interspecies divergence to correct for varying mutation rates, enabling robust and reliable estimation of Ne ratios across species.

## Dependencies

To run the scripts, you need the following dependencies:

- Python 3.7+
- NumPy
- SciPy
- Biopython

You can install them using the following command:
pip install -r requirements.txt

## Scripts

- `calculate_ne.py`: This script takes VCF files as input, computes nucleotide diversity (π) across sliding windows, calculates sequence divergence, and then estimates Ne ratios.
- `popgenWindows.py`: This script calculates nucleotide diversity using non-overlapping 50 kb sliding windows.

## Usage

1. Place your VCF file(s) in the `/data` folder or provide the path to your file.
2. Run the `calculate_ne.py` script as follows:
python calculate_ne.py --vcf /path/to/your_file.vcf
This will output the Ne ratios for the various genomic compartments.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
