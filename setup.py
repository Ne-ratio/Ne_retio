python
from setuptools import setup, find_packages
import os

# Read version from package
with open(os.path.join("ne_calculator", "__version__.py")) as f:
    exec(f.read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="ne_ratio_calculator",
    version=__version__,  # noqa: F821
    author="Liangyu Lu, Wennan Dai, Yuhao Pan, Zheng Yan",
    author_email="yanz@lzu.edu.cn",
    description="A framework for estimating effective population size ratios across genomic compartments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Ne-ratio/Ne_retio",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "ne-calculate=scripts.run_ne_analysis:main",
            "ne-process-vcf=scripts.process_vcf:main",
        ],
    },
    include_package_data=True,
)
