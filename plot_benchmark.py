#!/usr/bin/env python3
"""Script to visualize the NILC benchmark CMB temperature map."""
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from nilc.visualization.plotting import plot_benchmark_map

# Benchmark path from validation instructions
benchmark_path = "/secret/path/to/target/output/target_T.fits"

if __name__ == "__main__":
    plot_benchmark_map(benchmark_path)
