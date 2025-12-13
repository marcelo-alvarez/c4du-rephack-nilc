#!/usr/bin/env python
"""Run the needlet-ILC pipeline for CMB extraction.

This script runs the end-to-end needlet-ILC pipeline described in
arXiv:2307.01258 to extract CMB maps from ACT 90 GHz and 150 GHz data.

Usage:
    python scripts/run_nilc.py [options]

Options:
    --output-dir DIR    Output directory (default: ./output)
    --component STR     Stokes component: I, Q, or U (default: I)
    --target STR        Target component: cmb or ksz (default: cmb)
    --lmin INT          Minimum ell for high-pass filter (default: 100)
    --save-intermediate Save intermediate needlet maps
"""

import sys
import os
import argparse
from datetime import datetime

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from nilc.needlets.pipeline import run_nilc_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Run needlet-ILC pipeline for CMB extraction"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./output",
        help="Output directory for maps (default: ./output)"
    )
    parser.add_argument(
        "--component",
        type=str,
        default="I",
        choices=["I", "Q", "U"],
        help="Stokes component to process (default: I)"
    )
    parser.add_argument(
        "--target",
        type=str,
        default="cmb",
        choices=["cmb", "ksz"],
        help="Target component to extract (default: cmb)"
    )
    parser.add_argument(
        "--lmin",
        type=int,
        default=100,
        help="Minimum ell for high-pass filter (default: 100)"
    )
    parser.add_argument(
        "--save-intermediate",
        action="store_true",
        help="Save intermediate needlet maps"
    )

    args = parser.parse_args()

    print()
    print("=" * 70)
    print("NEEDLET-ILC CMB EXTRACTION")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    print()
    print(f"Configuration:")
    print(f"  Output directory:    {args.output_dir}")
    print(f"  Stokes component:    {args.component}")
    print(f"  Target component:    {args.target}")
    print(f"  lmin filter:         {args.lmin}")
    print(f"  Save intermediate:   {args.save_intermediate}")
    print()

    try:
        cmb_map = run_nilc_pipeline(
            output_dir=args.output_dir,
            component=args.component,
            lmin_cut=args.lmin,
            target_component=args.target,
            save_intermediate=args.save_intermediate,
        )

        print()
        print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print()
        return 0

    except Exception as e:
        print()
        print(f"ERROR: Pipeline failed with exception:")
        print(f"  {type(e).__name__}: {e}")
        print()
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
