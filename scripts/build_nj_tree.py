#!/usr/bin/env python3

import sys
import argparse
import subprocess
from pathlib import Path
from Bio import Phylo

def main():
    parser = argparse.ArgumentParser(
        description="Build a phylogenetic tree from an aligned FASTA using FastTree."
    )
    parser.add_argument("alignment", help="Input aligned FASTA (.afa)")
    parser.add_argument("output", help="Output Newick tree (.nwk)")
    parser.add_argument(
        "--model",
        default="blosum62",
        help="Distance model for DistanceCalculator (ignored, FastTree uses JTT by default)",
    )
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    print(f"[INFO] Running FastTree on {args.alignment}", file=sys.stderr)
    cmd = ["FastTree", "-out", args.output, args.alignment]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"[ERROR] FastTree failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    tree = Phylo.read(args.output, "newick")
    
    for clade in tree.find_clades():
        if clade.branch_length is not None and clade.branch_length < 0:
            clade.branch_length = 0.0

    tree.root_at_midpoint()
    tree.ladderize()

    for clade in tree.get_terminals():
        if clade.name and "|" in clade.name:
            clade.name = clade.name.split("|")[0]
            
    Phylo.write(tree, args.output, "newick")

    print(f"[INFO] Tree written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
