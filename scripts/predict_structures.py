#!/usr/bin/env python3

import sys
import os
import torch
from pathlib import Path
from Bio import SeqIO
from transformers import EsmForProteinFolding

def main():
    if len(sys.argv) < 3:
        sys.exit(1)
        
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        sys.exit(0)
        
    print("[INFO] Loading ESMFold v1...", file=sys.stderr)
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    
    model.trunk.set_chunk_size(64)
    model.eval()
    
    if torch.cuda.is_available():
        vram = torch.cuda.get_device_properties(0).total_memory
        if vram >= 14 * 1024**3:
            model = model.cuda()
        elif vram >= 4 * 1024**3:
            model = model.half().cuda()
            print("[INFO] Using FP16 to fit model in limited VRAM.", file=sys.stderr)
    
    for rec in records:
        seq = str(rec.seq)
        if len(seq) > 1000:
            print(f"[WARN] Skipping {rec.id} (length {len(seq)} > 1000)", file=sys.stderr)
            continue
            
        out_path = os.path.join(output_dir, f"{rec.id}.pdb")
        if os.path.exists(out_path):
            continue
            
        print(f"[INFO] Predicting {rec.id}...", file=sys.stderr)
        with torch.no_grad():
            try:
                pdb_str = model.infer_pdb(seq)
                with open(out_path, "w") as f:
                    f.write(pdb_str)
            except Exception as e:
                print(f"[ERROR] Failed on {rec.id}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
