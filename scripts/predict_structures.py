#!/usr/bin/env python3

import sys
import os
import torch
from pathlib import Path
from Bio import SeqIO
from transformers import AutoTokenizer, EsmForProteinFolding, pipeline

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
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    
    model.trunk.set_chunk_size(64)
    model.eval()
    
    device_id = -1
    if torch.cuda.is_available():
        vram = torch.cuda.get_device_properties(0).total_memory
        if vram >= 14 * 1024**3:
            device_id = 0
            model = model.cuda()
        elif vram >= 4 * 1024**3:
            device_id = 0
            model = model.half().cuda()
            print("[INFO] Using FP16 to fit model in limited VRAM.", file=sys.stderr)
    
    fold_pipe = pipeline("protein-folding", model=model, tokenizer=tokenizer, device=device_id)
    
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
                res = fold_pipe(seq)
                if isinstance(res, list):
                    pdb_str = res[0]
                else:
                    pdb_str = res
                    
                if isinstance(pdb_str, dict) and "pdb" in pdb_str:
                    pdb_str = pdb_str["pdb"]
                    
                with open(out_path, "w") as f:
                    f.write(pdb_str)
            except Exception as e:
                print(f"[ERROR] Failed on {rec.id}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
