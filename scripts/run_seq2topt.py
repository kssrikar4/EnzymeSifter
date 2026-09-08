#!/usr/bin/env python3

import argparse
import math
import sys
import os
import contextlib
from pathlib import Path
import pandas as pd
import torch
import warnings
from Bio import SeqIO

warnings.filterwarnings("ignore",
    message="Setting attributes on ParameterList is not supported.")

def load_model(seq2topt_dir, model_type="topt", device="cpu"):
    code_dir = os.path.join(seq2topt_dir, "code")
    if code_dir not in sys.path:
        sys.path.insert(0, code_dir)

    from model import MultiAttModel

    dim = 320
    window = 3
    n_head = 4
    n_RD = 4

    model = MultiAttModel(dim, window, n_head, n_RD).to(device)

    weights_dir = os.path.join(seq2topt_dir, "weights")
    if model_type == "topt":
        weight_file = os.path.join(weights_dir, "model_topt_window.3_r2.0.57.pth")
        scale_factor = 120.0
    elif model_type == "tm":
        weight_file = os.path.join(weights_dir, "model_tm_window.3_r2.0.76.pth")
        scale_factor = 100.0
    else:
        raise ValueError(f"Unknown model_type: {model_type}")

    if not os.path.isfile(weight_file):
        sys.exit(f"[ERROR] Weight file not found: {weight_file}")

    state_dict = torch.load(weight_file, map_location=device, weights_only=False)
    model.load_state_dict(state_dict)
    model.eval()

    return model, scale_factor

def get_esm_model(device="cpu"):
    import esm
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    return model, batch_converter

def predict_sequences(fasta_path, seq2topt_dir, do_topt=False, do_tm=False):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        df_cols = ["seq_id"]
        if do_topt: df_cols.append("predicted_topt_C")
        if do_tm: df_cols.append("predicted_tm_C")
        return pd.DataFrame(columns=df_cols)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[INFO] Using device: {device}", file=sys.stderr)

    print(f"[INFO] Loading ESM-2 model...", file=sys.stderr)
    esm_model, batch_converter = get_esm_model(device=device)

    model_topt, scale_topt = None, None
    model_tm, scale_tm = None, None

    if do_topt:
        print(f"[INFO] Loading Seq2Topt model...", file=sys.stderr)
        model_topt, scale_topt = load_model(seq2topt_dir, "topt", device=device)
    if do_tm:
        print(f"[INFO] Loading Seq2Tm model...", file=sys.stderr)
        model_tm, scale_tm = load_model(seq2topt_dir, "tm", device=device)

    seq_ids = [rec.id for rec in records]
    seq_strs = [str(rec.seq) for rec in records]

    preds_topt = []
    preds_tm = []
    
    batch_size = 64 if device.type == "cuda" else 4

    for i in range(math.ceil(len(seq_ids) / batch_size)):
        batch_ids = seq_ids[i * batch_size: (i + 1) * batch_size]
        batch_seqs = seq_strs[i * batch_size: (i + 1) * batch_size]
        inputs = [(batch_ids[j], batch_seqs[j]) for j in range(len(batch_ids))]

        try:
            _, _, batch_tokens = batch_converter(inputs)
            batch_tokens = batch_tokens.to(device)

            with torch.no_grad():
                autocast_ctx = torch.autocast(device_type=device.type, dtype=torch.float16) if device.type == "cuda" else contextlib.nullcontext()
                with autocast_ctx:
                    results_esm = esm_model(batch_tokens, repr_layers=[6], return_contacts=False)
                    emb = results_esm["representations"][6].transpose(1, 2)
                    
                    if do_topt:
                        ptopt = model_topt(emb).cpu().numpy().reshape(-1).tolist()
                        preds_topt.extend(ptopt)
                    if do_tm:
                        ptm = model_tm(emb).cpu().numpy().reshape(-1).tolist()
                        preds_tm.extend(ptm)

        except Exception as e:
            print(f"[WARN] Batch {i} failed: {e}", file=sys.stderr)
            if do_topt: preds_topt.extend([float("nan")] * len(batch_ids))
            if do_tm: preds_tm.extend([float("nan")] * len(batch_ids))

    out_dict = {"seq_id": seq_ids}
    if do_topt:
        out_dict["predicted_topt_C"] = [v * scale_topt for v in preds_topt]
    if do_tm:
        out_dict["predicted_tm_C"] = [v * scale_tm for v in preds_tm]

    return pd.DataFrame(out_dict)

def main():
    p = argparse.ArgumentParser(
        description="Run Seq2Topt/Seq2Tm predictions on a FASTA file."
    )
    p.add_argument("--fasta", required=True,
                   help="Input multi-FASTA file")
    p.add_argument("--seq2topt_dir", required=True,
                   help="Path to cloned Seq2Topt repository")
    p.add_argument("--output_topt", default=None,
                   help="Output TSV for Topt predictions")
    p.add_argument("--output_tm", default=None,
                   help="Output TSV for Tm predictions")
    args = p.parse_args()

    if not args.output_topt and not args.output_tm:
        sys.exit("[ERROR] At least one of --output_topt or --output_tm required.")

    fasta = Path(args.fasta).resolve()
    seq2topt_dir = Path(args.seq2topt_dir).resolve()

    if not fasta.exists():
        sys.exit(f"[ERROR] FASTA file not found: {fasta}")

    df = predict_sequences(
        fasta, seq2topt_dir, 
        do_topt=bool(args.output_topt), 
        do_tm=bool(args.output_tm)
    )

    if args.output_topt:
        out = Path(args.output_topt)
        out.parent.mkdir(parents=True, exist_ok=True)
        df[["seq_id", "predicted_topt_C"]].to_csv(out, sep="\t", index=False)
        print(f"[INFO] Wrote {len(df)} Topt predictions to {out}", file=sys.stderr)

    if args.output_tm:
        out = Path(args.output_tm)
        out.parent.mkdir(parents=True, exist_ok=True)
        df[["seq_id", "predicted_tm_C"]].to_csv(out, sep="\t", index=False)
        print(f"[INFO] Wrote {len(df)} Tm predictions to {out}", file=sys.stderr)

if __name__ == "__main__":
    main()
