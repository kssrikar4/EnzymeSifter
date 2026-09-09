# EnzymeSifter

A two-stage Snakemake pipeline for sifting through protein sequences and identifying promising enzyme candidates for downstream characterisation from raw FASTA all the way to ranked, clade-representative enzymes annotated with predicted biochemical values.

---

## Contents

- [Overview](#overview)
- [Pipeline architecture](#pipeline-architecture)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#Usage)
- [License](#license)

---

## Overview

The pipeline runs in two stages with a structure-prediction step (carried out externally by the user). Users are free to choose any of the optional filters:

1. **Stage 1 - filtering of sequences.** Filter by catalytic-residue motif, one or more Pfam families, and/or CLEAN-predicted EC number(s), then cluster at a user-defined identity threshold using MMseqs2.
2. **Structure prediction (ESMFold).** Generate PDB structures for Stage 1 filtered sequences.
3. **Stage 2 - structural screening.** Confirm enzymatic activity with EnzyMM, predict solubility/usability (NetSolP), optimal pH (pHoptNN), and optimal/melting temperatures (Seq2Topt), build an NJ tree, optionally partition it into clades, and select the best-scoring representative per clade according to your filtering criteria.

---

## Pipeline architecture

```
                          ┌────────────────────────────────┐
                          │           Input FASTA          │
                          └───────────────┬────────────────┘
                                          │
                  ╔═══════════════════════▼═══════════════════════╗
                  ║                    STAGE 1                    ║
                  ║  motif → Pfam (HMMER) → EC (CLEAN) → MMseqs2  ║
                  ╚═══════════════════════╤═══════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │ data/stage1/nonredundant.fasta │
                          └───────────────┬────────────────┘
                                          │
                  ╔═══════════════════════▼═══════════════════════╗
                  ║              predict 3D structures (ESMFold)  ║
                  ║                run_pdb_prediction.sh          ║
                  ╚═══════════════════════╤═══════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │     directory of .pdb files    │
                          └───────────────┬────────────────┘
                                          │
   ╔══════════════════════════════════════▼══════════════════════════════════════╗
   ║                                  STAGE 2                                    ║
   ║   EnzyMM ── filter to enzymatic hits ──┬── NetSolP   (solubility, usable)   ║
   ║                                        ├── pHoptNN   (pH optimum)           ║
   ║                                        ├── Seq2Topt  (T optimum)            ║
   ║                                        ├── Seq2Tm    (T melting)            ║
   ║                                        └── MUSCLE ─► NJ tree ─► clades ──►  ║
   ║                                                       representatives       ║
   ╚══════════════════════════════════════╤══════════════════════════════════════╝
                                          │
                          ┌───────────────▼────────────────┐
                          │  predictions_output/*.tsv      │
                          │  data/trees/nj_tree_clades.png │
                          └────────────────────────────────┘
```
---

## Requirements

- [Conda](https://docs.conda.io/en/latest/miniconda.html)

EnzymeSifter creates and manages its own conda environments automatically. You do not need to install any of the underlying tools yourself.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/Dar-Omar/EnzymeSifter.git
cd EnzymeSifter
```

Make the run scripts executable:

```bash
chmod +x run_stage1.sh run_stage2.sh run_pdb_prediction.sh
```

On the first invocation of either stage, a single unified conda environment (`enzymesifter`) will be created automatically, and the relevant setup scripts (`scripts/setup_*.sh`) will fetch external databases and model weights into `external/`.

---

## Usage


### Stage 1

```
./run_stage1.sh path/to/fasta [options]
```
#### Example
```bash
./run_stage1.sh ~/protein.fasta -residues GDSGGP -pfam PF00089 -identity 90
```

### Between the stages - structure prediction

Stage 2 needs a directory of PDB files of the filtered sequences. You can generate these locally using the provided ESMFold script:

```bash
./run_pdb_prediction.sh data/stage1/nonredundant.fasta data/predicted_pdbs
```
*Note: This script automatically scales across hardware. It uses chunking and FP16 to run efficiently on low-VRAM GPUs (>=4GB), and will seamlessly fall back to CPU memory mapping if needed.*

### Stage 2

```
./run_stage2.sh /path/to/pdbs [options]
```
#### Example
```bash
./run_stage2.sh data/predicted_pdbs/ -solubility 0.69 -phopt 7:10 -topt 30:45 -tm 55 -clades 13
```

Users can use -threads n at either stage to specify the number of CPU cores to use. If not used, the tool default to all available cores ($(nproc)).
  See [tutorial](tutorial.md) for complete features of the pipeline.

---

## License

- **Source code**: MIT License (see `LICENSE`)
- **Predicted structures** (`/pdbs/`): subject to the
  [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms).
  See `/pdbs/TERMS.txt` for details. Non-commercial use only.
  



---

## Performance Optimizations
This version of EnzymeSifter has been heavily optimized for speed and throughput compared to the original repository:
- **GPU Acceleration**: Neural network components including `CLEAN`, `pHoptNN`, `Seq2Topt`, and `NetSolP` have been updated to utilize CUDA-enabled PyTorch and ONNX Runtime if a compatible GPU is present. Mixed precision (`torch.autocast`) and batched feature extraction are used for maximum throughput.
- **Batched Feature Extraction**: `Seq2Topt` now extracts ESM-2 embeddings once in memory and passes them through both Topt and Tm heads simultaneously, eliminating duplicate sequence processing.
- **Hardware Concurrency**: Sub-processes (e.g. `mmseqs2`, `hmmsearch`, and `muscle`) are now automatically supplied with explicit thread directives matching your `--threads` arguments for full utilization of multicore CPUs.
- **Fast Tree Building**: BioPython's pure-Python implementation of distance matrices and NJ-tree building has been replaced with `FastTree` for massive speedups on large alignments.

---

## Acknowledgements
This repository is an optimized fork originally developed by the Bashton-Lab team. All credit for the pipeline concept, the overall architecture, and the foundational scripts goes to the original authors. Please visit the original [Bashton-Lab/EnzymeSifter](https://github.com/Bashton-Lab/EnzymeSifter) repository for the source project.
