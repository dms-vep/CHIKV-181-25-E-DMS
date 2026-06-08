"""Compute per-site Shannon entropy (nats) of the natural protein alignment.

Uses the overall protein alignment (protein_ungapped.fa). Gaps ('-') and ambiguous
residues ('X') are ignored; any other non-standard character fails fast. Run via the
Snakemake `script:` directive.
"""

import sys
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log

IGNORED = {"-", "X"}  # gaps and ambiguous residues


def read_fasta(path):
    recs = {}
    for record in SeqIO.parse(path, "fasta"):
        if record.id in recs:
            raise ValueError(f"duplicate sequence id {record.id} in {path}")
        recs[record.id] = str(record.seq)
    return recs


amino_acids = list(snakemake.params.amino_acids)
aa_set = set(amino_acids)

mapping = pd.read_csv(snakemake.input.mapping)
seqs = list(read_fasta(snakemake.input.protein).values())
lengths = {len(s) for s in seqs}
if len(lengths) != 1:
    raise ValueError(f"protein alignment has unequal sequence lengths: {lengths}")
aln_len = lengths.pop()

rows = []
for _, m in mapping.iterrows():
    col = int(m["alignment_column"])
    if not 1 <= col <= aln_len:
        raise ValueError(f"alignment column {col} out of range 1..{aln_len}")
    counts = Counter(s[col - 1] for s in seqs)
    unexpected = set(counts) - aa_set - IGNORED
    if unexpected:
        raise ValueError(f"column {col} has unexpected characters {unexpected}")
    aa_counts = np.array([counts.get(aa, 0) for aa in amino_acids], dtype=float)
    total = aa_counts.sum()
    if total < 1:
        raise ValueError(f"column {col} has no standard amino acids")
    p = aa_counts / total
    nz = p[p > 0]
    entropy = float(-(nz * np.log(nz)).sum())
    rows.append(
        {
            "sequential_site": int(m["sequential_site"]),
            "region": m["region"],
            "alignment_entropy": entropy,
            "n_nongap": int(total),
        }
    )

pd.DataFrame(rows).to_csv(snakemake.output[0], index=False)
print(f"wrote alignment entropy for {len(rows)} sites")
