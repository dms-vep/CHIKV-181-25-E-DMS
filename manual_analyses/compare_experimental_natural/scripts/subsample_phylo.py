"""Subsample sequences broadly across the phylogeny by farthest-point selection.

Greedily selects tips that maximize the minimum patristic distance to the already
selected set (max-min), which spreads the sample across the tree and down-weights
over-represented clades. Deterministic: the reference strain is the seed, and ties
are broken by taxon label. Run via the Snakemake `script:` directive.
"""

import sys

import dendropy
from Bio import SeqIO

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log


def read_fasta(path):
    recs = {}
    for record in SeqIO.parse(path, "fasta"):
        if record.id in recs:
            raise ValueError(f"duplicate sequence id {record.id} in {path}")
        recs[record.id] = str(record.seq)
    return recs


def farthest_point(pdm, taxa, seed_taxon, n):
    """Greedy max-min selection of n taxa starting from seed_taxon."""
    selected = [seed_taxon]
    selected_set = {seed_taxon}
    min_dist = {t: pdm.patristic_distance(seed_taxon, t) for t in taxa}
    while len(selected) < n:
        best = max(
            (t for t in taxa if t not in selected_set),
            key=lambda t: (min_dist[t], t.label),
        )
        selected.append(best)
        selected_set.add(best)
        for t in taxa:
            d = pdm.patristic_distance(best, t)
            if d < min_dist[t]:
                min_dist[t] = d
    return selected


n = snakemake.params.n
ref_strain = snakemake.params.ref_strain

seqs = read_fasta(snakemake.input.fasta)
tree = dendropy.Tree.get(
    path=snakemake.input.tree, schema="newick", preserve_underscores=True
)
taxa = list(tree.taxon_namespace)
if n > len(taxa):
    raise ValueError(f"requested {n} > {len(taxa)} available sequences")

seed = next((t for t in taxa if t.label == ref_strain), None)
if seed is None:
    raise ValueError(f"reference {ref_strain} not found in tree {snakemake.input.tree}")

pdm = tree.phylogenetic_distance_matrix()
selected = farthest_point(pdm, taxa, seed, n)
labels = [t.label for t in selected]
if ref_strain not in labels:
    raise ValueError("reference strain was not retained in the subsample")
if len(set(labels)) != n:
    raise ValueError("subsample contains duplicate taxa")

with open(snakemake.output.fasta, "w") as f:
    for label in labels:
        if label not in seqs:
            raise ValueError(f"tree taxon {label} not found in alignment")
        f.write(f">{label}\n{seqs[label]}\n")
with open(snakemake.output.ids, "w") as f:
    f.write("\n".join(labels) + "\n")
print(f"subsampled {len(labels)} of {len(taxa)} sequences (seed {ref_strain})")
