"""Extract the E-region codon alignment (sequential sites 2..N) from codon_ungapped.

Slices the E region out of the structural-polyprotein codon alignment, verifies the
181/25 reference translates to the DMS wildtype, and writes the codon alignment plus
a column -> sequential_site map. Run via the Snakemake `script:` directive.
"""

import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log


def read_fasta(path):
    recs = {}
    for record in SeqIO.parse(path, "fasta"):
        if record.id in recs:
            raise ValueError(f"duplicate sequence id {record.id} in {path}")
        recs[record.id] = str(record.seq)
    return recs


def dms_wildtype_string(df, n_sites):
    wt = df[["sequential_site", "wildtype"]].drop_duplicates()
    nuniq = wt.groupby("sequential_site").wildtype.nunique()
    bad = list(nuniq[nuniq != 1].index)
    if bad:
        raise ValueError(f"sites with inconsistent wildtype residue: {bad}")
    s = wt.set_index("sequential_site").wildtype.sort_index()
    if list(s.index) != list(range(1, n_sites + 1)):
        raise ValueError(f"sequential_site values are not the range 1..{n_sites}")
    return "".join(s.values)


def region_for_site(seqsite, region_bounds):
    for region, (lo, hi) in region_bounds.items():
        if lo <= seqsite <= hi:
            return region
    raise ValueError(f"sequential_site {seqsite} is not in any region")


ref_strain = snakemake.params.ref_strain
capsid_len = snakemake.params.capsid_len
n_sites = snakemake.params.n_sites
region_bounds = snakemake.params.region_bounds

dms = pd.read_csv(snakemake.input.dms)
e_wt = dms_wildtype_string(dms, n_sites)[1:]  # sequential sites 2..N
n_e = len(e_wt)

codon = read_fasta(snakemake.input.codon)
start = capsid_len * 3
end = (capsid_len + n_e) * 3

if ref_strain not in codon:
    raise ValueError(f"reference {ref_strain} not in {snakemake.input.codon}")
ref_e = codon[ref_strain][start:end]
if str(Seq(ref_e).translate()) != e_wt:
    raise ValueError("reference codon translation does not match DMS wildtype")

# phydms rejects ambiguous (e.g. "GNN") and partial-gap codons, so only sequences
# whose entire E region is fully resolved ACGT (no ambiguous nucleotide, no gap) are
# kept. The reference is always fully resolved (verified above); fail fast otherwise.
resolved = set("ACGT")
n_written = n_dropped = 0
written_names = set()
with open(snakemake.output.fasta, "w") as f:
    for name, seq in codon.items():
        if len(seq) < end:
            raise ValueError(f"sequence {name} shorter than expected E region")
        e = seq[start:end]
        if set(e) <= resolved:
            f.write(f">{name}\n{e}\n")
            written_names.add(name)
            n_written += 1
        else:
            n_dropped += 1
if ref_strain not in written_names:
    raise ValueError(
        f"reference {ref_strain} unexpectedly dropped as not fully resolved"
    )

rows = [
    {
        "codon_column": i + 1,  # 1-based within the E codon alignment
        "sequential_site": i + 2,  # E region begins at sequential site 2
        "region": region_for_site(i + 2, region_bounds),
    }
    for i in range(n_e)
]
pd.DataFrame(rows).to_csv(snakemake.output.mapping, index=False)
print(
    f"wrote E codon alignment: kept {n_written} fully-resolved seqs, "
    f"dropped {n_dropped} with ambiguous/gap codons ({n_e} codons)"
)
