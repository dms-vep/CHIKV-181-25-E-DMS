"""Verify the reference alignment matches the DMS wildtype and build the column map.

Checks that the 181/25 row of the overall protein alignment (and the codon
alignment) equals the DMS wildtype over the E region, then writes the mapping from
alignment column -> sequential_site. Fails fast on any mismatch. Run via the
Snakemake `script:` directive.
"""

import sys

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log


def read_fasta(path):
    """Return dict id -> sequence; fails fast on duplicate ids."""
    recs = {}
    for record in SeqIO.parse(path, "fasta"):
        if record.id in recs:
            raise ValueError(f"duplicate sequence id {record.id} in {path}")
        recs[record.id] = str(record.seq)
    return recs


def dms_wildtype_string(df, n_sites):
    """Wildtype amino-acid string for contiguous sequential sites 1..n_sites."""
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


def report_mismatch(label, observed, expected):
    diffs = [
        (i + 2, e, o)  # +2: first E-region position is sequential site 2
        for i, (e, o) in enumerate(zip(expected, observed))
        if e != o
    ]
    raise ValueError(
        f"{label} does not match DMS wildtype; first mismatches "
        f"(seqsite, dms_wt, alignment) = {diffs[:10]} (total {len(diffs)})"
    )


ref_strain = snakemake.params.ref_strain
capsid_len = snakemake.params.capsid_len
n_sites = snakemake.params.n_sites
region_bounds = snakemake.params.region_bounds

dms = pd.read_csv(snakemake.input.dms)
full_wt = dms_wildtype_string(dms, n_sites)  # sites 1..N
# The polyprotein alignments lack the E3 leading Met (site 1), so the alignment E
# region corresponds to sequential sites 2..N.
e_wt = full_wt[1:]
n_e = len(e_wt)

# --- protein alignment ---
prot = read_fasta(snakemake.input.protein)
if ref_strain not in prot:
    raise ValueError(f"reference {ref_strain} not in {snakemake.input.protein}")
ref_prot = prot[ref_strain]
expected_len = capsid_len + n_e
if len(ref_prot) != expected_len:
    raise ValueError(f"protein alignment length {len(ref_prot)} != {expected_len}")
prot_e = ref_prot[capsid_len : capsid_len + n_e]
if prot_e != e_wt:
    report_mismatch("protein alignment", prot_e, e_wt)

# --- codon alignment ---
codon = read_fasta(snakemake.input.codon)
if ref_strain not in codon:
    raise ValueError(f"reference {ref_strain} not in {snakemake.input.codon}")
ref_codon = codon[ref_strain]
if len(ref_codon) != expected_len * 3:
    raise ValueError(f"codon alignment length {len(ref_codon)} != {expected_len * 3}")
codon_e_nt = ref_codon[capsid_len * 3 : (capsid_len + n_e) * 3]
codon_e_aa = str(Seq(codon_e_nt).translate())
if codon_e_aa != e_wt:
    report_mismatch("codon alignment (translated)", codon_e_aa, e_wt)

# --- build column -> sequential_site map (1-based columns) ---
mapping = pd.DataFrame(
    [
        {
            "alignment_column": capsid_len + 1 + i,
            "sequential_site": 2 + i,
            "region": region_for_site(2 + i, region_bounds),
        }
        for i in range(n_e)
    ]
)
mapping.to_csv(snakemake.output.mapping, index=False)

verification = pd.DataFrame(
    [
        {"check": "reference_strain", "value": ref_strain},
        {"check": "protein_alignment_length", "value": len(ref_prot)},
        {"check": "codon_alignment_length", "value": len(ref_codon)},
        {"check": "e_region_sites", "value": f"2..{n_sites}"},
        {"check": "n_sites_mapped", "value": len(mapping)},
        {"check": "protein_matches_dms_wildtype", "value": True},
        {"check": "codon_translation_matches_dms_wildtype", "value": True},
        {"check": "site_1_has_alignment_column", "value": False},
    ]
)
verification.to_csv(snakemake.output.verification, index=False)
print(f"verified reference matches DMS wildtype; mapped {len(mapping)} sites")
