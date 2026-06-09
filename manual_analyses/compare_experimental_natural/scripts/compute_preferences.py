"""Convert DMS functional scores for one phenotype into amino-acid preferences.

For site r and amino acid a, the functional score f_{r,a} is 0 for the wildtype,
the measured score for measured non-wildtype mutations, and (for missing/NaN
non-wildtype mutations) the mean of the other measured non-wildtype scores at that
site. Preferences are pi_{r,a} = exp(f_{r,a}) / sum_a' exp(f_{r,a'}). Run via the
Snakemake `script:` directive.
"""

import sys

import numpy as np
import pandas as pd

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log


def dms_wildtype_by_site(df, n_sites):
    wt = df[["sequential_site", "wildtype"]].drop_duplicates()
    nuniq = wt.groupby("sequential_site").wildtype.nunique()
    bad = list(nuniq[nuniq != 1].index)
    if bad:
        raise ValueError(f"sites with inconsistent wildtype residue: {bad}")
    s = wt.set_index("sequential_site").wildtype.sort_index()
    if list(s.index) != list(range(1, n_sites + 1)):
        raise ValueError(f"sequential_site values are not the range 1..{n_sites}")
    return s


def region_for_site(seqsite, region_bounds):
    for region, (lo, hi) in region_bounds.items():
        if lo <= seqsite <= hi:
            return region
    raise ValueError(f"sequential_site {seqsite} is not in any region")


column = snakemake.params.column
amino_acids = list(snakemake.params.amino_acids)
n_sites = snakemake.params.n_sites
region_bounds = snakemake.params.region_bounds
slug = snakemake.wildcards.phenotype

dms = pd.read_csv(snakemake.input.dms)
if column not in dms.columns:
    raise ValueError(f"phenotype column {column!r} not in DMS CSV")
wt_by_site = dms_wildtype_by_site(dms, n_sites)
aa_set = set(amino_acids)

rows = []
for site, wildtype in wt_by_site.items():
    if wildtype not in aa_set:
        raise ValueError(f"site {site} wildtype {wildtype!r} is not a standard aa")
    sub = dms[dms.sequential_site == site]

    measured = {}
    for mutant, score in zip(sub.mutant, sub[column]):
        if mutant == wildtype:
            continue
        if mutant not in aa_set:
            raise ValueError(f"site {site} has non-standard mutant {mutant!r}")
        if pd.notna(score):
            if mutant in measured:
                raise ValueError(f"site {site} mutant {mutant} appears twice")
            measured[mutant] = float(score)

    if not measured:
        raise ValueError(
            f"site {site} has no measured non-wildtype {slug} score; cannot impute"
        )
    mean_nonwt = float(np.mean(list(measured.values())))

    f = {}
    for aa in amino_acids:
        if aa == wildtype:
            f[aa] = 0.0
        elif aa in measured:
            f[aa] = measured[aa]
        else:
            f[aa] = mean_nonwt  # imputed
    n_imputed = 19 - len(measured)  # 19 non-wildtype amino acids

    exp = {aa: np.exp(f[aa]) for aa in amino_acids}
    denom = sum(exp.values())
    prefs = {aa: exp[aa] / denom for aa in amino_acids}
    if not np.isclose(sum(prefs.values()), 1.0):
        raise ValueError(f"site {site} preferences do not sum to 1")

    row = {
        "sequential_site": site,
        "region": region_for_site(site, region_bounds),
        "wildtype": wildtype,
        "n_measured_mutations": len(measured),
        "n_imputed_mutations": n_imputed,
    }
    row.update(prefs)
    rows.append(row)

out = pd.DataFrame(
    rows,
    columns=[
        "sequential_site",
        "region",
        "wildtype",
        "n_measured_mutations",
        "n_imputed_mutations",
        *amino_acids,
    ],
)
out.to_csv(snakemake.output[0], index=False)
print(
    f"{slug}: wrote preferences for {len(out)} sites; "
    f"imputed {int(out.n_imputed_mutations.sum())} mutation scores total"
)
