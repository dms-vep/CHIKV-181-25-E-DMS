"""Write amino-acid preferences in phydms format, numbered by codon-alignment column.

phydms expects a preferences file whose "site" column matches the alignment column
numbering (1..L). The E codon alignment column c corresponds to sequential site c+1.
Run via the Snakemake `script:` directive.
"""

import sys

import numpy as np
import pandas as pd

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log

amino_acids = list(snakemake.params.amino_acids)

prefs = pd.read_csv(snakemake.input.prefs).set_index("sequential_site")
codon_map = pd.read_csv(snakemake.input.codon_map)

rows = []
for _, m in codon_map.iterrows():
    seqsite = int(m["sequential_site"])
    if seqsite not in prefs.index:
        raise ValueError(f"no preferences for sequential_site {seqsite}")
    pref_row = prefs.loc[seqsite, amino_acids].astype(float)
    if not np.isclose(pref_row.sum(), 1.0):
        raise ValueError(f"preferences for site {seqsite} do not sum to 1")
    row = {"site": int(m["codon_column"])}
    row.update({aa: pref_row[aa] for aa in amino_acids})
    rows.append(row)

out = pd.DataFrame(rows, columns=["site", *amino_acids])
if len(out) != len(codon_map):
    raise ValueError("number of preference rows != number of codon columns")
out.to_csv(snakemake.output[0], index=False)
print(f"wrote phydms preferences for {len(out)} sites")
