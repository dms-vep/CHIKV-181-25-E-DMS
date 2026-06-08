"""Build comparison figures and the "more conserved than expected" table.

For one entry phenotype, compares natural-alignment entropy with preference entropy
over the length of the protein, shows relative solvent accessibility, shows the
phydms per-site omega (ExpCM and YNGKP_M0), and flags sites that are more conserved
in nature than the measured functional constraint predicts. Run via the Snakemake
`script:` directive.
"""

import sys

import altair as alt
import numpy as np
import pandas as pd

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log

alt.data_transformers.disable_max_rows()

amino_acids = list(snakemake.params.amino_acids)
region_bounds = snakemake.params.region_bounds
delta_threshold = snakemake.params.delta
slug = snakemake.wildcards.phenotype


def preference_entropy(prefs):
    """Per-site Shannon entropy (nats) of the preferences."""
    p = prefs[amino_acids].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(p > 0, -p * np.log(p), 0.0)
    return pd.Series(terms.sum(axis=1), index=prefs.sequential_site.values)


def read_omegabysite(path):
    """Parse a phydms *_omegabysite.txt file into site/omega/P/dLnL columns."""
    df = pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
    df.columns = [c.strip() for c in df.columns]
    needed = {"site", "omega", "P", "dLnL"}
    if not needed.issubset(df.columns):
        raise ValueError(f"{path} missing columns {needed - set(df.columns)}")
    # phydms 'site' is the codon-alignment column; sequential site is column + 1
    df["sequential_site"] = df["site"].astype(int) + 1
    return df[["sequential_site", "omega", "P", "dLnL"]]


def region_bands():
    return pd.DataFrame(
        [
            {"region": r, "start": lo - 0.5, "end": hi + 0.5}
            for r, (lo, hi) in region_bounds.items()
        ]
    )


prefs = pd.read_csv(snakemake.input.prefs)
aln_ent = pd.read_csv(snakemake.input.aln_entropy)
rsa = pd.read_csv(snakemake.input.rsa)
expcm = read_omegabysite(snakemake.input.omega_expcm).rename(
    columns={"omega": "omega_ExpCM", "P": "P_ExpCM", "dLnL": "dLnL_ExpCM"}
)
m0 = read_omegabysite(snakemake.input.omega_m0).rename(
    columns={"omega": "omega_M0", "P": "P_M0", "dLnL": "dLnL_M0"}
)

pref_ent = preference_entropy(prefs).rename("preference_entropy")
data = prefs[["sequential_site", "region", "wildtype"]].copy()
data = data.merge(
    pref_ent.rename_axis("sequential_site").reset_index(), on="sequential_site"
)
data = data.merge(
    aln_ent[["sequential_site", "alignment_entropy"]],
    on="sequential_site",
    how="left",
)
data = data.merge(rsa[["sequential_site", "rsa"]], on="sequential_site", how="left")
data = data.merge(expcm, on="sequential_site", how="left")
data = data.merge(m0, on="sequential_site", how="left")
data["entropy_delta"] = data.preference_entropy - data.alignment_entropy

# --- entropy over the protein, with region bands ---
long = data.melt(
    id_vars=["sequential_site", "region"],
    value_vars=["alignment_entropy", "preference_entropy"],
    var_name="entropy_type",
    value_name="entropy",
).dropna(subset=["entropy"])
bands = (
    alt.Chart(region_bands())
    .mark_rect(opacity=0.12)
    .encode(
        x="start:Q",
        x2="end:Q",
        color=alt.Color("region:N", title="region"),
    )
)
entropy_lines = (
    alt.Chart(long)
    .mark_line(opacity=0.8)
    .encode(
        x=alt.X("sequential_site:Q", title="sequential site"),
        y=alt.Y("entropy:Q", title="Shannon entropy (nats)"),
        color=alt.Color("entropy_type:N", title=""),
    )
)
entropy_panel = (
    (bands + entropy_lines)
    .properties(
        width=900,
        height=200,
        title=f"entropy: natural alignment vs {slug} preferences",
    )
    .resolve_scale(color="independent")
)

rsa_panel = (
    bands
    + alt.Chart(data.dropna(subset=["rsa"]))
    .mark_point(size=12, filled=True, opacity=0.6)
    .encode(
        x=alt.X("sequential_site:Q", title="sequential site"),
        y=alt.Y("rsa:Q", title="relative solvent accessibility"),
        color=alt.Color("region:N", title="region"),
    )
).properties(width=900, height=140, title="relative solvent accessibility (E2, E1)")

omega_long = data.melt(
    id_vars=["sequential_site", "region"],
    value_vars=["omega_ExpCM", "omega_M0"],
    var_name="model",
    value_name="omega",
).dropna(subset=["omega"])
omega_panel = (
    bands
    + alt.Chart(omega_long)
    .mark_point(size=12, filled=True, opacity=0.6)
    .encode(
        x=alt.X("sequential_site:Q", title="sequential site"),
        y=alt.Y("omega:Q", title="phydms omega"),
        color=alt.Color("model:N", title="model"),
    )
).properties(width=900, height=180, title="phydms per-site omega")

# --- scatter: alignment vs preference entropy, with flagged sites ---
scatter_data = data.dropna(subset=["alignment_entropy", "preference_entropy"]).copy()
scatter_data["more_conserved_than_expected"] = (
    scatter_data.entropy_delta >= delta_threshold
)
diag = pd.DataFrame({"x": [0, float(np.log(20))]})
scatter = (
    alt.Chart(scatter_data)
    .mark_circle(size=30, opacity=0.6)
    .encode(
        x=alt.X("alignment_entropy:Q", title="natural alignment entropy (nats)"),
        y=alt.Y("preference_entropy:Q", title="preference entropy (nats)"),
        color=alt.Color("more_conserved_than_expected:N", title="flagged"),
        tooltip=[
            "sequential_site",
            "region",
            "wildtype",
            "alignment_entropy",
            "preference_entropy",
            "entropy_delta",
            "rsa",
            "omega_ExpCM",
        ],
    )
    .properties(width=400, height=400, title="more conserved in nature than expected")
)
diag_line = (
    alt.Chart(diag).mark_line(color="black", strokeDash=[4, 4]).encode(x="x:Q", y="x:Q")
)
scatter_panel = diag_line + scatter

chart = alt.vconcat(entropy_panel, rsa_panel, omega_panel, scatter_panel).resolve_scale(
    color="independent"
)
chart.save(snakemake.output.html)

flagged = scatter_data[scatter_data.more_conserved_than_expected].sort_values(
    "entropy_delta", ascending=False
)
flagged = flagged[
    [
        "sequential_site",
        "region",
        "wildtype",
        "alignment_entropy",
        "preference_entropy",
        "entropy_delta",
        "rsa",
        "omega_ExpCM",
        "P_ExpCM",
        "omega_M0",
    ]
]
flagged.to_csv(snakemake.output.flagged, index=False)
print(
    f"{slug}: wrote figure and {len(flagged)} flagged sites "
    f"(entropy_delta >= {delta_threshold})"
)
