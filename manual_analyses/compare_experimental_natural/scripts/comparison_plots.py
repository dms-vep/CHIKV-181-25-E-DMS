"""Build comparison figures and the "more conserved than expected" table.

For one entry phenotype, plots three site-aligned tracks over the protein: natural
alignment entropy vs preference entropy (independent y-scales), relative solvent
accessibility, and the phydms ExpCM per-site significance as a signed -log10(P) (negative
= more conserved in nature than expected, positive = more variable). A linked hover
highlights the same site across all three panels. Also writes a per-site summary table
(entropy, RSA, and ExpCM/YNGKP_M0 omega and P for every site). Run via the Snakemake
`script:` directive.

The x-axis uses the original DMS `site` label (e.g. "-1(E3)") sorted by `sequential_site`,
and every point carries a tooltip with both the site label and the sequential site.
"""

import sys

import altair as alt
import numpy as np
import pandas as pd

snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log

alt.data_transformers.disable_max_rows()

amino_acids = list(snakemake.params.amino_acids)
slug = snakemake.wildcards.phenotype

C_ALN = "#1f77b4"  # alignment entropy (left axis)
C_PREF = "#d62728"  # preference entropy (right axis)
C_RSA = "#2ca02c"
NUM_FORMAT = ".3g"  # at most 3 significant digits everywhere numbers are shown
SLP_CAP = (
    5  # clamp the *plotted* signed -log10 P magnitude so it doesn't blow out the scale
)
# colors for the signed-significance direction
DIR_CONSERVED = "more conserved in nature"
DIR_VARIABLE = "more variable in nature"
DIR_NEUTRAL = "neutral (ω=1)"


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


def x_site():
    """Shared x encoding: DMS site label, ordered by sequential_site, no label overlap."""
    return alt.X(
        "site:N",
        sort=alt.SortField("sequential_site", order="ascending"),
        axis=alt.Axis(labelAngle=-90, labelOverlap="greedy"),
        title="site",
    )


def yscale():
    """Quantitative y-scale that just spans the data with a little padding."""
    return alt.Scale(nice=False, padding=10)


def tt_q(field, title=None):
    """Quantitative tooltip field shown with at most 3 significant digits."""
    kw = {"title": title} if title else {}
    return alt.Tooltip(field, type="quantitative", format=NUM_FORMAT, **kw)


# hover that links the same site across every panel: mousing over a point sets the shared
# `site` value, enlarging and outlining the matching point in all panels at once. (No
# `nearest`: it compiles to a voronoi whose tooltip drops the selection's `site` field.)
hover = alt.selection_point(
    fields=["site"], on="mouseover", empty=False, clear="mouseout"
)


def hover_enc(base=30, large=250):
    """size/stroke encodings that grow + outline the hovered site (shared selection)."""
    return dict(
        size=alt.condition(hover, alt.value(large), alt.value(base)),
        stroke=alt.condition(hover, alt.value("black"), alt.value(None)),
        strokeWidth=alt.condition(hover, alt.value(2.5), alt.value(0)),
    )


prefs = pd.read_csv(snakemake.input.prefs)
aln_ent = pd.read_csv(snakemake.input.aln_entropy)
rsa = pd.read_csv(snakemake.input.rsa)
# the DMS `site` labels (e.g. "-1(E3)") live only in the DMS summary csv
dms = pd.read_csv(snakemake.input.dms)
site_map = dms[["sequential_site", "site"]].drop_duplicates()
expcm = read_omegabysite(snakemake.input.omega_expcm).rename(
    columns={"omega": "omega_ExpCM", "P": "P_ExpCM", "dLnL": "dLnL_ExpCM"}
)
m0 = read_omegabysite(snakemake.input.omega_m0).rename(
    columns={"omega": "omega_M0", "P": "P_M0", "dLnL": "dLnL_M0"}
)

pref_ent = preference_entropy(prefs).rename("preference_entropy")
data = prefs[["sequential_site", "region", "wildtype"]].copy()
data = data.merge(site_map, on="sequential_site", how="left")
if data.site.isna().any():
    missing = data.loc[data.site.isna(), "sequential_site"].tolist()
    raise ValueError(f"no DMS site label for sequential sites {missing}")
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

# signed phydms significance: magnitude |log10(P)|, signed by whether omega is >1 (more
# variable in nature than expected) or <1 (more conserved); omega == 1 -> 0.
with np.errstate(divide="ignore", invalid="ignore"):
    data["signed_logP_ExpCM"] = -np.log10(data.P_ExpCM) * np.sign(data.omega_ExpCM - 1)
data["direction"] = np.select(
    [data.omega_ExpCM > 1, data.omega_ExpCM < 1],
    [DIR_VARIABLE, DIR_CONSERVED],
    default=DIR_NEUTRAL,
)
data.loc[data.omega_ExpCM.isna(), "direction"] = None
# plotted value is clamped to +/- SLP_CAP; the tooltip still shows the real magnitude
data["signed_logP_ExpCM_plot"] = data["signed_logP_ExpCM"].clip(-SLP_CAP, SLP_CAP)

# --- panel 1: entropy, points only, independent y-scale per series ---
ent_tooltip = [
    "site",
    tt_q("sequential_site"),
    tt_q("alignment_entropy"),
    tt_q("preference_entropy"),
]
aln_y = alt.Y(
    "alignment_entropy:Q",
    title="alignment entropy (nats)",
    scale=yscale(),
    axis=alt.Axis(titleColor=C_ALN, format=NUM_FORMAT),
)
pref_y = alt.Y(
    "preference_entropy:Q",
    title="preference entropy (nats)",
    scale=yscale(),
    axis=alt.Axis(titleColor=C_PREF, orient="right", format=NUM_FORMAT),
)
entropy_panel = (
    alt.layer(
        alt.Chart(data)
        .mark_point(filled=True, color=C_ALN, opacity=0.7)
        .encode(x=x_site(), y=aln_y, tooltip=ent_tooltip, **hover_enc())
        .add_params(hover),
        alt.Chart(data)
        .mark_point(filled=True, color=C_PREF, opacity=0.7)
        .encode(x=x_site(), y=pref_y, tooltip=ent_tooltip, **hover_enc())
        .add_params(hover),
    )
    .resolve_scale(y="independent")
    .properties(
        width=900,
        height=200,
        title=f"entropy: natural alignment vs {slug} preferences",
    )
)

# --- panel 2: relative solvent accessibility (same site domain; gaps where missing) ---
rsa_panel = (
    alt.Chart(data)
    .mark_point(filled=True, color=C_RSA, opacity=0.7)
    .encode(
        x=x_site(),
        y=alt.Y(
            "rsa:Q",
            title="relative solvent accessibility",
            scale=yscale(),
            axis=alt.Axis(format=NUM_FORMAT),
        ),
        tooltip=["site", tt_q("sequential_site"), tt_q("rsa")],
        **hover_enc(),
    )
    .add_params(hover)
    .properties(width=900, height=140, title="relative solvent accessibility (E2, E1)")
)

# --- panel 3: ExpCM signed -log10(P) ---
SLP_H = 180
zero_rule = alt.Chart(pd.DataFrame({"y": [0]})).mark_rule(color="gray").encode(y="y:Q")
slp_points = (
    alt.Chart(data)
    .mark_point(filled=True, opacity=0.8)
    .encode(
        x=x_site(),
        y=alt.Y(
            "signed_logP_ExpCM_plot:Q",
            title=f"signed −log10 P (ExpCM, capped ±{SLP_CAP})",
            scale=yscale(),
            axis=alt.Axis(format=NUM_FORMAT),
        ),
        color=alt.Color(
            "direction:N",
            title="natural selection",
            scale=alt.Scale(
                domain=[DIR_CONSERVED, DIR_VARIABLE, DIR_NEUTRAL],
                range=["#2166ac", "#b2182b", "#999999"],
            ),
        ),
        tooltip=[
            "site",
            tt_q("sequential_site"),
            tt_q("omega_ExpCM"),
            tt_q("P_ExpCM"),
            tt_q("signed_logP_ExpCM"),
        ],
        **hover_enc(),
    )
    .add_params(hover)
)
note_df = pd.DataFrame({"_": [0]})
note_top = (
    alt.Chart(note_df)
    .mark_text(align="left", baseline="top", dx=3, dy=3, fontSize=10, color="#555")
    .encode(
        x=alt.value(0),
        y=alt.value(0),
        text=alt.value("↑ more variable in nature than expected"),
    )
)
note_bottom = (
    alt.Chart(note_df)
    .mark_text(align="left", baseline="bottom", dx=3, dy=-3, fontSize=10, color="#555")
    .encode(
        x=alt.value(0),
        y=alt.value(SLP_H),
        text=alt.value("↓ more conserved in nature than expected"),
    )
)
signed_logP_panel = alt.layer(zero_rule, slp_points, note_top, note_bottom).properties(
    width=900, height=SLP_H, title="phydms ExpCM signed significance"
)

chart = alt.vconcat(entropy_panel, rsa_panel, signed_logP_panel).resolve_scale(
    color="independent"
)
chart.save(snakemake.output.html)

# --- per-site summary table (one row per site; region is the protein, e.g. E2) ---
summary = data.sort_values("sequential_site")[
    [
        "site",
        "sequential_site",
        "region",
        "alignment_entropy",
        "preference_entropy",
        "rsa",
        "omega_ExpCM",
        "P_ExpCM",
        "omega_M0",
        "P_M0",
    ]
]
summary.to_csv(snakemake.output.summary, index=False, float_format="%.3g")
print(f"{slug}: wrote figure and per-site summary for {len(summary)} sites")
