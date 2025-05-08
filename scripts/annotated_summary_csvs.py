"""Make annotated summary CSVs."""

import sys

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

data = pd.read_csv(snakemake.input.data)
addtl_annotations = pd.read_csv(snakemake.input.addtl_annotations)

shared_cols = sorted(set(data.columns).intersection(addtl_annotations.columns))
assert shared_cols, f"{data.columns=}, {addtl_annotations.columns=}"

value_cols = [
    c
    for c in data.columns
    if c not in {"site", "wildtype", "mutant", "sequential_site", "region"}
]

mut = data.merge(addtl_annotations, on=shared_cols, how="outer", validate="m:1").sort_values("sequential_site")
mut.to_csv(snakemake.output.mut, index=False, float_format="%.4g")

site_mean = (
    mut
    [(~mut["mutant"].isin(["*", "-"])) & (mut["mutant"] != mut["wildtype"])]
    .groupby([c for c in mut.columns if c not in ["mutant"] + value_cols], as_index=False)
    .aggregate(**{c: pd.NamedAgg(c, "mean") for c in value_cols})
    .sort_values("sequential_site")
)
if ("protein_site" in site_mean.columns) and all(site_mean["protein_site"] == site_mean["protein_site"].astype(int)):
    site_mean["protein_site"] = site_mean["protein_site"].astype(int)
site_mean.to_csv(snakemake.output.site_mean, index=False, float_format="%.4f")
