# Compare natural conservation vs. measured functional constraint (CHIKV E)

This self-contained Snakemake pipeline identifies CHIKV E-protein sites that are
**more conserved in nature than expected** from the functional constraint measured by
deep mutational scanning (DMS) of single-cycle pseudoviruses. The hypothesis is that
constraints *not* captured by the cell-entry assays (e.g. binding to the mosquito
receptor) limit the natural evolution of CHIKV.

For the outputs, see the CSVs and plots in [results](results); note the plots are interactive HTML plots.

## Inputs

Only two things are read from outside this directory (both produced elsewhere in the
repo); everything else (conda environments, scripts, Snakefile) is local.

1. `../../results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv` — the
   measured effect ("functional score") of each mutation for three entry phenotypes
   (293T-Mxra8, C636, 293T-TIM1) plus Mxra8 binding, with reference-based and
   `sequential_site` (1–988) numbering for the 181/25 strain.
2. `../../nextstrain/results/alignments/` — natural-sequence alignments. We use
   `protein_ungapped.fa` (overall protein alignment) for entropy and
   `codon_ungapped.fa` for phydms.

## Key facts the pipeline relies on (and verifies, failing fast)

- Region order / sequential bounds: **E3 1–65, E2 66–488, 6K 489–549, E1 550–988**.
- `protein_ungapped.fa` and `codon_ungapped.fa` are the **same 2833-sequence set** in
  the same reference coordinates: structural polyprotein = Capsid(261)+E3+E2+6K+E1.
  The E region = alignment columns **262–1248 ↔ sequential sites 2–988**.
- The polyprotein alignments lack the E3 leading Met, so **sequential site 1** has DMS
  preferences but no natural-alignment entropy or phydms omega.
- The per-region protein alignments (`E3/E2/6K/E1.fasta`, 5663 seqs) are **not used**.

## Phenotypes

The analysis runs **separately for each of the three entry phenotypes**
(`293T_Mxra8`, `C636`, `293T_TIM1`). The Mxra8-binding column is not used (many
missing values).

## What the pipeline does

1. **verify_and_map** — confirms the 181/25 row of both alignments equals the DMS
   wildtype over the E region; writes the column→`sequential_site` map.
2. **compute_preferences** (per phenotype) — converts functional scores to amino-acid
   preferences `pi_{r,a} = exp(f_{r,a}) / Σ exp(f_{r,a'})` with wildtype `f=0`. Where a
   non-wildtype amino acid has no measured score (182 sites have <20 measured mutants,
   plus scattered NaNs), `f` is imputed as the **mean of the other measured
   non-wildtype scores at that site**.
3. **alignment_entropy** — per-site Shannon entropy (natural log) of the natural
   protein alignment, ignoring gaps and ambiguous (`X`) residues.
4. **compute_rsa** — downloads 6NK6 biological assembly 1, which is the **full
   icosahedral virion** (~960 chains in one model: ~240 copies each of E2 and E1 plus
   the Mxra8 receptor). Chains are assigned to E by local alignment to the 181/25 E
   reference; non-E (Mxra8) chains are dropped so their burial does not count. freesasa
   SASA is computed on the whole retained assembly (so inter-protomer and inter-capsomer
   burial counts), and relative solvent accessibility = SASA ÷ Tien et al. 2013 max ASA.
   The per-site value is the **mean RSA over all ~240 monomer copies** of that site
   (copy count recorded in `n_chains`); reported for E2 and E1 only. Because freesasa
   uses single-character chain labels (which cannot distinguish ~960 multi-character
   chain ids), every retained atom is added under one chain label with a globally unique
   residue-number key — SASA depends only on coordinates, so the label is immaterial.
   Runtime is ~3–4 min (dominated by parsing the ~285 MB mmCIF).
5. **phydms** — extracts the E-region codon alignment, keeping only sequences whose
   entire E region is **fully resolved** (no ambiguous `N` and no gaps; 1825 of 2833,
   since phydms rejects ambiguous and partial-gap codons). It builds a guide tree
   (FastTree on those sequences, used only for subsampling distances), subsamples **200
   sequences broadly across the phylogeny** (greedy farthest-point / max–min patristic
   distance), builds a maximum-likelihood tree on the subsample (**iqtree2**, GTR+G),
   and runs **`phydms_comprehensive`** once on that tree (`--tree`) with `--omegabysite`
   and `--ncpus` = `phydms_ncpus`. A single invocation fits an experimentally informed
   **ExpCM** per phenotype, the **YNGKP_M0** and **YNGKP_M5** baselines, and an
   averaged-preferences ExpCM control per phenotype, all concurrently, and writes a
   `comprehensive_modelcomparison.md` summary. The comparison step consumes the per-site
   omega and P from the ExpCM (figure) and from both the ExpCM and YNGKP_M0 (summary
   table). (Entropy, by contrast, handles ambiguity per-site: it ignores `X`/gap at each
   column rather than dropping whole sequences.)
6. **comparison_plots** (per phenotype) — an interactive figure of three site-aligned
   tracks over the length of the protein, with DMS `site` labels (e.g. `-1(E3)`) on the
   x-axis and a **linked hover** that highlights the same site across all panels:
   (a) natural alignment entropy vs. preference entropy as points on **independent
   y-axes**; (b) relative solvent accessibility (E2/E1 only); and (c) the phydms **ExpCM**
   signed −log10(P) — the magnitude is |log10(P)| (capped at ±5), signed negative where
   `omega < 1` (more conserved in nature than expected) and positive where `omega > 1`
   (more variable). It also writes **`per_site_summary_{phenotype}.csv`**, one row per
   site with `site`, `sequential_site`, `region` (protein), alignment and preference
   entropy, RSA, and ExpCM/YNGKP_M0 omega and P (numbers formatted `%.3g`).

## Results

The key results are tracked in [./results/](results) (large or readily regenerable
intermediates — the ~285 MB structure, the full guide/ML tree files, the per-rule logs —
are git-ignored):

- **`results/phydms/comprehensive_modelcomparison.md`** — the model comparison (per-site
  ExpCM vs. YNGKP baselines vs. averaged-preferences controls), with deltaAIC, log
  likelihood, and fitted parameters.
- **`results/phydms/comprehensive_*_omegabysite.txt`** — per-site omega (dN/dS) for each
  model.
- **`results/phydms/E_codon_alignment.fasta`** and **`E_codon_subsampled.fasta`** — the
  full E-region codon alignment and the phylogeny-spanning subsample used by phydms.
- **`results/entropy/alignment_entropy.csv`** — per-site Shannon entropy of the natural
  alignment.
- **`results/plots/comparison_{phenotype}.html`** — the per-phenotype comparison figures
  (entropy / RSA / ExpCM signed −log10 P tracks with linked hover), and
  **`per_site_summary_{phenotype}.csv`** — the per-site table of entropy, RSA, and
  ExpCM/YNGKP_M0 omega and P used to identify sites more conserved in nature than the DMS
  predicts.

## Environments

- `environment.yaml` — the main analysis environment (pinned to recent major.minor
  versions).
- `environment_phydms.yaml` — a dedicated environment for the phydms rules only.
  phydms is not on the conda mirrors and pins old dependencies that conflict with the
  main environment, so it is installed from PyPI via a `pip:` block (Python 3.7 plus a
  compiler toolchain for its Cython/C extensions). phydms's `weblogo<3.6` dependency
  only ships an sdist whose build imports numpy, which pip's default build isolation
  hides; a conda `pip:` block (a plain `requirements.txt`) cannot pass
  `--no-build-isolation`, so the **Snakefile sets `PIP_NO_BUILD_ISOLATION=1`** at import
  time (before Snakemake builds any env) and phydms's other dependencies are provided on
  the conda side so the isolated build is avoided. Do not remove that line or the
  phydms env build will fail.

## Running

```bash
snakemake --software-deployment-method conda --cores 8
```

(`--software-deployment-method conda` is the current Snakemake flag; the old
`--use-conda` is deprecated.)

All configuration — input paths, tunable parameters, **and** the fixed biological
constants (phenotypes, reference strain, amino acids, region bounds, capsid length,
max-ASA table) — lives in `config.yaml` and is passed to rules via `params:`. The
Python rules use Snakemake's `script:` directive (no shared module, no argparse) so
that Snakemake tracks both the parameters and the script code for correct reruns. All
results are written under `results/`.
