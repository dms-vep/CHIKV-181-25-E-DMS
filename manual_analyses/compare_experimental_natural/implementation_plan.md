# Implementation plan: compare natural conservation vs. measured functional constraint

## 1. Context and goal

We want to identify CHIKV E-protein sites that are **more conserved in nature than we
would expect** given the functional constraint measured by deep mutational scanning
(DMS) in single-cycle pseudoviruses. The scientific hypothesis is that constraints
*not* captured by the cell-entry assays (for example, binding to the mosquito
receptor) limit the natural evolution of CHIKV. We bring together two data sources:

- **Measured mutation effects ("functional scores")** for each mutation in three
  entry phenotypes plus a Mxra8-binding phenotype.
- **Natural-sequence alignments** of the CHIKV E proteins.

The deliverables are (a) per-site amino-acid **preferences** derived from the DMS,
(b) a comparison of **site entropy** from natural alignments vs. from the
preferences, (c) **relative solvent accessibility (RSA)** of E2/E1 sites, and
(d) a formal **phydms** per-site `omega` analysis using an experimentally informed
codon model. Outputs include figures over the length of the protein and tables that
flag sites more conserved in nature than expected.

## 2. Hard constraints (must follow)

- **Self-contained.** The ONLY files read from parent directories are the two named
  inputs (below). All conda environments, scripts, `Snakefile`, and `config.yaml`
  live in this directory. Do **not** read parent environments, other parent files,
  or git history.
- **Two permitted inputs** (read-only, from parent dirs):
  1. `../../results/summaries/entry_293T-Mxra8_C636_293T-TIM1_Mxra8-binding.csv`
  2. `../../nextstrain/results/alignments/` (the alignment FASTA files)
- **Fail fast** with clear error messages on any unexpected input — alignment or
  numbering mismatches, missing strain rows, gap-only columns, preferences not
  summing to 1, missing structure chains, codon-translation mismatches, etc. No
  silent handling; no `dict.get` with masking defaults.
- Follow the repo conventions in `../../CLAUDE.md`: every rule has `log:` →
  `results/logs/{rule}_{wildcards}.txt`, a `conda:` directive, and writes output to a
  `results/` subdirectory; keep code concise; define shared constants in one place;
  keep `README.md` updated; run `snakefmt`, `snakemake --lint`, `black`, `ruff` after
  changes.

## 3. Decisions confirmed with the user

| Topic | Decision |
|---|---|
| Phenotypes | Run the whole analysis **separately for each of the three entry phenotypes**: `entry in 293T_Mxra8 cells`, `entry in C636 cells`, `entry in 293T_TIM1 cells`. The `binding to mouse Mxra8` column is **not** used for preferences (many missing values). |
| Alignment for entropy | Use the **overall protein alignment** `protein_ungapped.fa` (2833 seqs) — the same sequence set as the codon alignment used by phydms. **The per-region protein alignments (`E3/E2/6K/E1.fasta`, 5663 seqs) are NOT used for anything.** |
| Missing mutation scores | Some sites lack a measured score for some non-wildtype amino acids (182 sites have <20 measured mutants; plus scattered NaNs). For any **missing non-wildtype** amino acid at a site, set its functional score to the **mean of all other non-wildtype measured scores at that site** (per phenotype). Wildtype stays at 0. |
| phydms sequence set | **Subsample to 250 sequences**, chosen to sample **broadly across the phylogeny** so over-represented clades are down-weighted (greedy farthest-point / max–min patristic distance on a guide tree). |
| 6NK6 | **Download** in a rule; compute RSA on the **biological assembly excluding Mxra8 chains**; SASA tool = `freesasa`. |

## 4. Facts established from the permitted inputs (and how they drive the design)

### 4.1 DMS CSV
- Columns: `site, wildtype, mutant, entry in 293T_Mxra8 cells, entry in C636 cells,
  entry in 293T_TIM1 cells, binding to mouse Mxra8, sequential_site, region`.
- 988 sequential sites, 20 amino acids. **Not** a complete grid: 19462 rows vs. a full
  988×20 = 19760, and 182 sites have <20 measured mutants (min 3). Per-phenotype NaNs
  also occur in present rows (≈8/35/10 non-wildtype NaNs for 293T-Mxra8 / C636 /
  293T-TIM1). ⇒ preference computation must impute missing non-wildtype scores
  (Section 8, rule 3).
- Region order / sequential bounds: **E3 = 1–65, E2 = 66–488, 6K = 489–549,
  E1 = 550–988**. `sequential_site` is the continuous 1–988 numbering for 181/25.
- E3 reference numbering starts at `-1(E3)` (the leading Met, seqsite 1); there is no
  `0(E3)`.

### 4.2 Overall protein alignment (used for entropy)
- `protein_ungapped.fa`: **2833 sequences**, **1248 columns**, reference-coordinate
  structural polyprotein = Capsid(261) + E3(64) + E2(423) + 6K(61) + E1(439).
- The reference strain `181-25_MW473668` is present. Verified: its E region (columns
  **262–1248**) equals the DMS wildtype for **seqsites 2–988**.
- So the alignment column → sequential_site map is a single contiguous mapping
  (col 262 → seqsite 2, …, col 1248 → seqsite 988). **Seqsite 1** (the E3 leading Met)
  has no alignment column and therefore no natural-alignment entropy (it still has DMS
  preferences). This is the same coordinate system as the codon alignment, so entropy
  and phydms use the **same 2833-sequence set** and the **same site mapping**.
- The gapped `protein.fa` (1260 cols) and the per-region alignments are **not used**.

### 4.3 Codon alignment (used for phydms)
- `codon_ungapped.fa`: **2833 sequences**, **1248 codons**, identical coordinate system
  and sequence set as `protein_ungapped.fa`. E region = **codon columns 262–1248 →
  seqsites 2–988**. The gapped `codon.fa` is not used.

## 5. Directory layout to create

```
compare_experimental_natural/
├── Snakefile
├── config.yaml                 # input paths, phenotypes, n_subsample, seed, pdb id, thresholds
├── environment.yaml            # main analysis env
├── environment_phydms.yaml     # phydms-only env
├── README.md
├── implementation_plan.md      # this file
├── data/                       # downloaded 6NK6 biological assembly
├── scripts/
│   ├── shared.py               # single source of truth for constants/helpers
│   ├── verify_and_map.py       # verify ref alignment == DMS wt; build col→seqsite map
│   ├── compute_preferences.py
│   ├── alignment_entropy.py
│   ├── compute_rsa.py
│   ├── extract_E_codon_alignment.py
│   ├── subsample_phylo.py
│   ├── prefs_to_phydms.py
│   └── comparison_plots.py
└── results/
    ├── logs/
    ├── mapping/                # alncol → seqsite map + verification summary
    ├── preferences/            # prefs_{phenotype}.csv
    ├── entropy/                # alignment_entropy.csv, merged entropy tables
    ├── rsa/                    # site_rsa.csv
    ├── phydms/                 # E codon alns, trees, omegabysite_{phenotype}.csv
    └── plots/                  # figures + "more_conserved_than_expected" tables
```

### `scripts/shared.py` (centralized constants — defined once)
- `REGIONS = ["E3", "E2", "6K", "E1"]` and their sequential-site bounds.
- `ENTRY_PHENOTYPES` = the three entry column names; short slugs for filenames
  (`293T_Mxra8`, `C636`, `293T_TIM1`).
- `AMINO_ACIDS` (the 20).
- `REF_STRAIN = "181-25_MW473668"`.
- `CAPSID_LEN = 261`, so the E region = polyprotein columns 262–1248 ↔ seqsites 2–988.
- `MAX_ASA` = Tien et al. (2013) theoretical maxima (for RSA).
- Helpers: load DMS CSV; build the 1..988 DMS wildtype string; FASTA I/O.

## 6. Config (`config.yaml`)
Keys: input CSV path; alignments dir; the three entry phenotypes (column name → slug);
`n_subsample: 250`; `random_seed`; `pdb_id: 6NK6` + assembly URL; "more conserved than
expected" flagging thresholds; phydms options.

## 7. Conda environments
- **`environment.yaml`** (main): `python`, `snakemake-minimal`, `snakefmt`, `black`,
  `ruff`, `pandas`, `numpy`, `biopython`, `altair`, `vl-convert-python`, `freesasa`,
  `dendropy`, `fasttree`, `requests`.
- **`environment_phydms.yaml`** (phydms only): `phydms` from bioconda with its own
  pinned dependencies. **Rationale:** phydms pins older biopython/pandas/numpy that
  conflict with `altair` and modern `pandas`; a separate env avoids an unsolvable
  solve. Deliberate, documented deviation from the repo's usual single-env principle.
  (If a single env solves cleanly during implementation, we collapse to one.)

## 8. Snakemake rules

Every rule has `log:`, `conda:`, and writes to `results/`. Wildcard `{phenotype}`
ranges over the three entry-phenotype slugs.

1. **`download_6nk6`** → `data/6NK6-assembly1.cif`
   Download biological assembly 1 from RCSB (`requests`); assert non-empty/parseable.

2. **`verify_and_map`** (`verify_and_map.py`) → `results/mapping/alncol_to_seqsite.csv`,
   `results/mapping/verification.csv`
   Read the `REF_STRAIN` row from `protein_ungapped.fa`; assert it has 1248 columns and
   that columns 262–1248 (after the 261-residue capsid) equal the DMS wildtype for
   seqsites 2–988 (**fail fast** on any mismatch or unexpected length). Also assert the
   `codon_ungapped.fa` `REF_STRAIN` row translates to the same E-region sequence. Emit
   the single contiguous column→seqsite map (with region labels from the seqsite bounds)
   and a verification summary. (Records that seqsite 1 has no alignment column.)

3. **`compute_preferences`** (`compute_preferences.py`, per phenotype) →
   `results/preferences/prefs_{phenotype}.csv`
   Build the functional-score vector `f_{r,a}` for every site `r` and all 20 amino
   acids `a`:
   - wildtype `a`: `f = 0`;
   - non-wildtype with a measured (non-NaN) score: `f = score`;
   - **non-wildtype that is missing** (no row, or NaN for this phenotype):
     `f = mean of all other non-wildtype measured scores at site r` (this phenotype).
   Assert each site has ≥1 measured non-wildtype score (else cannot impute → fail fast).
   Then `pi_{r,a} = exp(f_{r,a}) / Σ_a' exp(f_{r,a'})`; assert all 20 aa present and
   prefs sum to 1 (within tolerance). Output: `sequential_site, region, wildtype,
   n_measured_mutations, <one column per amino acid>` (record how many were imputed).

4. **`alignment_entropy`** (`alignment_entropy.py`) → `results/entropy/alignment_entropy.csv`
   Over all 2833 sequences of `protein_ungapped.fa`, take the E-region columns
   (262–1248), compute per-column amino-acid frequencies ignoring gaps (assert ≥1
   non-gap residue per column), Shannon entropy (natural log). Map columns →
   sequential_site via the rule-2 map. Output: `sequential_site, region,
   alignment_entropy, n_nongap`. (Seqsite 1 absent by construction.)

5. **`compute_rsa`** (`compute_rsa.py`) → `results/rsa/site_rsa.csv`
   Parse the 6NK6 assembly (Biopython). Identify chains by sequence: keep chains whose
   sequence matches the E1 or E2 181/25 wildtype, **drop the Mxra8 chain(s)** (anything
   not E1/E2). **Fail fast** if E1 or E2 chains are not found. Run `freesasa` on the
   retained assembly (so inter-protomer burial counts but Mxra8 burial does not).
   RSA = residue SASA / `MAX_ASA[wildtype aa]`. Map structure residues →
   sequential_site for **E2 and E1 only**; assert structure residues agree with the DMS
   wildtype at mapped sites. Output: `sequential_site, region, wildtype, rsa`.

6. **`extract_E_codon_alignment`** (`extract_E_codon_alignment.py`) →
   `results/phydms/E_codon_alignment.fasta`, `results/phydms/codon_site_map.csv`
   Slice codon columns 262–1248 of `codon_ungapped.fa` (all 2833 seqs). Translate the
   `REF_STRAIN` row and **assert it equals the DMS wildtype for seqsites 2–988**
   (fail fast). Write the E-region codon alignment and the column→seqsite map.

7. **`guide_tree`** → `results/phydms/guide_tree.nwk`
   FastTree (GTR nt) on the full 2833-seq E codon alignment (guides subsampling only).

8. **`subsample_phylo`** (`subsample_phylo.py`) →
   `results/phydms/E_codon_subsampled.fasta`, `results/phydms/subsampled_ids.txt`
   Greedy farthest-point selection of 250 tips by patristic distance (dendropy): start
   from a deterministic seed tip, repeatedly add the tip maximizing the minimum distance
   to the selected set. Down-weights dense clades, spreads broadly across the phylogeny.
   Deterministic (fixed seed). Force-include `REF_STRAIN`; assert it is present.

9. **`phydms_tree`** → `results/phydms/phydms_tree.nwk`
   FastTree on the 250-seq alignment (phydms re-optimizes branch lengths under the codon
   model, so a FastTree starting topology suffices).

10. **`prefs_to_phydms`** (`prefs_to_phydms.py`, per phenotype) →
    `results/phydms/prefs_phydms_{phenotype}.csv`
    Convert preferences to phydms format, indexed 1..N over the codon-alignment columns
    (seqsites 2–988). Assert every alignment column has a preference row.

11. **`run_phydms`** (`environment_phydms.yaml`, per phenotype) →
    `results/phydms/omegabysite_{phenotype}.csv`
    Run `phydms` with the **ExpCM** (experimentally informed codon model) and
    `--omegabysite` on the 250-seq alignment + tree + preferences, plus a `YNGKP_M0
    --omegabysite` baseline. Parse per-site omega (and P-values) into a tidy CSV keyed
    by sequential_site.

12. **`comparison_plots`** (`comparison_plots.py`, per phenotype + combined) using
    **altair** → `results/plots/...`
    - **Preference entropy:** Shannon entropy of `pi_{r,·}` per site per phenotype.
    - **Track figure over protein length** (x = sequential_site; E3/E2/6K/E1 annotated):
      alignment entropy vs. preference entropy, plus an RSA track.
    - **Scatter** alignment entropy vs. preference entropy with the `y=x` diagonal; flag
      sites **more conserved in nature than expected** (low alignment entropy at high
      preference entropy) → `results/plots/more_conserved_than_expected_{phenotype}.csv`.
    - **phydms omega-by-site** over protein length (purifying sites, omega < 1,
      FDR-flagged).

13. **`all`** → aggregates final tables and figures for all phenotypes.

## 9. README.md
Describe purpose, the two inputs, the single 2833-seq set used for both entropy
(`protein_ungapped.fa`) and phydms (`codon_ungapped.fa`→250), that the per-region
alignments are unused, that seqsite 1 lacks alignment data, the preference imputation
rule, each rule's role, the two conda environments and why phydms is separate, and how
to run (`snakemake --use-conda --cores N`).

## 10. Verification after implementation
- `snakemake -n` builds the full DAG without error.
- `results/mapping/verification.csv` confirms the reference E region matches DMS
  wildtype for seqsites 2–988.
- Spot-check preferences sum to 1, imputed-count column is sensible, entropies ∈
  `[0, ln 20]`.
- Confirm `compute_rsa` dropped exactly the Mxra8 chain(s) and RSA ∈ ~`[0, 1.4]`.
- Confirm phydms produced omega-by-site for all sites for each phenotype.
- Run `snakefmt`, `snakemake --lint`, `black`, `ruff`; fix all findings.
- Sanity-review the comparison figures and the "more conserved than expected" tables.

## 11. Additional decisions confirmed with the user
- **phydms models:** run **both ExpCM and YNGKP_M0**, each with `--omegabysite`, so we
  contrast experimentally informed vs. naive per-site omega.
- **Subsampling:** use **greedy farthest-point (max–min patristic distance)** selection
  (deterministic, no extra dependencies).
- **Entropy units:** **natural log (nats)** for all Shannon entropies.
