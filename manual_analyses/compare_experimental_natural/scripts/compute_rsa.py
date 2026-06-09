"""Compute relative solvent accessibility (RSA) of E2 and E1 sites from 6NK6.

6NK6 biological assembly 1 is the full icosahedral virion (~960 chains in one model):
~240 copies each of E2 and E1 plus the Mxra8 receptor chains. SASA is computed on the
whole assembly with the Mxra8 (and any other non-E) chains excluded, so burial between
E protomers and between capsomers counts but Mxra8 burial does not. RSA = residue SASA
/ Tien et al. (2013) max ASA, and the per-site value is the mean over all monomer copies
of that site in the virion. Run via the Snakemake `script:` directive.

freesasa uses single-character chain labels, which cannot distinguish the ~960
multi-character chain ids (e.g. "A-59") in this assembly. Because SASA depends only on
atom coordinates, every retained atom is added under one chain label with a globally
unique residue-number key, and that key is mapped back to its sequential site.
"""

import sys

import freesasa
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

try:  # biopython renamed/removed three_to_one across versions
    from Bio.PDB.Polypeptide import three_to_one as _three_to_one

    def aa3to1(resname):
        return _three_to_one(resname)

except ImportError:
    from Bio.PDB.Polypeptide import protein_letters_3to1

    def aa3to1(resname):
        return protein_letters_3to1[resname]


snakemake = snakemake  # noqa: F821  (injected by the Snakemake script directive)
sys.stdout = sys.stderr = open(snakemake.log[0], "w")  # capture all output to the log

# E region of the polyprotein corresponds to sequential sites 2..N (the E3 leading
# Met, site 1, is absent from the alignments and is not in the structure either).
FIRST_E_SITE = 2


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


def chain_residues(chain):
    """Return ordered (residue, one-letter aa) for standard amino acids in a chain."""
    out = []
    for res in chain:
        if not (is_aa(res, standard=True) and res.id[0] == " "):
            continue
        if res.id[2] != " ":
            raise ValueError(f"insertion code at {chain.id}{res.id} not supported")
        out.append((res, aa3to1(res.resname)))
    return out


def map_chain_to_sites(chain_seq, ref_full, aligner):
    """Local-align chain to the E reference; return (identity, coverage, idx->ref_idx)."""
    aln = aligner.align(ref_full, chain_seq)[0]
    ref_blocks, query_blocks = aln.aligned
    idx_to_ref = {}
    identical = aligned = 0
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        for k in range(r1 - r0):
            ref_i, q_i = r0 + k, q0 + k
            idx_to_ref[q_i] = ref_i
            aligned += 1
            if ref_full[ref_i] == chain_seq[q_i]:
                identical += 1
    identity = identical / aligned if aligned else 0.0
    coverage = aligned / len(chain_seq) if chain_seq else 0.0
    return identity, coverage, idx_to_ref


n_sites = snakemake.params.n_sites
region_bounds = snakemake.params.region_bounds
max_asa = snakemake.params.max_asa
min_identity = snakemake.params.min_identity
min_coverage = snakemake.params.min_coverage

dms = pd.read_csv(snakemake.input.dms)
wt_by_site = dms_wildtype_by_site(dms, n_sites)
# reference covering sequential sites 2..N (index 0 -> site 2)
ref_full = "".join(wt_by_site.values)[1:]

aligner = PairwiseAligner()
aligner.mode = "local"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

structure = MMCIFParser(QUIET=True).get_structure("6NK6", snakemake.input.assembly)
models = list(structure)
if len(models) != 1:
    raise ValueError(
        f"expected 1 model in the assembly, found {len(models)}; the biological "
        "assembly should contain all chains in one coordinate frame"
    )
model = models[0]

# Identify E chains by local alignment to the E reference; everything else (Mxra8,
# etc.) is dropped. For each retained chain, map its residues to sequential sites.
retained = {}  # chain_id -> {residue object: (seqsite, aa)}
n_kept = n_dropped = 0
kept_idents = []
for chain in model:
    residues = chain_residues(chain)
    if not residues:
        n_dropped += 1
        continue
    chain_seq = "".join(aa for _, aa in residues)
    identity, coverage, idx_to_ref = map_chain_to_sites(chain_seq, ref_full, aligner)
    if identity >= min_identity and coverage >= min_coverage:
        res_to_site = {
            res: (idx_to_ref[idx] + FIRST_E_SITE, aa)
            for idx, (res, aa) in enumerate(residues)
            if idx in idx_to_ref
        }
        retained[chain.id] = res_to_site
        n_kept += 1
        kept_idents.append(identity)
    else:
        n_dropped += 1

if not retained:
    raise ValueError("no chains matched the CHIKV E reference; cannot compute RSA")
print(
    f"kept {n_kept} E chains (identity {min(kept_idents):.2f}-{max(kept_idents):.2f}), "
    f"dropped {n_dropped} non-E chains (e.g. Mxra8)"
)

# Build a freesasa structure for the whole retained (Mxra8-free) assembly. Every atom
# uses one chain label with a globally unique residue-number key, so the ~480 E chains
# do not collide under freesasa's single-character chain model.
fs = freesasa.Structure()
key_to_site = {}  # unique residue key -> (seqsite, aa)
ridx = 0
for cid, res_to_site in retained.items():
    for res, (seqsite, aa) in res_to_site.items():
        rkey = str(ridx)
        for atom in res:
            if atom.element == "H":
                continue
            x, y, z = (float(c) for c in atom.get_coord())
            fs.addAtom(atom.get_fullname(), res.get_resname(), rkey, "A", x, y, z)
        key_to_site[rkey] = (seqsite, aa)
        ridx += 1

result = freesasa.calc(fs)
chain_areas = result.residueAreas()["A"]

records = []
for rkey, (seqsite, aa) in key_to_site.items():
    if region_for_site(seqsite, region_bounds) not in ("E2", "E1"):
        continue  # RSA reported for E2 and E1 only
    if rkey not in chain_areas:
        raise ValueError(f"freesasa missing area for residue key {rkey}")
    records.append(
        {"sequential_site": seqsite, "rsa": chain_areas[rkey].total / max_asa[aa]}
    )

df = pd.DataFrame(records)
summary = (
    df.groupby("sequential_site")
    .agg(rsa=("rsa", "mean"), n_chains=("rsa", "size"))
    .reset_index()
)
summary["region"] = summary.sequential_site.map(
    lambda s: region_for_site(s, region_bounds)
)
summary["wildtype"] = summary.sequential_site.map(wt_by_site)
summary = summary[["sequential_site", "region", "wildtype", "rsa", "n_chains"]]
summary = summary.sort_values("sequential_site")
summary.to_csv(snakemake.output.csv, index=False)
e2 = (summary.region == "E2").sum()
e1 = (summary.region == "E1").sum()
print(f"wrote RSA for {len(summary)} sites (E2={e2}, E1={e1})")
