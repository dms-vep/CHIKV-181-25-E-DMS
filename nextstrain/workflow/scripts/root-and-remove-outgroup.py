#!/usr/bin/env python3
import argparse
from Bio import SeqIO, Phylo


def parse_arguments():
    parser = argparse.ArgumentParser(description="Remove the outgroup from a tree and alignments.")
    parser.add_argument(
        "--tree", 
        required=True,
        help="Path to raw newick tree file."
    )
    parser.add_argument(
        "--protein",
        required=True,
        help="Path to a protein-level alignment file."
    )
    parser.add_argument(
        "--codon",
        required=True,
        help="Path to a codon-phased alignment file."
    )
    parser.add_argument(
        "--outgroup",
        required=True,
        type=str,
        help="Accession of the outgroup to be removed from the tree and alignments."
    )
    parser.add_argument(
        "--output_tree",
        required=True,
        help="A newick tree file with the outgroup removed and re-rooted."
    )
    parser.add_argument(
        "--output_protein",
        required=True,
        help="A protein-level alignment file with the outgroup removed."
    )
    parser.add_argument(
        "--output_codon",
        required=True,
        help="A codon-phased alignment file with the outgroup removed."
    )
    return parser.parse_args()


def root_and_remove_outgroup(tree, outgroup):
    """
    Root tree and remove outgroup from the newick tree file.

    Parameters
    ----------
    tree : str
        Path to the input newick tree file.
    outgroup : str
        Accession of the outgroup to be removed from the tree.

    Returns
    -------
    
        Writes a new tree file with the outgroup removed and re-rooted.
    """

    # Check that outgroup is in tree
    n_init = len(tree.get_terminals())
    if not any(leaf.name == outgroup for leaf in tree.get_terminals()):
        raise ValueError(f"Outgroup '{outgroup}' not found in tree.")


    # Root tree with outgroup
    tree.root_with_outgroup(outgroup)
    tree.root = tree.root.clades[0]

    # Ensure that the outgroup is removed
    n_final = len(tree.get_terminals())
    assert not any(clade.name == outgroup for clade in tree.get_terminals())
    assert n_final == n_init - 1

    return tree


def remove_sequence_from_alignment(alignment, outgroup):
    """
    Remove a sequence from and alignment file.

    Parameters
    ----------
    alignment : str
       Alignment file.
    outgroup : str
        Accession of the outgroup to be removed from the alignment.

    Returns
    -------
        Writes a new alignment file with the outgroup removed.
    """

    # Check sequence is in alignment
    n_seqs_init = len(alignment)
    if not any(s.id == outgroup for s in alignment):
        raise ValueError("Outgroup not found in alignment.")

    # Remove sequence
    alignment = [s for s in alignment if s.id != outgroup]
    assert n_seqs_init == len(alignment) + 1
    return alignment


def main():
    """
    Main entry point of the script.
    
    Based on code from Caleb Carr.
    """
    # Parse arguments
    args = parse_arguments()

    # Root tree and remove outgroup
    print(f"Rooting and removing outgroup ({args.outgroup}) from tree...")
    tree = Phylo.read(args.tree, "newick")
    rooted_tree = root_and_remove_outgroup(tree, args.outgroup)
    Phylo.write(rooted_tree, args.output_tree, "newick")
    print(f"Writing tree to file: {args.output_tree}\n")

    # Remove outgroup sequence from the alignments
    print(f"Removing outgroup ({args.outgroup}) from alignments...")
    protein_alignment = list(SeqIO.parse(args.protein, "fasta"))
    protein_alignment_no_outgroup = remove_sequence_from_alignment(protein_alignment, args.outgroup)
    SeqIO.write(protein_alignment_no_outgroup, args.output_protein, "fasta")
    print(f"Writing protein alignment to file: {args.output_protein}")
    codon_alignment = list(SeqIO.parse(args.codon, "fasta"))
    codon_alignment_no_outgroup = remove_sequence_from_alignment(codon_alignment, args.outgroup)
    SeqIO.write(codon_alignment_no_outgroup, args.output_codon, "fasta")
    print(f"Writing codon alignment to file: {args.output_codon}\n")

if __name__ == "__main__":
    main()