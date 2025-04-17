#!/usr/bin/env python3
import argparse
from Bio import SeqIO, Phylo

def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--tree", 
        required=True,
        help=""
    )
    parser.add_argument(
        "--protein",
        required=True,
        help=""
    )
    parser.add_argument(
        "--codon",
        required=True,
        help=""
    )
    parser.add_argument(
        "--outgroup",
        required=True,
        help=""
    )
    parser.add_argument(
        "--output_tree",
        required=True,
        help=""
    )
    parser.add_argument(
        "--output_protein",
        required=True,
        help=""
    )
    parser.add_argument(
        "--output_codon",
        required=True,
        help=""
    )
    return parser.parse_args()


def root_and_remove_outgroup(input_tree_file, output_tree_file, outgroup):
    """Function to root tree and remove outgroup"""

    # Initialize tree
    tree = Phylo.read(input_tree_file, "newick")

    # Check that outgroup is in tree
    n_init = len(tree.get_terminals())
    assert any(clade.name == outgroup for clade in tree.get_terminals())

    # Root tree
    tree.root_with_outgroup(outgroup)
    tree.root = tree.root.clades[0]

    # Remove outgroup and check removal
    n_final = len(tree.get_terminals())
    assert not any(clade.name == outgroup for clade in tree.get_terminals())
    assert n_final == n_init - 1

    # Write new tree
    _ = Phylo.write(tree, output_tree_file, "newick")


def remove_sequence_from_alignment(input_alignment, output_alignment, outgroup):
    """Function to remove a sequence from an alignment"""

    # Initialize alignment and check sequence is in alignment
    alignment = list(SeqIO.parse(input_alignment, "fasta"))
    n_seqs_init = len(alignment)
    assert any(s.id == outgroup for s in alignment)

    # Remove sequence
    alignment = [s for s in alignment if s.id != outgroup]
    assert n_seqs_init == len(alignment) + 1

    # Write new alignment
    _ = SeqIO.write(alignment, output_alignment, "fasta")


def main():
    """Main entry point of the script."""
    # Parse arguments
    args = parse_arguments()

    # Root tree without group and remove outgroup
    root_and_remove_outgroup(args.tree, args.output_tree, args.outgroup)

    # Remove outgroup from alignments
    remove_sequence_from_alignment(args.protein, args.output_protein, args.outgroup)
    remove_sequence_from_alignment(args.codon, args.output_codon, args.outgroup)


if __name__ == "__main__":
    main()