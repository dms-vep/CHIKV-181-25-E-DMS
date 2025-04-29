#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create codon alignment based on protein alignment.")
    parser.add_argument(
        "--alignment", 
        required=True,
        help="Path to an alignment file of the protein sequences"
    )
    parser.add_argument(
        "--sequences", 
        required=True,
        help="Path to a fasta file of the untranslated sequences"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to codon alignment"
    )
    return parser.parse_args()


def codon_alignment(alignment, sequences):
    """
    Create a codon alignment based on protein alignment.
    Use the untranslated sequences to get the codons.

    Parameters
    ----------
    alignment: str
        Path to an alignment file of the protein sequences
    sequences: str
        Path to a fasta file of the untranslated sequences

    Returns
    -------
    list
        A list of aligned codon seq records
    """
    codons = []
    for accession, translation in alignment.items():
        sequence = sequences[accession]
        codon_sequence = ""
        codon_index = 0
        for amino_acid in translation:
            if amino_acid == "-":
                codon_sequence += "---"
            else:
                codon_sequence += sequence[codon_index:codon_index + 3]
                codon_index += 3
        codons.append(codon_sequence)
    return codons


def main():
    """Main entry point of the script."""
    args = parse_arguments()
    print(f"Creating codon alignment based on protein alignment {args.alignment} and untranslated sequences {args.sequences}\n")

    # Read in the fasta files
    alignment = SeqIO.to_dict(SeqIO.parse(args.alignment, "fasta"))
    sequences = SeqIO.to_dict(SeqIO.parse(args.sequences, "fasta"))
    assert (sequences.keys() == alignment.keys()), "Keys in sequences and alignment dictionaries don't match"
    
    # Create the codon alignment
    with open(args.output, "w") as handle:
        SeqIO.write(codon_alignment(alignment, sequences), handle, "fasta")

    print(f"\nFinished. Results exported to {args.output}")


if __name__ == "__main__":

    main()