#!/usr/bin/env python3
import argparse
from Bio import AlignIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Remove gaps from an alignment relative to a reference.")
    parser.add_argument(
        "--alignment", 
        required=True,
        help="Path to an alignment in FASTA format."
    )
    parser.add_argument(
        "--reference", 
        required=True,
        help="Strain name of the reference in the alignment."
    )
    parser.add_argument(
        "--output", 
        required=True,
        help="Path to save the ungapped alignment in FASTA format."
    )
    return parser.parse_args()


def ungap_alignment(alignment, reference):
    """
    Remove gaps from the alignment relative to the reference sequence.

    Parameters:
    -----------
    alignment : Bio.Align.MultipleSeqAlignment
        The alignment containing sequences with gaps.
    reference : str
        The strain name of the reference sequence in the alignment.

    Returns:
    --------
    Bio.Align.MultipleSeqAlignment
        The ungapped alignment with gaps removed relative to the reference.
    """
    # Get the reference sequence from the alignment
    reference_sequence = None
    for seq in alignment:
        if reference in seq.description or reference == seq.id:
            reference_sequence = str(seq.seq)
            break
    if reference_sequence is None:
        raise ValueError(f"Reference strain '{reference}' not found in the alignment.")

    # Remove gaps from the alignment relative to the reference
    num_gaps_removed = 0
    for index, char in enumerate(reference_sequence):
        if char == "-":
            alignment = alignment[:, :(index - num_gaps_removed)] + alignment[:, (index - num_gaps_removed + 1):]
            num_gaps_removed += 1

    print(f"Removed {num_gaps_removed} gaps from the alignment relative to the reference '{reference}'.\n")

    return alignment


def main():
    """
    Based on code from Caleb Carr
    """
    args = parse_arguments()
    print(f"Removing gaps from alignment {args.alignment} relative to the reference '{args.reference}'.\n")
    # Load in the alignment file
    alignment = AlignIO.read(args.alignment, "fasta")
    if not any(record.id == args.reference for record in alignment):
        raise ValueError(f"Reference strain '{args.reference}' not found in the alignment.")
    # Remove gaps from the alignment relative to the reference
    ungapped_alignment = ungap_alignment(alignment, args.reference)
    # Save the ungapped alignment to the output file
    AlignIO.write(ungapped_alignment, args.output, "fasta")
    print(f"Wrote new alignment to {args.output}.")

if __name__ == "__main__":
    main()