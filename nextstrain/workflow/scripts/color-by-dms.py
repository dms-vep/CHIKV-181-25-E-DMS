#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import AlignIO 
from augur.utils import write_json
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Make an Auspice color scheme for traits of interest.")
    parser.add_argument(
        "--data",
        required=True,
        help="Path table with data to plot on the tree."
    )
    parser.add_argument(
        "--alignment",
        required=True,
        help="Path to the final protein-level alignment."
    )
    parser.add_argument(
        "--reference",
        required=True,
        type=str,
        help="Name of the reference strain."
    )
    parser.add_argument(
        "--column",
        required=True,
        type=str,
        help="Name of the data column to color the tree by."
    )
    parser.add_argument(
        "--data-offset",
        required=True,
        type=int,
        help="Skip 'data-offset' positions in the data that aren't in the reference."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output JSON for augur."
    )
    return parser.parse_args()


def main():
    """Main entry point of the script."""
    args = parse_arguments()

    # Read in the data file
    print(f"Reading in data file {args.data} and making tree data for {args.column}:\n")
    data_df = pd.read_csv(args.data)

    # Replace the spaces in all column names with underscores
    data_df.columns = data_df.columns.str.replace(" ", "_")

    # Read in the alignment file
    alignment_dict = dict()
    for sequence in AlignIO.read(args.alignment, "fasta"):
        alignment_dict[str(sequence.id)] = str(sequence.seq)

    # Get the sequence of the reference (library) strain
    assert args.reference in alignment_dict, f"{args.reference} not in alignment_dict"
    reference_seq = alignment_dict[args.reference]

    # Get the mutations for each sequence relative to the reference
    mutant_dict = dict()
    for id, seq in alignment_dict.items():
        # Check that this is a valid alignment
        if len(seq) != len(reference_seq):
            raise ValueError(f"Length of {seq} does not match length of reference sequence {reference_seq}, alignment is corrupted")
        # Find mutations relative to the reference
        mutations = [(ref, i, alt) for i, (ref, alt) in enumerate(zip(reference_seq, seq)) if ref != alt]
        mutant_dict[id] = mutations

    # Find the index in the reference sequence where the data starts
    data_seq = "".join(data_df[data_df['mutant'] == data_df['wildtype']].reset_index(drop=True).wildtype)
    assert data_seq[args.data_offset:] in reference_seq, "Reference sequence does not contain the data sequence."
    reference_start = reference_seq.find(data_seq[args.data_offset:])

    # Update the mutant_dict to match the positions in the data
    data_mutant_dict = dict()
    for id, mutations in mutant_dict.items():
        data_mutants = []
        for ref, pos, alt in mutations:
            if pos >= reference_start and alt != "-":
                data_mutants.append((ref, (((pos - reference_start) + 1) + args.data_offset), alt))
        data_mutant_dict[id] = data_mutants

    # Score each sequence based on its mutations
    lookup_dict = {}
    for _, row in data_df.iterrows():
        key = (row['wildtype'], row['sequential_site'], row['mutant'])
        lookup_dict[key] = row[args.column]
    score_dict = {}
    for id, mutations in data_mutant_dict.items():
        score = sum(lookup_dict.get((ref, pos, alt), 0) for ref, pos, alt in mutations)
        score_dict[id] = score

    # Calculate the min and max scores
    print(f"Min score: {min(score_dict.values())}")
    print(f"Max score: {max(score_dict.values())}")
    
    # Output a JSON with the node data
    output_dict = {
            "nodes": defaultdict(dict)
        }
    for id, score in score_dict.items():
        output_dict["nodes"][id] = {
            args.column: score
        }
    write_json(output_dict, args.output)


if __name__ == "__main__":
    main()