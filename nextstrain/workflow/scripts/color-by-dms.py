#!/usr/bin/env python3
import argparse
import json
import os
import pandas as pd
from Bio import AlignIO 
from augur.utils import write_json
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Make an Auspice color scheme for traits of interest.")
    parser.add_argument(
        "--alignment",
        required=True,
        help="Path to the final protein-level alignment."
    )
    parser.add_argument(
        "--template-auspice",
        required=True,
        help="Path to the template auspice configuration file."
    )
    parser.add_argument(
        "--reference",
        required=True,
        type=str,
        help="Name of the reference strain."
    )
    parser.add_argument(
        "--dms-config",
        type=json.loads,
        required=True,
        help="Dictionary of information about DMS data to color the tree with."
    )
    parser.add_argument(
        "--auspice-config",
        required=True,
        help="Path to output JSON for augur."
    )
    parser.add_argument(
        "--dms-scores",
        nargs="+",
        required=True,
        help="Path to output JSON for augur."
    )
    return parser.parse_args()


def main():
    """Main entry point of the script."""
    args = parse_arguments()

    print(f"Determining mutations in the alignment file {args.alignment}:\n")
    # Determine the mutations in the alignment file
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

    # Read in the template auspice configuration file
    print(f"Reading in template auspice configuration file {args.template_auspice}:\n")
    with open(args.template_auspice) as f:
        auspice_config = json.load(f)

    # Make a JSON with the node data for each metric
    for output_json in args.dms_scores:
        metric = os.path.basename(output_json).split(".")[0]
        print(f"Creating JSON file ({output_json}) with the node data for {metric}:\n")
        
        # Read in the DMS data and replace the spaces in all column names with underscores
        data = args.dms_config[metric]["file"]
        offset = args.dms_config[metric]["offset"]
        print(f"Reading in data file {data}...")
        data_df = pd.read_csv(data)
        data_df.columns = data_df.columns.str.replace(" ", "_")

        # Find the index in the reference sequence where the data starts
        data_seq = "".join(data_df[data_df['mutant'] == data_df['wildtype']].reset_index(drop=True).wildtype)
        assert data_seq[offset:] in reference_seq, "Reference sequence does not contain the data sequence."
        reference_start = reference_seq.find(data_seq[offset:])

        # Update the mutant_dict to match the positions in the data
        data_mutant_dict = dict()
        for id, mutations in mutant_dict.items():
            data_mutants = []
            for ref, pos, alt in mutations:
                if pos >= reference_start and alt != "-":
                    data_mutants.append((ref, (((pos - reference_start) + 1) + offset), alt))
            data_mutant_dict[id] = data_mutants

        # Score each sequence based on its mutations
        lookup_dict = {}
        for _, row in data_df.iterrows():
            key = (row['wildtype'], row['sequential_site'], row['mutant'])
            lookup_dict[key] = row[metric]
        score_dict = {}
        for id, mutations in data_mutant_dict.items():
            score = sum(lookup_dict.get((ref, pos, alt), 0) for ref, pos, alt in mutations)
            score_dict[id] = score

        # Update the auspice configuration file for this metric
        data_config = {
            "key": metric,
            "title": args.dms_config[metric]["title"],
            "type": "continuous",
            "scale": [
                [min(score_dict.values()), args.dms_config[metric]['scale'][0]],
                [max(score_dict.values()), args.dms_config[metric]['scale'][1]]
            ]
        }
        auspice_config['colorings'].append(data_config)
        
        # Output a JSON with the node data
        output_dict = {
                "nodes": defaultdict(dict)
            }
        for id, score in score_dict.items():
            output_dict["nodes"][id] = {
                metric: score
            }
        write_json(output_dict, output_json)
        print(f"Output written to {output_json}.\n")

    # Write the template auspice configuration file
    print(f"Writing out the auspice configuration file {args.auspice_config}:\n")
    with open(args.auspice_config, 'w') as outfile:
        json.dump(auspice_config, outfile, indent=2)

if __name__ == "__main__":
    main()