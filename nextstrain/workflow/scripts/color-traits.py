#!/usr/bin/env python3
import json
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Make an Auspice color scheme for traits of interest.")
    parser.add_argument(
        "--metadata",
        required=True,
        help="Path containing traits for each sequence."
    )
    parser.add_argument(
        "--schemes", 
        required=True,
        help="Path to a file containing color schemes."
    )
    parser.add_argument(
        "--traits", 
        nargs="+",
        required=True,
        help="List of metadata columns to treat as traits (space-separated)"
    )
    parser.add_argument(
        "--ordering",
        type=json.loads,
        default=None,
        help="Group key columns by value columns before sorting"
    )
    parser.add_argument(
        "--output", 
        required=True,
        help="Path to output the color scheme file."
    )
    return parser.parse_args()


def main():
    """Main entry point of the script."""
    args = parse_arguments()

    # Read the metadata file
    metadata_df = pd.read_csv(args.metadata, sep="\t")

    # Read the color schemes file
    schemes = {}
    with open(args.schemes) as f:
        for row, line in enumerate(f.readlines()):
            array = line.lstrip().rstrip().split("\t")
            schemes[row] = array

    # Assign colors to the various traits
    trait_mappings = {}
    for trait in args.traits:
        if trait not in metadata_df.columns:
            raise ValueError(f"Trait '{trait}' not found in metadata columns: {metadata_df.columns.tolist()}")
        
        print(f"\nMaking color scheme for trait '{trait}':\n--------------------------------------")
        
        # Sort geographical traits by their hierarchy
        if trait in args.ordering:
            print(f"Ordering '{trait}' by {args.ordering[trait]}.")
            if args.ordering[trait] not in metadata_df.columns:
                raise ValueError(f"Ordering column '{args.ordering[trait]}' not found in metadata columns: {metadata_df.columns.tolist()}")
            groups = list(metadata_df[args.ordering[trait]].unique())
            trait_values = []
            for group in groups:
                values = sorted([t for t in metadata_df[metadata_df[args.ordering[trait]] == group][trait].unique() if t != "?"])
                trait_values.extend(values)
        else:
            trait_values = sorted([t for t in metadata_df[trait].unique() if t != "?"])
        
        # Print the unique trait values
        for value in trait_values:
            print(f"\t- {value}")

        # Map each trait value to a color from the scheme
        trait_scheme = schemes[len(trait_values)]
        color_mapping = {value: color for value, color in zip(trait_values, trait_scheme)}
        trait_mappings[trait] = color_mapping

    # Combine all trait mappings into a single dictionary
    output_rows = []
    for trait, mapping in trait_mappings.items():
        for value, color in mapping.items():
            output_rows.append([trait, value, color])
    output_df = pd.DataFrame(output_rows)
    output_df.to_csv(args.output, sep="\t", index=False, header=False)
    print(f"\nFinished. Results saved to {args.output}")


if __name__ == "__main__":
    main()