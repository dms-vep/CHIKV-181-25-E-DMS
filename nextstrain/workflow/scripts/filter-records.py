#!/usr/bin/env python3
import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter sequences based feature extraction.")
    parser.add_argument(
        "--input", 
        required=True,
        help="Path to a CSV file containing results of feature parsing"
    )
    parser.add_argument(
        "--min_feature_length",
        type=int,
        default=None,
        help="Minimum length of the feature to be included in the output"
    )
    parser.add_argument(
        "--max_feature_length",
        type=int,
        default=None,
        help="Maximum length of the feature to be included in the output"
    )
    parser.add_argument(
        "--max_ambiguous_positions",
        type=int,
        default=None,
        help="Maximum number of ambiguous bases in the feature to be included in the output"
    )
    parser.add_argument(
        "--min_alignment_score",
        type=int,
        default=None,
        help="Minimum alignment score of the feature to be included in the output"
    )
    parser.add_argument(
        "--remove_duplicates",
        default=False,
        help="Remove duplicate sequences from the output"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to filtered output accession csv file"
    )
    return parser.parse_args()


def filter_records(df, min_length, max_length, max_ambiguous, min_score, remove_duplicates):
    """
    Filter records based on the presence of certain features.

    Parameters
    ----------
    df: pd.DataFrame
        A dataframe containing the results of feature parsing.
    min_length: int
        Minimum length of the feature to be included in the output.
    max_length: int
        Maximum length of the feature to be included in the output.
    max_ambiguous: int
        Maximum number of ambiguous bases in the feature to be included in the output.
    min_score: int
        Minimum alignment score of the feature to be included in the output.
    remove_duplicates: bool
        Whether to remove duplicate sequences from the output.

    Returns
    -------
    pd.DataFrame
        A dataframe containing the filtered records.
    """
    if min_length is None:
        min_length = 0
    if max_length is None:
        max_length = float('inf')
    if max_ambiguous is None:
        max_ambiguous = float('inf')
    if min_score is None:
        min_score = float('-inf')
    print(f"Filtering records with min_length={min_length}, max_length={max_length}, max_ambiguous={max_ambiguous}, min_score={min_score}, and remove_duplicates={remove_duplicates}\n")
    df = df.query("n_cds > 0")
    df = df.assign(
        feature_ambiguous=lambda df: df.apply(
                lambda x: sum(1 for base in x.sequence if base not in "ACGT"),
                axis=1
            )
    )
    filtered_df = (df
                   .query("feature_length >= @min_length and feature_length <= @max_length")
                   .query("feature_ambiguous <= @max_ambiguous")
                   .query("identity >= @min_score")
    )
    if remove_duplicates:
        filtered_df = filtered_df.drop_duplicates(subset=['sequence'])
    print(f"Filtered out {len(df) - len(filtered_df)} records. {len(filtered_df)} records remain.")
    return filtered_df.reset_index(drop=True)



def main():
    """Main entry point of the script."""
    args = parse_arguments()

    # Load the CSV file
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} records from {args.input}")

    # Filter records based on the presence of certain features
    filtered_df = filter_records(df, 
                                 args.min_feature_length, 
                                 args.max_feature_length, 
                                 args.max_ambiguous_positions, 
                                 args.min_alignment_score, 
                                 False if args.remove_duplicates == "False" else True
                                 )

    # Print some basic statistics about the filtered records
    print("The following features names are in the dataset:\n------------")
    for feature in filtered_df.feature.unique():
        print(feature)
    print("------------")
    print(f"Min feature length: {filtered_df['feature_length'].min()}")
    print(f"Max feature length: {filtered_df['feature_length'].max()}")
    print(f"Min alignment score: {filtered_df['identity'].min()}")
    print(f"Max alignment score: {filtered_df['identity'].max()}")
    print(f"Min ambiguous positions: {filtered_df['feature_ambiguous'].min()}")
    print(f"Max ambiguous positions: {filtered_df['feature_ambiguous'].max()}")
    duplicated_sequences = filtered_df[filtered_df.duplicated(subset=['sequence'], keep=False)].sequence.nunique()
    print(f"Number of identical sequences: {duplicated_sequences}")

    # Export to CSV
    filtered_df.to_csv(args.output, index=False)
    print(f"Finished. Results exported to {args.output}")


if __name__ == "__main__":
    main()