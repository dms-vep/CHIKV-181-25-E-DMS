#!/usr/bin/env python3
import argparse
import pickle
import pandas as pd
from Bio.SeqFeature import ExactPosition


def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter sequences based feature extraction.")
    parser.add_argument(
        "--features", 
        required=True,
        help="Path to a CSV file containing results of feature parsing"
    )
    parser.add_argument(
        "--records", 
        required=True,
        help="Path to a pickle file containing the GenBank records"
    )
    parser.add_argument(
        "--include", 
        default=None,
        help="Path to a file containing the accessions to include"
    )
    parser.add_argument(
        "--outgroup", 
        type=str,
        default=None,
        help="Accession of the outgroup sequence"
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
        type=float,
        default=None,
        help="Maximum fraction of ambiguous bases in the feature to be included in the output"
    )
    parser.add_argument(
        "--min_alignment_score",
        type=float,
        default=None,
        help="Minimum alignment score of the feature to be included in the output"
    )
    parser.add_argument(
        "--remove_duplicates",
        type=str,
        choices=["True", "False"],
        default="False",
        help="Remove duplicate sequences from the output"
    )
    parser.add_argument(
        "--remove_partial",
        type=str,
        choices=["True", "False"],
        default="True",
        help="Remove partial coding sequences from the output"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to filtered output accession csv file"
    )
    return parser.parse_args()


def filter_records(
        df,
        records, 
        include, 
        outgroup, 
        min_length, 
        max_length, 
        max_ambiguous, 
        min_score, 
        remove_duplicates, 
        remove_partial
    ):
    """
    Filter records based on the presence of certain features.

    Parameters
    ----------
    df: pd.DataFrame
        A dataframe containing the results of feature parsing.
    records: dict
        A dict of GenBank records.
    include: list
        A list of accessions to ensure are included in the output.
    outgroup: str
        Accession of the outgroup sequence.
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
    remove_partial: bool
        Whether to remove partial coding sequences from the output.

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

    print(f"""
            Filtering records with the following parameters:
            Minimum feature length: {min_length}
            Maximum feature length: {max_length}
            Maximum number of ambiguous positions: {max_ambiguous}
            Minimum alignment/identity score: {min_score}
            Remove duplicates: {remove_duplicates}
            Remove partial coding sequences: {remove_partial}
            Check for 'included' accessions: {"True" if include else "False"}
            Outgroup: {outgroup if outgroup else "None"}
            \n
          """
        )
    # Filter accessions based on quantitative parameters
    df = df.query("n_cds > 0")
    df = df.assign(
        feature_ambiguous=lambda df: df.apply(
                lambda x: sum(1 for base in x.sequence if base not in "ACGT") / len(x.sequence),
                axis=1
            )
    )
    filtered_df = (df
                   .query("feature_length >= @min_length and feature_length <= @max_length")
                   .query("feature_ambiguous <= @max_ambiguous")
                   .query("identity >= @min_score")
    )
    # Filter accessions based on boolean parameters
    if remove_duplicates:
        filtered_df = filtered_df.drop_duplicates(subset=['sequence'])
    if remove_partial:
        partial_accessions = remove_partial_cds(filtered_df, records)
        filtered_df = filtered_df[~filtered_df['accession'].isin(partial_accessions)]
    # Ensure that the 'include' accessions are present in the output
    if include:
        missing_accessions = set(include) - set(filtered_df['accession'].unique())
        if missing_accessions: 
            print(f"\nWarning: filtering parameters removed {len(missing_accessions)} 'include' sequences.")
            for accession in missing_accessions:
                print(f"\t> {accession}")
    if outgroup:
        if outgroup not in filtered_df['accession'].unique():
            filtered_df = pd.concat([filtered_df, df[df['accession'] == outgroup]])
            print("Warning: filtering parameters removed the outgroup sequence. Added it back.")

    print(f"\nFiltered out {len(df) - len(filtered_df)} records. {len(filtered_df)} records remain.")
    return filtered_df.reset_index(drop=True)


def remove_partial_cds(df, records):
    """
    Remove sequences with partial coding regions:
    - Partial coordinates in the GenBank record.
    - Coding regions are those that are missing a Start or Stop codons.
    - A sequence that isn't divisible by 3 with a remainder of 0.
    - A sequence whose translation doesn't match the annotated translation.

    Parameters
    ----------
    df: pd.DataFrame
        A dataframe containing the results of feature parsing.
    records: dict
        A dict of GenBank records.

    Returns
    -------
    list
        A list of accessions with partial coordinates in the GenBank record.
    """
    failed = 0
    failed_accessions = []
    for row in df.itertuples():
        accession = row.accession
        sequence = row.sequence
        for cds in records[accession].coding_regions:
            if cds.sequence == sequence:
                start_is_exact = isinstance(cds.feature.location.start, ExactPosition)
                end_is_exact = isinstance(cds.feature.location.end, ExactPosition)
                if not start_is_exact or not end_is_exact:
                    failed += 1
                    failed_accessions.append(accession)
                elif len(cds.sequence) %3 != 0:
                    failed += 1
                    failed_accessions.append(accession)
                elif cds.sequence[:3] not in ["ATG", "AUG"] or cds.sequence[-3:] not in ["TAA", "TAG", "TGA", "UAA", "UAG", "UGA"]:
                    failed += 1
                    failed_accessions.append(accession)
                elif not cds.check_translation():
                    print(f"Warning: sequence {accession} has a different translation than the annotated translation despite having exact positions and a start and stop codon.")
                    failed += 1
                    failed_accessions.append(accession)
    print(f"Found {failed} partial coding sequences.")
    return failed_accessions


def main():
    """Main entry point of the script."""
    args = parse_arguments()

    # Load the CSV file
    df = pd.read_csv(args.features)
    print(f"Loaded {len(df)} records from {args.features}")

    # Load in the included accessions
    if args.include:
        with open(args.include, 'r') as f:
            include = [line.strip() for line in f.readlines()]
        print(f"Loaded {len(include)} accessions from {args.include}")
    else:
        include = None
    # Load in the GenBank records
    with open(args.records, 'rb') as f:
        records = pickle.load(f)

    # Filter records based on the presence of certain features
    filtered_df = filter_records(df, 
                                 records,
                                 include,
                                 args.outgroup,
                                 args.min_feature_length, 
                                 args.max_feature_length, 
                                 args.max_ambiguous_positions, 
                                 args.min_alignment_score, 
                                 True if args.remove_duplicates == "True" else False,
                                 True if args.remove_partial == "True" else False
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