#!/usr/bin/env python3
import argparse
import time
import pickle
import pandas as pd
from Bio import Align
from Genbank import GenBankRecord


def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse a sequence feature from a dictionary GenBank records.")
    parser.add_argument(
        "--input", 
        required=True,
        help="Path to a pickle object file containing the GenBank records"
    )
    parser.add_argument(
        "--feature",
        required=True,
        help="The feature to extract from the GenBank records"
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="The reference accession used to identify the feature"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="The path to the output csv file"
    )
    return parser.parse_args()


def align_feature(feature, queries, type='nucleotide'):
    """
    Align a set of query coding sequences to a reference coding sequence.
    Either the nucleotide or protein sequences based on specified type.

    Parameters
    ----------
    feature : Genbank.CDS
        The reference feature to align against.
    queries : list of Genbank.CDS
        A list of query sequences to align.
    type : str
        The type of alignment to perform ('nucleotide' or 'protein').

    Returns
    -------
    tuple
        The best matching query sequence and its score.
    """
    aligner = Align.PairwiseAligner(scoring="blastn" if type == 'nucleotide' else "blastp")
    best_score = float('-inf')
    best_query = None
    for query in queries:
        if type == 'nucleotide':
            alignments = aligner.align(feature.sequence, query.sequence)
        else:
            alignments = aligner.align(feature.translation, query.translation)
        # Select the best alignment
        alignment = alignments[0]
        if alignment is None:
            continue
        if alignment.score > best_score:
            best_score = alignment.score
            best_query = query
    return best_query, best_score


def parse_features(records, reference, feature_name, type='nucleotide'):
    """
    Parse features from GenBank records and align them to a reference feature.

    This will only work for cases where the sequence you're parsing is a CDS.

    An 'identity' score is calculated based on the alignment. This is relative
    to the alignment of the reference to itself. It's close to percent identity
    for nucleotide sequences and a relative score for protein sequences.

    Parameters
    ----------
    records : dict
        A dictionary of GenBank records keyed by accession.
    reference : Genbank.CDS
        The reference feature to align against.
    feature_name : str
        The feature to extract from the GenBank records.
    type : str
        The type of alignment to perform ('nucleotide' or 'protein').

    Returns
    -------
    pd.DataFrame
        A table containing the accession, identity, feature name, sequence, and translation.
    """
    try:
        feature = [cds for cds in reference.coding_regions if cds.product == feature_name][0]
    except IndexError:
        raise ValueError(f"No feature in the reference found with name '{feature_name}'")
   
    matches = []
    for accession, record in records.items():
        queries = record.coding_regions
        match, score = align_feature(feature, queries, type)
        if type == 'nucleotide':
            identity = ((score / 2) / len(feature.sequence)) * 100 if score != float('-inf') and score > 0 else 0
        else:
            _, reference_score = align_feature(feature, feature, type)
            identity = (score / reference_score) * 100 if score != float('-inf') and score > 0 else 0
        name = match.product if match else "no coding sequence"
        sequence = match.sequence if match else ""
        translation = match.translation if match else ""
        matches.append((accession, identity, name, len(queries), len(sequence), len(translation), sequence, translation))
    
    return pd.DataFrame(matches, columns=['accession', 'identity', 'feature', 'n_cds',  'feature_length', 'feature_translation_length', 'sequence', 'translation'])


def main():
    """Main entry point of the script."""
    args = parse_arguments()
    start_time = time.time()

    # Read in the pickle file containing the GenBank records
    with open(args.input, 'rb') as f:
        records = pickle.load(f)
    
    # Check that the reference is present in the records
    if args.reference not in records:
        raise ValueError(f"Reference accession '{args.reference}' not found in the records.")

    print(f"Parsing features from {len(records)} records...")
    feature_df = parse_features(records, records[args.reference], args.feature, type='nucleotide')

    # Export the results to a CSV file
    feature_df.to_csv(args.output, index=False)
    end_time = time.time()
    print(f"Parsed features from {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
    print(f"Finished. Results exported to {args.output}")


if __name__ == "__main__":
    main()