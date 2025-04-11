#!/usr/bin/env python3
import argparse
import time
import pandas as pd
from Bio import Align
from Genbank import GenBankRecord


def parse_arguments():
    parser = argparse.ArgumentParser(description="Download and parse a list of GenBank sequences.")
    parser.add_argument(
        "--accessions", 
        required=True,
        help="Path to a file containing the accessions to download"
    )
    parser.add_argument(
        "--feature",
        required=True,
        help="The feature to extract from the GenBank records"
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="The reference accession used for feature alignment"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="The path to the output csv file"
    )
    parser.add_argument(
        "--include", 
        nargs="+",
        default=[],
        help="List of accessions to include (space-separated)"
    )
    parser.add_argument(
        "--exclude", 
        nargs="+",
        default=[],
        help="List of accessions to exclude (space-separated)"
    )
    return parser.parse_args()


def fetch_records(accessions):
    """
    Fetch GenBank records for a list of accessions.

    Parameters
    ----------
    accessions : list of str
        List of GenBank accessions to fetch.

    Returns
    -------
    dict
        A dictionary of GenBank records keyed by accession.
    """
    start_time = time.time()
    records = {}
    for accession in accessions:
        records[accession] = GenBankRecord(accession).fetch()
    end_time = time.time()
    print(f"Fetched {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
    return records


def algin_feature(feature, queries, type='nucleotide'):
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
    start_time = time.time()
    try:
        feature = [cds for cds in reference.coding_regions if cds.product == feature_name][0]
    except IndexError:
        raise ValueError(f"No feature found with name '{feature_name}'")
   
    matches = []
    for accession, record in records.items():
        queries = record.coding_regions
        match, score = algin_feature(feature, queries, type)
        if type == 'nucleotide':
            identity = ((score / 2) / len(feature.sequence)) * 100 if score != float('-inf') and score > 0 else 0
        else:
            _, reference_score = algin_feature(feature, feature, type)
            identity = (score / reference_score) * 100 if score != float('-inf') and score > 0 else 0
        name = match.product if match else "no coding sequence"
        matches.append((accession, identity, name, len(queries), match.sequence, match.translation))
    
    end_time = time.time()
    print(f"Parsed features from {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
    return pd.DataFrame(matches, columns=['accession', 'identity', 'feature', 'n_cds', 'sequence', 'translation'])


def extract_metadata(records):
    """
    Extract metadata from GenBank records.

    Parameters
    ----------
    records: dict of GenBankRecords
        A dictionary of GenBank records keyed by accession.

    Returns
    -------
    pd.DataFrame
        A dataframe of the metadata for each record
    """
    start_time = time.time()
    base_properties = [
        'accession',
        'url',
        'authors',
        'title',
        'journal',
        'paper_link',
        'submission',
        'strain',
        'organism',
        'host',
        'date',
        'location',
        'ambiguous',
        'length'
    ]
    geographic_properties = ['country', 'local', 'region', 'subregion']
    properties = base_properties + geographic_properties
    
    metadata = []
    for accession, record in records.items():
        # Get base properties
        record_metadata = [getattr(record, prop, None) for prop in base_properties]
        # Fetch the geographic properties from a REST API
        geographic_information = record.fetch_geographic_information()
        if geographic_information is None:
            print(f"Warning: Unable to fetch geographic information for {accession}")
            record_metadata.extend([
                getattr(record, 'country', "Unknown"),
                getattr(record, 'local', "Unknown"),
                "Unknown",  # region
                "Unknown"   # subregion
            ])
        else:
            record_metadata.extend([geographic_information.get(prop, "Unknown") for prop in geographic_properties])
        # Save the metadata for this record 
        metadata.append(record_metadata)

    end_time = time.time()
    print(f"Extracted metadata {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
    return pd.DataFrame(metadata, columns=properties)


def main():
    """Main entry point of the script."""
    args = parse_arguments()

    # Prepare the accessions to fetch
    accession_df = pd.read_csv(args.accessions, header=None, names=['accession'])
    accessions = list(accession_df.accession.unique())
    # Add the reference accession if it's not there
    if args.reference not in accessions:
        accessions.append(args.reference)
    # Add the accessions to include
    if args.include:
        for acc in args.include:
            if acc not in accessions:
                accessions.append(acc)
    # Remove the accessions the exclude
    if args.exclude:
        for acc in args.exclude:
            if acc not in accessions:
                accessions.remove(acc)

    # Fetch the records for each accession
    print(f"Fetching {len(accessions)} records from GenBank...")
    records = fetch_records(accessions)

    # Parse the features from the records
    print(f"Parsing features from {len(accessions)} records...")
    feature_df = parse_features(records, records[args.reference], args.feature, type='nucleotide')
    
    # Parse the metadata from the records
    print(f"Extracting metadata for {len(accessions)} records...")
    metadata_df = extract_metadata(records)

    # Export to CSV
    joined_df = pd.merge(metadata_df, feature_df, on='accession', how='inner')
    joined_df.to_csv(args.output, index=False)
    print(f"Finished. Results exported to {args.output}")


if __name__ == "__main__":
    main()