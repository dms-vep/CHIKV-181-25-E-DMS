#!/usr/bin/env python3
import argparse
import time
import pickle
import pandas as pd
from Genbank import GenBankRecord


def parse_arguments():
    parser = argparse.ArgumentParser(description="Download and a list of GenBank sequences.")
    parser.add_argument(
        "--accessions", 
        required=True,
        help="Path to a file containing the accessions to download"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="The path to the output pickle file"
    )
    parser.add_argument(
        "--include", 
        help="Path to a file containing the accessions to ensure are included"
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
    records = {}
    for accession in accessions:
        records[accession] = GenBankRecord(accession).fetch()
    return records


def main():
    """Main entry point of the script."""
    args = parse_arguments()

    # Prepare the accessions to fetch
    accession_df = pd.read_csv(args.accessions, header=None, names=['accession'])
    accessions = list(accession_df.accession.unique())
   
    # Add the accessions to include
    if args.include:
        include_df = pd.read_csv(args.include, header=None, names=['accession'])
        include_accessions = list(include_df.accession.unique())
        for acc in include_accessions:
            if acc not in accessions:
                accessions.append(acc)
    
    # Remove the accessions the exclude
    if args.exclude:
        for acc in args.exclude:
            if acc in accessions:
                accessions.remove(acc)
    
    # Fetch the records for each accession
    start_time = time.time()
    print(f"Fetching {len(accessions)} records from GenBank...")
    records = fetch_records(accessions)
    end_time = time.time()
    print(f"Fetched {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
        
    # Save the records dictionary as a pickle file
    with open(args.output, 'wb') as f:
        pickle.dump(records, f)

    print(f"Finished. Results exported to {args.output}")


if __name__ == "__main__":
    main()