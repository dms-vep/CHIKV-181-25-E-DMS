#!/usr/bin/env python3
import argparse
import json
import time
import pickle
import pandas as pd
from Genbank import GenBankRecord

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract the relevant metadata from GenBank records.")
    parser.add_argument(
        "--records", 
        required=True,
        help="Path to a pickle object file containing the GenBank records"
    )
    parser.add_argument(
        "--filtered", 
        required=True,
        help="Path to a CSV containing the filtered GenBank records"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="The path to the output metadata csv file"
    )
    parser.add_argument(
        "--country_mapping",
        type=json.loads,
        default=None,
        help="Custom country name mappings"
    )
    return parser.parse_args()


def extract_metadata(records, country_mapping=None):
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
        geographic_information = record.fetch_geographic_information(
            country_mapping=country_mapping
        )
        if geographic_information is None:
            if record.country == "Unknown":
                print(f"Warning: Country is Unknown for {accession}")
            else:
                print(f"Warning: Unable to fetch geographic information for {record.country}: {accession}")
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
    return pd.DataFrame(metadata, columns=properties)


def main():
    """Main entry point of the script."""
    args = parse_arguments()
    start_time = time.time()

    # Read in the pickle file containing the GenBank records
    with open(args.records, 'rb') as f:
        records = pickle.load(f)

    # Read in the filtered records
    filtered_df = pd.read_csv(args.filtered)
    filtered_accessions = list(filtered_df.accession.unique())

    # Filter the records to only include the filtered accessions
    records = {acc: record for acc, record in records.items() if acc in filtered_accessions}

    # Parse the metadata from the records
    print(f"Extracting metadata for {len(records)} records...")
    metadata_df = extract_metadata(records, country_mapping=args.country_mapping)

    # Export the metadata to a CSV file
    metadata_df.to_csv(args.output, index=False)
    end_time = time.time()
    print(f"Extracted metadata from {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
    print(f"Finished. Results exported to {args.output}")


if __name__ == "__main__":
    main()