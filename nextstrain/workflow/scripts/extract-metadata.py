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
        "--virus", 
        type=str,
        default=None,
        help="Abbreviation of the virus name for strain name (e.g. 'CHIKV' for Chikungunya virus)"
    )
    parser.add_argument(
        "--country_mapping",
        type=json.loads,
        default=None,
        help="Custom country name mappings"
    )
    parser.add_argument(
        "--host_mapping",
        type=json.loads,
        default=None,
        help="Custom host name mappings"
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
            if record.country == "?":
                print(f"Warning: Country is Unknown for {accession}")
            else:
                print(f"Warning: Unable to fetch geographic information for {record.country}: {accession}")
            record_metadata.extend([
                getattr(record, 'country', "?"),
                getattr(record, 'local', "?"),
                "?",  # region
                "?"   # subregion
            ])
        else:
            record_metadata.extend([geographic_information.get(prop, "?") for prop in geographic_properties])
        # Save the metadata for this record 
        metadata.append(record_metadata)
    return pd.DataFrame(metadata, columns=properties)


def make_unique_id(df, virus):
    """
    Make a unique 'strain' identifier for each record in the dataframe.
    Check if the strains are unique in the dataframe.

    If they aren't unique, print a warning.
    Make a new columns called 'name'

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe containing the metadata.

    Returns
    -------
    pd.DataFrame
        The dataframe with the new 'name' column.
    """
    # Check if the strains are unique
    if df.strain.duplicated().any():
        print("\nWarning: Strains are not unique in the dataframe:")
        print("A number will be appended to each strain duplicate to ensure uniqueness.")
        duplicated_strains = df[df.strain.duplicated()].strain.unique()
        # Loop through the duplicated strains and print the accessions for each
        for strain in duplicated_strains:
            print(f"\nStrain: {strain}\n----------")
            duplicate_strain_accessions = df[df.strain == strain].accession.unique()
            for i, acc in enumerate(duplicate_strain_accessions):
                print(f"\tAccession: {acc}")
                df.loc[df.accession == acc, 'strain'] = f"{strain}_{i+1}"
    # Make a new column called 'name'
    if virus:
        df['name'] = df.apply(lambda row: f"{virus}/{row['strain']}/{row['host']}/{row['date'][:4]}", axis=1)
    else:
        df['name'] = df.apply(lambda row: f"{row['strain']}/{row['host']}/{row['date'][:4]}", axis=1)
    return df


def check_accessions(df):
    """
    Ensure that all accessions are unique in the dataframe.

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe containing the metadata.
    """
    if df.accession.duplicated().any():
        print("Warning: Accessions are not unique in the dataframe.")
        print(df[df.accession.duplicated()].accession.unique())


def check_hosts(df, host_mapping=None):
    """
    Print the unique host species in the dataframe.

    Parameters
    ----------
    df: pd.DataFrame
        The dataframe containing the metadata.
    host_mapping: dict, optional
        A mapping of host species to their correct names.

    Returns
    -------
    pd.DataFrame
        The dataframe with the host species fixed.
    """
    # Print the unique species in the dataframe
    print("Unique species in the dataframe before formatting:")
    for host in df.host.unique():
        print(host)

    # If a host mapping is provided, use it to fix the host species
    if host_mapping:
        df['host'] = df['host'].map(lambda x: host_mapping.get(x, x))

    # Print the unique species in the dataframe after formatting
    print("\nUnique species in the dataframe after formatting:")
    for host in df.host.unique():
        print(host)
    
    return df


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

    # Check the strains, accessions, and species
    metadata_df = make_unique_id(metadata_df, args.virus)
    metadata_df = check_hosts(metadata_df, host_mapping=args.host_mapping)
    check_accessions(metadata_df)

    # Export the metadata to a CSV file
    metadata_df.to_csv(args.output, index=False)
    end_time = time.time()
    print(f"Extracted metadata from {len(records)} records from GenBank in {end_time - start_time:.4f} seconds.\n")
    print(f"Finished. Results exported to {args.output}")


if __name__ == "__main__":
    main()