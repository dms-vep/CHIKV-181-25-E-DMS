#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Genbank import GenBankRecord


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract the relevant metadata from GenBank records.")
    parser.add_argument(
        "--reference",
        type=str, 
        required=True,
        help="Accession of the reference GenBank record"
    )
    parser.add_argument(
        "--product", 
        type=str, 
        required=True,
        help="Name of the CDS product to extract for the reference"
    )
    parser.add_argument(
        "--name",
        type=str,
        required=True,
        help="New name for the coding sequence"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to the output GenBank file"
    )
    return parser.parse_args()


def genbank_from_cds(accession, product, name, include_mat_peptides=True):
    """
    Make a new GenBank reference for a subset of the original GenBank record.

    Parameters
    ----------
    accession : str
        The accession of the reference GenBank record.
    product : str
        The name of the CDS product to extract.
    name : str
        The name of the new GenBank record.
    include_mat_peptides : bool, optional
        Whether to include the mat_peptide features in the new record. Default is True.
    """
    # Format the name
    name = name.replace(" ", "_")
    # Fetch the GenBank record
    reference = GenBankRecord(accession).fetch()
    # Initialize a list of features for the new record
    features = []
    # Find the main CDS with the specified product
    for cds in reference.coding_regions:
        if cds.product == product:
            features.append(cds.feature)
    if len(features) == 0:
            raise ValueError(f"CDS with product {product} not found in record {accession}")
    # Get position information for the coding sequence
    start = int(features[0].location.start)
    end = int(features[0].location.end)
    strand = features[0].location.strand
    length = end - start
    # Collect a list of mature peptide features within the cds region
    if include_mat_peptides:
        for feature in reference.record.features:
            if feature.type == "mat_peptide":
                if feature.location.start >= start and feature.location.end <= end:
                    features.append(feature)
    # Create a new GenBank record
    record = SeqRecord(
        features[0].location.extract(reference.record).seq, 
        id=name,
        name=name,
        description="",
        annotations={
            "molecule_type": reference.record.annotations["molecule_type"]
        }
    )
    # Add the source feature to the new record
    source = (
        SeqFeature(
            FeatureLocation(start=0, end=length, strand=strand),
            type="source",
            qualifiers=features[0].qualifiers
        )
    )
    source.qualifiers["product"] = name
    record.features.append(source)
    # Add the other features to the new record
    for i, feature in enumerate(features):
        relative_start = feature.location.start - start
        relative_end = feature.location.end - start
        if i == 0:
            gene_name = name
            translation = feature.extract(reference.record.seq).translate(to_stop=True)
            feature.qualifiers["product"] = gene_name
            feature.qualifiers["gene"] = gene_name
            feature.qualifiers["locus_tag"] = gene_name
        else:
            gene_name = feature.qualifiers.get("product", [""])[0].replace(" ", "_")
            translation = feature.extract(reference.record.seq).translate(to_stop=True)
            feature.qualifiers["gene"] = gene_name
            feature.qualifiers["locus_tag"] = gene_name
        if "translation" not in feature.qualifiers:
            feature.qualifiers["translation"] = [str(translation)]
        cds = SeqFeature(FeatureLocation(start=relative_start, end=relative_end), type="CDS", qualifiers=feature.qualifiers)
        record.features.append(cds)

    return record


def main():
    """Main entry point of the script."""
    args = parse_arguments()
    record = genbank_from_cds(args.reference, args.product, args.name)
    print(f"Created new GenBank record with {len(record.features)} features at {args.output}:")
    for feature in record.features:
        print(f"{feature.type}: {feature.qualifiers.get('gene')}")
    with open(args.output, "w") as handle:
        SeqIO.write(record, handle, "genbank")


if __name__ == "__main__":
    main()