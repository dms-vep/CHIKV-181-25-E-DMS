#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Create codon alignment based on protein alignment.")
    parser.add_argument(
        "--orfs", 
        required=True,
        help="Path to the FASTA output from EMBOSS getorf"
    )
    parser.add_argument(
        "--sequences", 
        required=True,
        help="Path to the FASTA file of the untranslated sequences"
    )
    parser.add_argument(
        "--protein", 
        required=True,
        help="Path to the output file of the protein sequences with strain headers"
    )
    parser.add_argument(
        "--codon", 
        required=True,
        help="Path to the output file of the codon sequences with strain headers"
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    print(f"Extracting protein and codon sequences from {args.orfs} and untranslated sequences {args.sequences}\n")
    # Read in the untranslated sequences into a dictionary by strain name
    sequences = {seq.name: seq for seq in SeqIO.parse(args.sequences, "fasta")}
    # Got through all of the ORFs and see if more than one is identified for a given accession.
    identified_orfs = defaultdict(list)
    for protein in SeqIO.parse(args.orfs, "fasta"):   
        description = protein.description
        # Extract coordinates of the ORF in the original sequence
        coordinates = re.search(r'\[(\d+)\s*-\s*(\d+)\]', description.strip())
        if coordinates:
            start = int(coordinates.group(1))
            end = int(coordinates.group(2))
        else:
            print(f"Warning: couldn't find coordinates in {protein.id}")
            continue
        # Extract the accession number and the ORF number 
        id = re.search(r'^(.+)_(\d+)$', description.split(" ")[0])
        if id:
            accession = id.group(1)
            number = id.group(2)
        else:
            print(f"Warning: couldn't find accession in {protein.id}")
            continue
        # Store the ORF in a dictionary
        identified_orfs[accession].append((number, start, end, protein.seq))

    # Extract and check the protein sequences and codon translations
    protein_records = []
    codon_records = []
    for accession, orfs in identified_orfs.items():
        # Check if there is exactly one ORF for the accession
        if len(orfs) > 1:
            print(f"Warning: {accession} has multiple ORFs: {orfs}")
            continue   
        # Check that the parsed accession matches the nucleotide sequence
        try:
            nucleotide = sequences[accession]
        except KeyError:
            print(f"Warning: {accession} not found in nucleotide sequences")
            continue
        # Check that coordinates in the nucleotide sequence match the protein sequence
        _, start, end, protein = orfs[0]
        if not str(nucleotide[start-1:end].translate().seq) == str(protein):
            print(f"Warning: {accession} the translation of the sequence of {start}-{end} does not match extracted protein sequence")
            for i, (a, b) in enumerate(zip(str(nucleotide[start-1:end].translate().seq), str(protein))):
                if a != b:
                    print(f"\tMismatch at position {i+1}: {a} != {b}")
            continue

        # Store the protein record
        protein_record = SeqIO.SeqRecord(
            protein,
            id=accession,
            name=accession,
            description=""
        )
        protein_records.append(protein_record)
        
        # Store the codon record
        codon_record = SeqIO.SeqRecord(
            nucleotide[start-1:end].seq,
            id=accession,
            name=accession,
            description=""
        )
        codon_records.append(codon_record)

    # Write the protien records to a FASTA file
    SeqIO.write(protein_records, args.protein, "fasta")
    # Write the codon records to a FASTA file
    SeqIO.write(codon_records, args.codon, "fasta")


if __name__ == "__main__":
    main()