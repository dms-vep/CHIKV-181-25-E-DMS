def single_ended(accession):
    """
    Returns `True` if the Library Layout is single-end reads.
    """
    # Get the Library Layout
    layout = pd.read_csv(config["fetch"]["sra_table"]).set_index('acc').at[accession, 'LibraryLayout']
    if layout == "SINGLE":
        return True
    else:
        return False


def get_processed_fastqs(wildcards):
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config['filter_dir'], "{accession}", "{accession}.fastq.gz"), accession=wildcards.accession)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config['filter_dir'], "{accession}", "{accession}_1.fastq.gz"),
                   join(config['filter_dir'], "{accession}", "{accession}_2.fastq.gz")], accession=wildcards.accession)
     

def get_minimap_preset(accession):
    platform = pd.read_csv(config["fetch"]["sra_table"]).set_index('acc').at[accession, 'Platform']
    preset_map = {
        'OXFORD_NANOPORE': 'map-ont',
        'PACBIO_SMRT': 'map-pb',
        'ILLUMINA': 'sr'
    }
    return preset_map.get(platform, 'sr') 