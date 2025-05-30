rule convert_sra_to_fastq_se:
    input: join(config['sra_dir'], "{accession}", "{accession}.sra")
    output:        
        temp(join(config['fastq_dir'], "{accession}", "{accession}.fastq.gz")),
    params:
        directory = join(config['fastq_dir'], "{accession}"),
    conda: 'intrahost'
    log: join(config['log_dir'], "convert_sra_to_fastq", "convert_sra_to_fastq_se_{accession}.log")
    shell:
        """
        fastq-dump \
            --gzip \
            --origfmt \
            --dumpbase \
            --skip-technical \
            --clip \
            --outdir {params.directory} \
            {input} \
            &>> {log}
        """


rule trim_all_reads_se:
    input:
        R1 = join(config['fastq_dir'], "{accession}", "{accession}.fastq.gz"),
    output:
        R1 = temp(join(config['trim_dir'], "{accession}", "{accession}.fastq.gz")),
        html = join(config['qc_dir'], "{accession}", "{accession}.fastp.html"),
        json = join(config['qc_dir'], "{accession}", "{accession}.fastp.json")
    conda: 'intrahost'
    threads: config['trim']['threads']
    log: join(config['log_dir'], "trim_all_reads_se", "trim_all_reads_se_{accession}.log")
    shell:
        """ 
        fastp \
            -w {threads} \
            -i {input.R1} \
            -o {output.R1} \
            --trim_poly_g \
            --trim_poly_x \
            --cut_front \
            --cut_tail \
            --cut_front_mean_quality 10 \
            --cut_tail_mean_quality 10 \
            --html {output.html} \
            --json {output.json} \
            &>> {log}
        """


rule filter_viral_reads_w_bbduk_se:
    input: 
        R1 = join(config['trim_dir'], "{accession}", "{accession}.fastq.gz"),
        reference = join(config['ref_dir'], "reference.fa"),
    output:
        R1 = temp(join(config['filter_dir'], "{accession}", "{accession}.fastq.gz")),
        stats = join(config['qc_dir'], "{accession}", "{accession}.filter.stats"),
    params:
        k = config['filter']['kmer_length'],
        hdist = config['filter']['hamming_distance'],
    conda: 'intrahost'
    threads: config['filter']['threads']
    log: join(config['log_dir'], "filter_viral_reads_w_bbduk_se", "filter_viral_reads_w_bbduk_se_{accession}.log")
    shell:
        """
        bbduk.sh \
            in={input.R1} \
            outm={output.R1} \
            ref={input.reference} \
            k={params.k} \
            hdist={params.hdist} \
            t={threads} \
            stats={output.stats} \
            &>> {log}
        """
