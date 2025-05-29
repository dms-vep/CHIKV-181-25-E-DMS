rule convert_sra_to_fastq_pe:
    input: join(config['sra_dir'], "{accession}", "{accession}.sra")
    output:        
        R1 = join(config['fastq_dir'], "{accession}", "{accession}_1.fastq.gz"),
        R2 = join(config['fastq_dir'], "{accession}", "{accession}_2.fastq.gz"),
    params:
        directory = join(config['fastq_dir'], "{accession}"),
    conda: 'intrahost'
    log: join(config['log_dir'], "convert_sra_to_fastq", "convert_sra_to_fastq_pe_{accession}.log")
    shell:
        """
        fastq-dump \
            --gzip \
            --dumpbase \
            --clip \
            --skip-technical \
            --readids \
            --origfmt \
            --split-spot \
            --split-3 \
            --outdir {params.directory} \
            {input} \
            &>> {log}
        """


rule trim_all_reads_pe:
    input:
        R1 = join(config['fastq_dir'], "{accession}", "{accession}_1.fastq.gz"),
        R2 = join(config['fastq_dir'], "{accession}", "{accession}_2.fastq.gz"),
    output:
        R1 = join(config['trim_dir'], "{accession}", "{accession}_1.fastq.gz"),
        R2 = join(config['trim_dir'], "{accession}", "{accession}_2.fastq.gz"),
        html = join(config['qc_dir'], "{accession}", "{accession}.fastp.html"),
        json = join(config['qc_dir'], "{accession}", "{accession}.fastp.json")
    conda: 'intrahost'
    threads: config['trim']['threads']
    log: join(config['log_dir'], "trim_all_reads_pe", "trim_all_reads_pe_{accession}.log")
    shell:
        """ 
        fastp \
            -w {threads} \
            -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
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


rule filter_viral_reads_w_bbduk_pe:
    input: 
        R1 = join(config['trim_dir'], "{accession}", "{accession}_1.fastq.gz"),
        R2 = join(config['trim_dir'], "{accession}", "{accession}_2.fastq.gz"),
        reference = join(config['ref_dir'], "reference.fa"),
    output:
        R1 = join(config['filter_dir'], "{accession}", "{accession}_1.fastq.gz"),
        R2 = join(config['filter_dir'], "{accession}", "{accession}_2.fastq.gz"),
        stats = join(config['qc_dir'], "{accession}", "{accession}.filter.stats"),
    params:
        k = config['filter']['kmer_length'],
        hdist = config['filter']['hamming_distance'],
    conda: 'intrahost'
    threads: config['filter']['threads']
    log: join(config['log_dir'], "filter_viral_reads_w_bbduk_pe", "filter_viral_reads_w_bbduk_pe_{accession}.log")
    shell:
        """
        bbduk.sh \
            in1={input.R1} \
            in2={input.R2} \
            outm1={output.R1} \
            outm2={output.R2} \
            ref={input.reference} \
            k={params.k} \
            hdist={params.hdist} \
            t={threads} \
            stats={output.stats} \
            &>> {log}
        """
