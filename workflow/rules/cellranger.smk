rule get_sra_10x:
    params:
        url = lambda wc: pep.get_sample(wc.sample).sra_accession,
        dir = "results/get-sra-10x/"
    output:
        temp(directory("results/get-sra-10x/{sample}"))
    shell:
        """
        wget -P {params.dir}/{wildcards.sample} {params.url}
        """

checkpoint fastq_dump_10x:
    input:
        rules.get_sra_10x.output
    output:
        temp(directory("results/fastq-dump-10x/{sample}/"))
    params:
        run_id = lambda wc: basename(pep.get_sample(wc.sample).sra_accession[0])
    threads:
        16
    resources:
        time=60,
        mem=config.get("FASTERQDUMP_MEM", 8000),
        cpus=16
    log:
        "results/logs/fastq-dump-10x/{sample}.log"
    conda:
        "../envs/sratools.yaml"
    shell:
       """
       fasterq-dump --mem {resources.mem}MB -s -S --include-technical -e {threads} -O {output} {input}/{params.run_id} 2> {log}
       """

def determine_layout(wildcards, dir = "results/fastq-dump-10x"):
    checkpoint_output = checkpoints.fastq_dump_10x.get(**wildcards).output[0]
    run_id = basename(pep.get_sample(wildcards.sample).sra_accession[0])
    base = "{d}/{s}/{r}".format(d=dir,s=wildcards.sample, r=run_id)
    glc = glob_wildcards(base + "_{i}.fastq")
    op = expand("{b}_{i}.fastq", b=base,i=glc.i)
    return(op)

rule rename_fqs_10x:
    input:
        determine_layout
    output:
        r1=temp("results/fastq-rename-10x/{sample}/{sample}_S1_L001_R1_001.fastq"),
        r2=temp("results/fastq-rename-10x/{sample}/{sample}_S1_L001_R2_001.fastq"),
    log:
        "results/logs/fastq-rename-10x/{sample}.log"
    run:
        s = wildcards.sample
        run_id = basename(pep.get_sample(s).sra_accession[0])
        if (len(input) == 2):
            shell("cp results/fastq-dump-10x/{s}/{r}_1.fastq {o}".format(s=s, r=run_id, o=output.r1)),
            shell("cp results/fastq-dump-10x/{s}/{r}_2.fastq {o}".format(s=s, r=run_id, o=output.r2))
        elif (len(input) == 3):
            shell("cp results/fastq-dump-10x/{s}/{r}_2.fastq {o}".format(s=s, r=run_id, o=output.r1)),
            shell("cp results/fastq-dump-10x/{s}/{r}_3.fastq {o}".format(s=s, r=run_id, o=output.r2))

rule cellranger_mkref:
    input:
        fa = custom_genome('results/custom-genome/combined.fasta'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        directory("results/cellranger/idx")
    params:
        mem=int(config.get("CELLRANGER_MKREF_MEM",16000)/1000.0)
    resources:
        time="18:00:00",
        mem=config.get("CELLRANGER_MKREF_MEM",16000),
        cpus=1
    shell:
        """
        cd results/cellranger
        ~/cellranger-3.1.0/cellranger mkref --genome=idx \
            --fasta={input.fa} \
	    --genes={input.gtf} \
	    --memgb={params.mem}
        """

rule cellranger_count:
    input:
        idx = rules.cellranger_mkref.output,
        fqs = rules.rename_fqs_10x.output,
    output:
        "results/cellranger/counts-{sample}/outs/filtered_feature_bc_matrix.h5"
    threads:
        config.get("CELLRANGER_COUNT_CPUS",32)
    params:
        ver = config.get("CELLRANGER_VERSION",'3.1.0'),
        mem = int((config.get("CELLRANGER_COUNT_MEM",128000)-5000)/1000.0)
    resources:
        time="18:00:00",
        mem=config.get("CELLRANGER_COUNT_MEM",128000),
        cpus=config.get("CELLRANGER_COUNT_CPUS",32),
    shell:
        """
        cd results/cellranger \
        rm -rf counts-{wildcards.sample} \
        ~/cellranger-{params.ver}/cellranger count --id=counts-{wildcards.sample} \
        --sample={wildcards.sample} \
        --fastqs=../fastq-rename-10x/{wildcards.sample} \
	    --transcriptome=idx \
	    --localcores={threads} \
	    --localmem={params.mem}
        """
