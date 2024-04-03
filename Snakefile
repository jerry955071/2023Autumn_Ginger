configfile: "configs/config.json"

# configure docker mounting options
# ==============================================================================
docker_mount_opt = ""
for volume in config["volumes"]:
    docker_mount_opt += "-v %s:%s:%s " % (
        volume["real"],
        volume["virtual"],
        volume["mode"]
    )
    if volume["is_workdir"]:
        docker_mount_opt += "-w %s " % volume["virtual"]

# define query function
# ==============================================================================
def query(**kwargs):
    for d in config["samples"]:
        for k, v in kwargs.items():
            if d[k] == v:
                return d
                
    raise ValueError("No such sample: %s" % kwargs)
    return None

# all
# ==============================================================================
rule all:
    input:
        # "outputs/rsem-dmat/all.genes.tpm.matrix",
        "outputs/rsem-dmat/all.rcounts.genes.matrix"

# GSEA
# ==============================================================================
# make blast database
rule makeblastdb:
    input:
        "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa"
    output:
        "references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa.psq"
    log:
        "logs/blastp/makeblastdb.log"
    shell:
        """
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name makeblastdb \
            biocontainers/blast:v2.2.31_cv2 \
                makeblastdb \
                    -in {input} \
                    -dbtype prot \
                    2> {log} \
                    1> {log}
        """


# blast ginger proteins against arabidopsis protein
rule blastp:
    threads: 24
    input:
        query="references/ginger_protein/ncbi_dataset/data/GCF_018446385.1/protein.faa",
        db="references/ath/Phytozome/PhytozomeV9/Athaliana/annotation/Athaliana_167_protein_primaryTranscriptOnly.fa"
    output:
        "outputs/blastp/ginger_vs_ath.tsv"
    log:
        "logs/blastp/ginger_vs_ath.log"
    shell:
        """
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name blastp \
            biocontainers/blast:v2.2.31_cv2 \
                blastp \
                    -query {input.query} \
                    -db {input.db} \
                    -out {output} \
                    -outfmt 6 \
                    -num_alignments 1 \
                    -num_threads {threads} \
                    2> {log} \
                    1> {log}
        """



# RSEM quantification
# ==============================================================================
# rule rsem_tpm_dmat:
#     conda:
#         "env/2023aut_ginger.yaml"
#     input:
#         expand(
#             "outputs/rsem-cal-expr/{id}/{id}.genes.results",
#             id=[d["id"] for d in config["samples"]]
#         )
#     output:
#         "outputs/rsem-dmat/all.genes.tpm.matrix"
#     log:
#         "logs/rsem/rsem_tpm_dmat.log"
#     shell:
#         """
#         rsem-generate-tpm-data-matrix {input} 1> {output} 2> {log}
#         """

rule rsem_dmat:
    conda:
        "env/2023aut_ginger.yaml"
    input:
        expand(
            "outputs/rsem-cal-expr/{id}/{id}.genes.results",
            id=[d["id"] for d in config["samples"]]
        )
    output:
        "outputs/rsem-dmat/all.rcounts.genes.matrix"
    log:
        "logs/rsem/rsem_dmat.log"
    shell:
        """
        rsem-generate-data-matrix {input} 1> {output} 2> {log}
        """

rule rsem_cal_exp:
    threads: 12
    conda:
        "env/2023aut_ginger.yaml"
    input:
        ref_dir="references/rsem/",
        fq1="outputs/afqc/{id}/{id}_1.fastq.gz",
        fq2="outputs/afqc/{id}/{id}_2.fastq.gz",
    output:
        dout=directory("outputs/rsem-cal-expr/{id}"),
        gene_res="outputs/rsem-cal-expr/{id}/{id}.genes.results"
    params:
        ref_prefix="references/rsem/rsem",
        output_prefix=lambda wildcards: f"outputs/rsem-cal-expr/{wildcards.id}/{wildcards.id}",
        unused_flags=" ".join(
            ["--strand-specific",]
        )
    log:
        "logs/rsem/rsem_cal_exp/{id}.log"
    shell:
        """
        # create output directory
        mkdir -p $(dirname {output})

        # calculate expression
        rsem-calculate-expression \
            --paired-end \
            --num-threads {threads} \
            --bowtie2 \
            --bowtie2-path $(dirname $(which bowtie2)) \
            --time \
            {input.fq1} \
            {input.fq2} \
            {params.ref_prefix} \
            {params.output_prefix} \
        2>&1 \
        > {log}
        """

rule rsem_prep_ref:
    threads: 4
    conda:
        "env/2023aut_ginger.yaml"
    input:
        gtf="references/ncbi_dataset/data/GCF_018446385.1/genomic_modified.gtf",
        fa="references/ncbi_dataset/data/GCF_018446385.1/GCF_018446385.1_Zo_v1.1_genomic.fna"
    output:
        dout=directory("references/rsem")
    params:
        ref_prefix="references/rsem/rsem"
    log:
        "logs/rsem/rsem_prep_ref.log"
    shell:
        """
        # create output directory
        mkdir -p {output.dout}

        # prepare reference
        rsem-prepare-reference \
            --num-threads {threads} \
            --bowtie2 \
            --bowtie2-path $(dirname $(which bowtie2)) \
            --gtf {input.gtf} \
            --polyA \
            {input.fa} \
            {params.ref_prefix} \
            2>&1 \
            > {log}
        """

# get ginger genome
# ==============================================================================
rule remove_ambigous:
    input:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic.gtf"
    output:
        "references/ncbi_dataset/data/GCF_018446385.1/genomic_modified.gtf"
    log:
        "logs/ncbi-dataset/remove_ambigous.log"
    shell:
        """
        # remove ambigous gene 'F6E76_pgp044'
        grep -v F6E76_pgp044 {input} > {output}

        # report ambigous gene
        echo "Ambigous gene:" > {log}
        grep F6E76_pgp044 {input} >> {log}
        """


rule ncbi_rehydrate:
    input:
        fin="references/ncbi_dataset.zip"
    output:
        dout=directory("references/ncbi_dataset/data")
    log:
        "logs/ncbi-dataset/ncbi_rehydrate.log"
    shell:
        """
        # unzip archive
        unzip -u -d $(dirname {input.fin}) {input.fin}

        # rehydrate
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_rehydrate \
            ccc/ncbi-datasets:20230926 \
                datasets rehydrate \
                    --directory $(dirname {input.fin}) \
                2> {log} \
                1> {log}
        """

rule ncbi_dload_dehydrated:
    output:
        fullname="references/ncbi_dataset.zip"
    log:
        "logs/ncbi-dataset/ncbi_dload_dehydrated.log"
    shell:
        """
        mkdir -p $(dirname {output.fullname})
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name ncbi_dload_dehydrated \
            ccc/ncbi-datasets:20230926 \
                datasets download genome accession GCF_018446385.1 \
                    --include gff3,gtf,genome,seq-report \
                    --filename {output.fullname} \
                    --dehydrated \
                2> {log} \
                1> {log}
        """
    


# QC
# ==============================================================================
rule fastp:
    threads: 4
    input:
        in1=lambda wildcards: query(id=wildcards.id)["fq1"],
        in2=lambda wildcards: query(id=wildcards.id)["fq2"],
    output:
        out1="outputs/afqc/{id}/{id}_1.fastq.gz",
        out2="outputs/afqc/{id}/{id}_2.fastq.gz",
        report_json="outputs/afqc/{id}/fastp.json",
        report_html="outputs/afqc/{id}/fastp.html"
    params:
        flags=" ".join([
            "--detect_adapter_for_pe",
            "--correction",
            "--cut_front",
            "--cut_tail",
            "--disable_trim_poly_g"
        ])
    log:
        "logs/fastp/{id}.log"
    shell:
        """
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            --name fastp_{wildcards.id} \
            biocontainers/fastp:v0.20.1_cv1 \
                fastp \
                    {params.flags} \
                    --in1 {input.in1} \
                    --in2 {input.in2} \
                    --out1 {output.out1} \
                    --out2 {output.out2} \
                    --json {output.report_json} \
                    --html {output.report_html} \
                1> {log} \
                2> {log} 
        """