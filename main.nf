nextflow.enable.dsl=2

workflow {
    // 1. preprocess step channels (run in parallel)
    // samples_ch = Channel.of(params.normal, params.tumor)

    // 2. preprocess step
    // PREPROCESS(samples_ch)

    // 3. merge step
    // only for normal
    // GENERATE_BAMS_LIST(PREPROCESS.out.collect())
    // MERGE_NORMAL(GENERATE_BAMS_LIST.out)

    normal_bam_ch  = Channel.fromPath("${params.results}/${params.normal}/${params.normal}.bam").collect()
    normal_bai_ch  = Channel.fromPath("${params.results}/${params.normal}/${params.normal}.bam.bai").collect()
    tumor_bam_ch   = Channel.fromPath("${params.results}/${params.tumor}/${params.tumor}.bam").collect()
    tumor_bai_ch   = Channel.fromPath("${params.results}/${params.tumor}/${params.tumor}.bam.bai").collect()
    genome_fai_ch  = Channel.fromPath("${params.genome}.fai").collect()

    // 4. Exclude / Include step
    // depends on step 3
    PREPARE_BEDS(normal_bam_ch, normal_bai_ch, genome_fai_ch)

    // 5. binaries (Lumpy, GATK, CaVEMan, ASCAT)
    // parallel run
    RUN_LUMPY_NORMAL(normal_bam_ch, normal_bai_ch, PREPARE_BEDS.out.exclude_bed)
    RUN_LUMPY_TUMOR(tumor_bam_ch, tumor_bai_ch, PREPARE_BEDS.out.exclude_bed)

    RUN_STATS_NORMAL(normal_bam_ch, normal_bai_ch)
    RUN_STATS_TUMOR(tumor_bam_ch, tumor_bai_ch)

    
    RUN_GATK(tumor_bam_ch, tumor_bai_ch, normal_bam_ch, normal_bai_ch, PREPARE_BEDS.out.include_bed)
    RUN_CAVEMAN(tumor_bam_ch, tumor_bai_ch, normal_bam_ch, normal_bai_ch, PREPARE_BEDS.out.include_bed, PREPARE_BEDS.out.filtered_fai)
    RUN_ASCAT(tumor_bam_ch, tumor_bai_ch, normal_bam_ch, normal_bai_ch)
}

/* ==========================================
   Preprocess
   ========================================== */

process PREPROCESS {
    tag "${sample_id}"

    input:
    val sample_id

    script:
    """
    bash ${params.scripts}/preprocess.sh ${sample_id} noERX noFastq1 noFastq2 ${params.genome}
    """
}

process GENERATE_BAMS_LIST {
    input:
    val ready_signals // ensures the preprocess finished

    output:
    path "bams_${params.normal}.txt"

    script:
    """
    ls ${params.results}/${params.normal}/*/*.dedup.bam > bams_${params.normal}.txt
    """
}

process MERGE_NORMAL {
    input:
    path bams_list

    output:
    val true // merge finished signal

    script:
    """
    bash ${params.scripts}/merge.sh ${params.normal} ${bams_list}
    """
}

process PREPARE_BEDS {
    input:
    path normal_bam
    path normal_bai
    path genome_fai

    output:
    path "exclude_${params.normal}.bed", emit: exclude_bed
    path "include_${params.normal}.bed", emit: include_bed
    path "filtered_genome.fai", emit: filtered_fai

    script:
    """
    module load SAMtools
    samtools view -H ${normal_bam} | \\
    awk -F'[\\t:]' '\$1=="@SQ" && \$3 !~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ {print \$3"\\t1\\t"\$5}' > exclude_${params.normal}.bed

    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ {print \$1"\\t1\\t"\$2}' ${params.genome}.fai > include_${params.normal}.bed

    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/' ${genome_fai} > filtered_genome.fai
    """
}

process RUN_LUMPY_NORMAL {
    conda "${params.lumpy_env}"

    input:
    path bam
    path bai
    path exclude_bed

    script:
    """
    module load miniconda
    bash ${params.scripts}/sv.sh ${params.normal} ${bam} ${params.genome} ${exclude_bed}
    """
}

process RUN_LUMPY_TUMOR {
    conda "${params.lumpy_env}"

    input:
    path bam
    path bai
    path exclude_bed

    script:
    """
    module load miniconda
    bash ${params.scripts}/sv.sh ${params.tumor} ${bam} ${params.genome} ${exclude_bed}
    """
}

process RUN_STATS_NORMAL {
    conda "${params.lumpy_env}"

    input:
    path bam
    path bai

    script:
    """
    module load miniconda
    export REF_PATH=${params.genome}
    bash ${params.scripts}/stats.sh ${params.normal} ${bam}
    """
}

process RUN_STATS_TUMOR {
    conda "${params.lumpy_env}"

    input:
    path bam
    path bai

    script:
    """
    module load miniconda
    export REF_PATH=${params.genome}
    bash ${params.scripts}/stats.sh ${params.tumor} ${bam}
    """
}

process RUN_GATK {
    publishDir "${params.results}/${params.tumor}", mode: 'copy'

    input:
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path include_bed

    output:
    path "gatk-${params.tumor}.vcf"

    script:
    def tumor_sm  = params.tumor_sample_id ?: params.tumor
    def normal_sm = params.normal_sample_id ?: params.normal
    """
    module load GATK
    gatk Mutect2 -R ${params.genome} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -L ${include_bed} \\
        -tumor ${tumor_sm} \\
        -normal ${normal_sm} \\
        --pair-hmm-implementation LOGLESS_CACHING \\
        -O "gatk-${params.tumor}.vcf"
    """
}

process RUN_CAVEMAN {
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${params.tumor}", mode: 'copy'

    output:
    path "caveman"

    input:
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path include_bed
    path filtered_fai

    script:
    """
    caveman.pl -o caveman \\
        -r ${filtered_fai} \\
        -tb ${tumor_bam} \\
        -nb ${normal_bam} \\
        -ig ~/hg38-blacklist.v2.bed -tc ~/empty.txt -td 2 -nc ~/empty.txt -nd 2 \\
        -s Human -sa GRCh38 -b ~/empty.txt \\
        -in ~/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \\
        -st genome -u ~/empty_dir -t 96 -noflag
    """
}

process RUN_ASCAT {
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${params.tumor}", mode: 'copy'

    input:
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai

    output:
    path "ascat"

    script:
    """
    ascat.pl -o ascat \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -r ${params.genome} -pr WGS -g XY -gc chrY \\
        -sg ~/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv -c 8
    """
}