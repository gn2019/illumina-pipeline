nextflow.enable.dsl=2

workflow {
    // 1. preprocess step channels (run in parallel)
    samples_ch = Channel.of(params.normal, params.tumor)

    // 2. preprocess step
    PREPROCESS(samples_ch)

    // 3. merge step
    // only for normal
    GENERATE_BAMS_LIST(PREPROCESS.out.collect())
    MERGE_NORMAL(GENERATE_BAMS_LIST.out)

    // 4. Exclude / Include step
    // depends on step 3
    PREPARE_BEDS(MERGE_NORMAL.out)

    // 5. binaries (Lumpy, GATK, CaVEMan, ASCAT)
    // parallel run
    RUN_LUMPY_NORMAL()
    RUN_LUMPY_TUMOR()
    RUN_STATS_NORMAL()
    RUN_STATS_TUMOR()
    
    // using bed files
    RUN_GATK(PREPARE_BEDS.out.include_bed)
    RUN_CAVEMAN()
    RUN_ASCAT()
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
    val normal_merged_signal

    output:
    path "exclude_${params.normal}.bed", emit: exclude_bed
    path "include_${params.normal}.bed", emit: include_bed

    script:
    """
    module load SAMtools
    samtools view -H ${params.results}/${params.normal}/${params.normal}.cram | \\
    awk -F'[\\t:]' '\$1=="@SQ" && \$3 !~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ {print \$3"\\t1\\t"\$5}' > exclude_${params.normal}.bed

    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ {print \$1"\\t1\\t"\$2}' ${params.genome}.fai > include_${params.normal}.bed
    """
}

process RUN_LUMPY_NORMAL {
    conda 'lumpy-sv'

    script:
    """
    module load miniconda
    bash ${params.scripts}/sv.sh ${params.normal} ${params.results}/${params.normal}/${params.normal}.bam ${params.genome} ${params.annotate}/exclude_${params.normal}.bed
    """
}

process RUN_LUMPY_TUMOR {
    conda 'lumpy-sv'

    script:
    """
    module load miniconda
    bash ${params.scripts}/sv.sh ${params.tumor} ${params.results}/${params.tumor}/${params.tumor}.bam ${params.genome} ${params.annotate}/exclude_${params.normal}.bed
    """
}

process RUN_STATS_NORMAL {
    conda 'lumpy-sv'

    script:
    """
    module load miniconda
    export REF_PATH=${params.genome}
    bash ${params.scripts}/stats.sh ${params.normal} ${params.results}/${params.normal}/${params.normal}.bam
    """
}

process RUN_STATS_TUMOR {
    conda 'lumpy-sv'

    script:
    """
    module load miniconda
    export REF_PATH=${params.genome}
    bash ${params.scripts}/stats.sh ${params.tumor} ${params.results}/${params.tumor}/${params.tumor}.bam
    """
}

process RUN_GATK {
    input:
    path include_bed

    script:
    """
    module load GATK
    gatk Mutect2 -R ${params.genome} \\
        -I ${params.results}/${params.tumor}/${params.tumor}.bam \\
        -I ${params.results}/${params.normal}/${params.normal}.bam \\
        -L ${include_bed} \\
        -tumor ${params.tumor} \\
        -normal ${params.normal} \\
        -O ${params.results}/${params.tumor}/gatk-${params.tumor}.vcf
    """
}

process RUN_CAVEMAN {
    container 'cgpwgs.sif'

    script:
    """
    caveman.pl -o ${params.results}/${params.tumor}/caveman \\
        -r ${params.genome}.fai \\
        -tb ${params.results}/${params.tumor}/${params.tumor}.bam \\
        -nb ${params.results}/${params.normal}/${params.normal}.bam \\
        -ig ~/hg38-blacklist.v2.bed -tc ~/empty.txt -td 2 -nc ~/empty.txt -nd 2 \\
        -s Human -sa GRCh38 -b ~/empty.txt \\
        -in ~/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \\
        -st genome -u ~/empty_dir -t 96 -noflag
    """
}

process RUN_ASCAT {
    container 'cgpwgs.sif'

    script:
    """
    ascat.pl -o ${params.results}/${params.tumor}/ascat \\
        -t ${params.results}/${params.tumor}/${params.tumor}.bam \\
        -n ${params.results}/${params.normal}/${params.normal}.bam \\
        -r ${params.genome} -pr WGS -g XY -gc chrY \\
        -sg ~/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv -c 8
    """
}