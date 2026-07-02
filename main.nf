// ==========================================
// main.nf
// ==========================================
nextflow.enable.dsl=2

// ==========================================
// WORKFLOW DEFINITION
// ==========================================
workflow {

    // 1. Define Meta Maps (Identities for our samples)
    def meta_tumor  = [id: params.tumor,  type: 'tumor']
    def meta_normal = [id: params.normal, type: 'normal']

    // 2. Expected final output BAM files after MERGE
    def tumor_final_bam  = file("${params.results}/${params.tumor}/${params.tumor}.bam")
    def tumor_final_bai  = file("${params.results}/${params.tumor}/${params.tumor}.bam.bai")
    def normal_final_bam = file("${params.results}/${params.normal}/${params.normal}.bam")
    def normal_final_bai = file("${params.results}/${params.normal}/${params.normal}.bam.bai")

    def ch_tumor_sample
    def ch_normal_sample

    // 3. Smart Skipping Logic: Check if final BAMs already exist on disk
    if (tumor_final_bam.exists() && normal_final_bam.exists()) {
        log.info "========================================================"
        log.info "Merged BAMs found on disk. Skipping PREPROCESS and MERGE."
        log.info "========================================================"

        // Feed the existing files directly into the channels
        ch_tumor_samples = Channel
            .fromList( params.tumors.tokenize(',') )
            .map { tumor_id ->
                def bam = file("${params.results}/${tumor_id}/${tumor_id}.bam")
                def bai = file("${params.results}/${tumor_id}/${tumor_id}.bam.bai")
                return [ [id: tumor_id, type: 'tumor'], bam, bai ]
            }
        ch_normal_sample = Channel.of([ meta_normal, normal_final_bam, normal_final_bai ])

    } else {
        log.info "========================================================"
        log.info "Merged BAMs missing. Starting PREPROCESS from FASTQs..."
        log.info "========================================================"

        // A. Read FASTQs from params.fastq_dir dynamically based on sample names
        // Expecting files like: {sample_id}_R1.fastq.gz
        ch_fastqs = Channel.fromFilePairs("${params.fastq_dir}/*{${params.tumor},${params.normal}}*R{1,2}*.{fastq,fastq.gz}", checkIfExists: true)
            .map { name, reads ->
                // Determine if this pair belongs to tumor or normal based on the filename
                def meta = name.contains(params.tumor) ? meta_tumor : meta_normal
                return [ meta, reads ]
            }

        // B. Run generic PREPROCESS
        ch_preprocessed = PREPROCESS(ch_fastqs)

        // C. Group preprocessed BAMs by type (Tumor together, Normal together)
        ch_grouped_bams = ch_preprocessed.groupTuple(by: 0)

        // D. Generate lists of BAMs
        ch_bams_lists = GENERATE_BAMS_LIST(ch_grouped_bams)

        // E. Run generic MERGE
        ch_merged_bams = MERGE_BAMS(ch_bams_lists)

        // F. Split the resulting merged BAMs into specific Tumor and Normal channels
        // to feed the rest of the pipeline
        ch_tumor_sample = ch_merged_bams.filter { meta, bam, bai -> meta.type == 'tumor' }
        ch_normal_sample = ch_merged_bams.filter { meta, bam, bai -> meta.type == 'normal' }
    }

    // 4. References Channels
    ch_genome_fa  = Channel.fromPath(params.genome).collect()
    ch_genome_fai = Channel.fromPath("${params.genome}.fai").collect()

    // 5. Downstream Analysis
    PREPARE_BEDS(
        ch_normal_sample.map { meta, bam, bai -> bam },
        ch_normal_sample.map { meta, bam, bai -> bai },
        ch_genome_fai
    )

    // Mix tumor and normal to run them in parallel through single-sample tools
    ch_all_samples = ch_tumor_samples.mix(ch_normal_sample)

    RUN_LUMPY(ch_all_samples, PREPARE_BEDS.out.exclude_bed.collect())
    RUN_STATS(ch_all_samples)
    RUN_AMPLICONARCHITECT(ch_tumor_sample)

    ch_paired_for_somatic = ch_tumor_samples.combine(ch_normal_sample)

    // Paired Somatic Analyses
    RUN_GATK(
        ch_paired_for_somatic,
        PREPARE_BEDS.out.include_bed.collect()
    )

    RUN_CAVEMAN(
        ch_paired_for_somatic,
        PREPARE_BEDS.out.include_bed.collect(),
        ch_genome_fa,
        PREPARE_BEDS.out.filtered_fai.collect()
    )

    RUN_ASCAT(
        ch_paired_for_somatic
    )
}

// ==========================================
// GENERIC MODULES & PROCESSES
// ==========================================

process PREPROCESS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.dedup.bam")

    script:
    // Modified to pass the fastq files dynamically to your shell script
    """
    bash ${params.scripts}/preprocess.sh ${meta.id} noERX ${reads[0]} ${reads[1]} ${params.genome} 2>&1
    """
}

process GENERATE_BAMS_LIST {
    tag "${meta.type}"

    input:
    tuple val(meta), path(bams)

    output:
    // Passes forward the metadata, the text file with the list, and the actual bams
    tuple val(meta), path("bams_${meta.type}.txt"), path(bams)

    script:
    """
    # Nextflow brings the bams into the working directory, so we just list them
    ls *.bam > bams_${meta.type}.txt
    """
}

process MERGE_BAMS {
    tag "${meta.type}"
    publishDir "${params.results}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bams_list), path(bams)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")

    script:
    """
    bash ${params.scripts}/merge.sh ${meta.id} ${bams_list} 2>&1
    """
}

process PREPARE_BEDS {
    tag "Prepare_Refs"

    input:
    path normal_bam
    path normal_bai
    path genome_fai

    output:
    path "exclude.bed", emit: exclude_bed
    path "include.bed", emit: include_bed
    path "filtered_genome.fai", emit: filtered_fai

    script:
    """
    module load SAMtools

    samtools view -H ${normal_bam} | \\
    awk -F'[\\t:]' '\$1=="@SQ" && \$3 !~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ {print \$3"\\t1\\t"\$5}' > exclude.bed

    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/ {print \$1"\\t1\\t"\$2}' ${genome_fai} > include.bed

    awk '\$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)\$/' ${genome_fai} > filtered_genome.fai
    """
}

process RUN_LUMPY {
    tag "${meta.id}_${meta.type}"
    conda "${params.lumpy_env}"
    publishDir "${params.results}/${meta.id}/lumpy", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path exclude_bed

    script:
    """
    module load miniconda
    bash ${params.scripts}/sv.sh ${meta.id} ${bam} ${params.genome} ${exclude_bed} 2>&1
    """
}

process RUN_STATS {
    tag "${meta.id}_${meta.type}"
    conda "${params.lumpy_env}"
    publishDir "${params.results}/${meta.id}/stats", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    script:
    """
    module load miniconda
    export REF_PATH=${params.genome}
    bash ${params.scripts}/stats.sh ${meta.id} ${bam} 2>&1
    """
}

process RUN_GATK {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    publishDir "${params.results}/${tumor_meta.id}/gatk", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai)
    tuple val(normal_meta), path(normal_bam), path(normal_bai)
    path include_bed

    output:
    path "gatk-${tumor_meta.id}.vcf"

    script:
    """
    module load GATK
    gatk Mutect2 -R ${params.genome} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -L ${include_bed} \\
        -tumor ${tumor_meta.id} \\
        -normal ${normal_meta.id} \\
        --pair-hmm-implementation LOGLESS_CACHING \\
        -O "gatk-${tumor_meta.id}.vcf" 2>&1
    """
}

process RUN_CAVEMAN {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${tumor_meta.id}/caveman", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai)
    tuple val(normal_meta), path(normal_bam), path(normal_bai)
    path include_bed
    path genome_fasta
    path filtered_fai

    output:
    path "caveman_results"

    script:
    """
    ln -s ${genome_fasta} local_genome.fa
    ln -s ${filtered_fai} local_genome.fa.fai

    mkdir -p caveman_results
    caveman.pl -o caveman_results \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} \\
        -nb ${normal_bam} \\
        -ig ${params.caveman_blacklist} -tc ~/empty.txt -td 2 -nc ~/empty.txt -nd 2 \\
        -s Human -sa GRCh38 -b ~/empty.txt \\
        -in ${params.caveman_indels} \\
        -st genome -u ~/empty_dir -t 96 -noflag 2>&1
    """
}

process RUN_ASCAT {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${tumor_meta.id}/ascat", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai)
    tuple val(normal_meta), path(normal_bam), path(normal_bai)

    output:
    path "ascat_results"

    script:
    """
    mkdir -p ascat_results
    ascat.pl -o ascat_results \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -r ${params.genome} -pr WGS -g XY -gc chrY \\
        -rs ${params.species} \\
        -ra ${params.assembly} \\
        -sg ${params.ascat_gc_correction} -c 8 2>&1
    """
}

process RUN_AMPLICONARCHITECT {
    tag "${meta.id}"
    conda "${params.ampsuite_env}"
    publishDir "${params.results}/${meta.id}/AmpliconSuite", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "*"

    script:
    """
    module load miniconda
    AmpliconSuite-pipeline.py -s ${meta.id} -t 16 --bam ${bam} --run_AA --run_AC 2>&1
    """
}
