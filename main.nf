// ==========================================
// main.nf
// ==========================================
nextflow.enable.dsl=2

// ==========================================
// WORKFLOW DEFINITION
// ==========================================
workflow {

    // 1. Parse lists of samples dynamically
    def tumor_ids = params.tumors.tokenize(',')
    def normal_id = params.normal

    // 2. Build explicit definitions for all target samples
    def tumor_definitions = tumor_ids.collect { id ->
        [ [id: id, type: 'tumor'], file("${params.results}/${id}/${id}.bam"), file("${params.results}/${id}/${id}.bam.bai") ]
    }
    def normal_definition = [ [id: normal_id, type: 'normal'], file("${params.results}/${normal_id}/${normal_id}.bam"), file("${params.results}/${normal_id}/${normal_id}.bam.bai") ]

    def all_sample_definitions = tumor_definitions + [normal_definition]

    // 3. Separate samples by checking physical file existence on disk
    def existing_samples = all_sample_definitions.findAll { meta, bam, bai -> bam.exists() && bai.exists() }
    def missing_samples  = all_sample_definitions.findAll { meta, bam, bai -> !bam.exists() || !bai.exists() }

    // 4. Initialize channels based on disk status
    ch_ready_samples = Channel.fromList(existing_samples)
    def ch_all_processed_samples

    if (missing_samples.size() > 0) {
        log.info "Found ${missing_samples.size()} samples missing their BAM files. Triggering PREPROCESS..."

        def missing_ids = missing_samples.collect { meta, bam, bai -> meta.id }

        // Fetch FASTQs only for the samples that actually need them
        ch_fastqs = Channel.fromFilePairs("${params.fastq_dir}/*{${missing_ids.join(',')}}*R{1,2}*.{fastq,fastq.gz}", checkIfExists: true)
            .map { name, reads ->
                def sample_id = missing_ids.find { name.contains(it) }
                def meta = missing_samples.find { m, b, bi -> m.id == sample_id }[0]
                return [ meta, reads ]
            }

        // Run Preprocess -> List Generation -> Merge for missing samples
        ch_preprocess_outputs = PREPROCESS(ch_fastqs)
        ch_grouped_bams = ch_preprocess_outputs.final_bam.groupTuple(by: 0)
        ch_bams_lists   = GENERATE_BAMS_LIST(ch_grouped_bams)
        ch_merged_bams  = MERGE_BAMS(ch_bams_lists)

        // Mix pre-existing static samples with newly generated ones
        ch_all_processed_samples = ch_ready_samples.mix(ch_merged_bams)
    } else {
        log.info "All BAM files found on disk. Skipping all preprocessing modules entirely."
        ch_all_processed_samples = ch_ready_samples
    }

    // 5. References Channels
    ch_genome_fa  = Channel.fromPath(params.genome).collect()
    ch_genome_fai = Channel.fromPath("${params.genome}.fai").collect()

    // 6. Split master channel back to specific sub-channels for analytical steps
    ch_tumor_samples = ch_all_processed_samples.filter { meta, bam, bai -> meta.type == 'tumor' }
    ch_normal_sample = ch_all_processed_samples.filter { meta, bam, bai -> meta.type == 'normal' }

    // 7. Dynamic Reference BED preparation (Triggers as soon as Normal BAM is resolved)
    PREPARE_BEDS(
        ch_normal_sample.map { meta, bam, bai -> bam },
        ch_normal_sample.map { meta, bam, bai -> bai },
        ch_genome_fai
    )

    // 8. Single-Sample Analyses (Parallel execution for all inputs via mix)
    ch_all_samples = ch_tumor_samples.mix(ch_normal_sample)

    RUN_LUMPY(ch_all_samples, PREPARE_BEDS.out.exclude_bed.collect())
    RUN_STATS(ch_all_samples)
    RUN_AMPLICONARCHITECT(ch_tumor_samples)

    // 9. Somatic Paired Analyses (Every tumor automatically mapped against the normal sample)
    ch_paired_somatic = ch_tumor_samples.combine(ch_normal_sample)

    RUN_GATK(ch_paired_somatic, PREPARE_BEDS.out.include_bed.collect())
    RUN_CAVEMAN(ch_paired_somatic, PREPARE_BEDS.out.include_bed.collect(), ch_genome_fa, PREPARE_BEDS.out.filtered_fai.collect())
    RUN_ASCAT(ch_paired_somatic)
}

// ==========================================
// GENERIC MODULES & PROCESSES
// ==========================================

process PREPROCESS {
    tag "${meta.id}"
    publishDir "${params.results}/${meta.id}/preprocess", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    // 1. The clean output needed for downstream tasks
    tuple val(meta), path("${meta.id}/*.dedup.bam"), emit: final_bam

    // 2. Catch-all for intermediate files (logs, sam, unsorted bams, etc.)
    // Using "*" captures every non-hidden file created in the work directory
    path("*"), emit: all_intermediates

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

    output:
    path "lumpy"

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

    output:
    path "stats"

    script:
    """
    module load miniconda
    export REF_PATH=${params.genome}
    bash ${params.scripts}/stats.sh ${meta.id} ${bam} 2>&1
    """
}

process RUN_GATK {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai)
    path include_bed

    output:
    path "gatk"

    script:
    """
    module load GATK
    mkdir -p gatk
    gatk Mutect2 -R ${params.genome} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -L ${include_bed} \\
        -tumor ${tumor_meta.id} \\
        -normal ${normal_meta.id} \\
        --pair-hmm-implementation LOGLESS_CACHING \\
        -O "gatk/${tumor_meta.id}.vcf" 2>&1
    """
}

process RUN_CAVEMAN {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai)
    path include_bed
    path genome_fasta
    path filtered_fai

    output:
    path "caveman"

    script:
    """
    ln -s ${genome_fasta} local_genome.fa
    ln -s ${filtered_fai} local_genome.fa.fai

    mkdir -p caveman
    caveman.pl -o caveman \\
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
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai)

    output:
    path "ascat"

    script:
    """
    mkdir -p ascat
    ascat.pl -o ascat \\
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
