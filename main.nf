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

    DOWNLOAD_REFS()

    if (missing_samples.size() > 0) {
        log.info "Found ${missing_samples.size()} samples missing their BAM files. Triggering PREPROCESS..."

        def missing_ids = missing_samples.collect { meta, bam, bai -> meta.id }

        // Fetch FASTQs only for the samples that actually need them
        ch_fastqs = Channel.fromFilePairs("${params.fastq_dir}/*{${missing_ids.join(',')}}*{1,2}.{fastq,fastq.gz}", checkIfExists: true)
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

    // 9. Somatic Paired Analyses
    ch_paired_somatic = ch_tumor_samples.combine(ch_normal_sample)

    // --- SCATTER PREPARATION (GATK only - Mutect2 benefits from BED-region chunking) ---
    SPLIT_BED_INTO_CHUNKS(PREPARE_BEDS.out.include_bed.collect())
    chunk_beds_ch = SPLIT_BED_INTO_CHUNKS.out.flatten()

    ch_scattered_somatic = ch_paired_somatic.combine(chunk_beds_ch)

    // --- RUN SCATTERED TOOLS ---
    RUN_GATK(ch_scattered_somatic)

    // CAVEMAN is NOT BED-chunked here. caveman.pl does its own internal
    // chromosome-level split/mstep/estep parallelization, so we drive that
    // explicitly as one Nextflow/LSF task per split index (mirrors the
    // manual bsub job-array script) instead of letting a single caveman.pl
    // invocation loop over everything internally with -t.
    CAVEMAN(
        ch_paired_somatic,
        ch_genome_fa,
        PREPARE_BEDS.out.filtered_fai.collect(),
        DOWNLOAD_REFS.out.caveman_blacklist.collect(),
        DOWNLOAD_REFS.out.caveman_indels.collect(),
        DOWNLOAD_REFS.out.caveman_indels_tbi.collect()
    )

    // --- GATHER RESULTS ---
    // merge GATK results
    gatk_per_sample_ch = RUN_GATK.out.groupTuple(by: 0)
    MERGE_GATK_VCFS(gatk_per_sample_ch)

    RUN_ASCAT(ch_paired_somatic, DOWNLOAD_REFS.out.ascat_gc_correction.collect())
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
    tuple val(meta), path("${meta.id}/*.dedup.bam"), emit: final_bam
    path("*"), emit: all_intermediates

    script:
    """
    bash ${params.scripts}/preprocess.sh ${meta.id} noERX ${reads[0]} ${reads[1]} ${params.genome} 2>&1
    """
}

process DOWNLOAD_REFS {
    output:
    path "caveman_blacklist.bed", emit: caveman_blacklist
    path "caveman_indels.vcf.gz", emit: caveman_indels
    path "caveman_indels.vcf.gz.tbi", emit: caveman_indels_tbi
    path "ascat_gc.txt", emit: ascat_gc_correction

    script:
    """
    module load BEDTools

    bash ${params.scripts}/download_refs.sh hg38

    ln -fs ${params.caveman_blacklist} caveman_blacklist.bed
    ln -fs ${params.caveman_indels} caveman_indels.vcf.gz
    ln -fs ${params.caveman_indels}.tbi caveman_indels.vcf.gz.tbi
    ln -fs ${params.ascat_gc_correction} ascat_gc.txt
    """
}

process GENERATE_BAMS_LIST {
    tag "${meta.type}"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("bams_${meta.type}.txt"), path(bams)

    script:
    """
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

// --- NEW SCATTER PROCESS ---
process SPLIT_BED_INTO_CHUNKS {
    tag "split_bed"

    input:
    path main_bed

    output:
    path "*.bed"

    script:
    """
	awk -v chunk=50000000 '{
        for (i=0; i<\$3; i+=chunk) {
            end = (i+chunk > \$3) ? \$3 : i+chunk
            print \$1"\\t"i"\\t"end > \$1"_"i".bed"
        }
    }' ${main_bed}
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
    tag "${tumor_meta.id}_vs_${normal_meta.id}_${chunk_bed.baseName}"

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), path(chunk_bed)

    output:
    tuple val(tumor_meta), path("gatk-${tumor_meta.id}-${chunk_bed.baseName}.vcf")

    script:
    """
    module load GATK
    gatk Mutect2 -R ${params.genome} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -L ${chunk_bed} \\
        -tumor ${tumor_meta.id} \\
        -normal ${normal_meta.id} \\
        --native-pair-hmm-threads 8 \\
        --pair-hmm-implementation LOGLESS_CACHING \\
        -O "gatk-${tumor_meta.id}-${chunk_bed.baseName}.vcf" 2>&1
    """
}

// ==========================================
// CAVEMAN SUBWORKFLOW
// ==========================================
// caveman.pl orchestrates setup -> split -> split_concat -> mstep -> merge ->
// estep -> merge_results -> add_ids, and every one of those steps reads and
// writes into the SAME shared output directory tree (it manages its own
// state there - it is not a clean isolated in/out tool). So instead of
// scattering by BED region and running caveman.pl end-to-end inside a single
// task (which forced its own internal split/mstep/estep to run serially,
// using only local -t threads), each step below is its own Nextflow process,
// fanned out one task per split index - exactly like the manual bsub
// job-array script. Every step targets a FIXED absolute outdir per
// tumor/normal pair so state persists across tasks; only lightweight
// context/index values flow through the channels.
workflow CAVEMAN {
    take:
    ch_pairs              // tuple(tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai)
    genome_fasta
    filtered_fai
    caveman_blacklist
    caveman_indels
    caveman_indels_tbi

    main:
    ch_pair_ctx = ch_pairs.map { tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
        def outdir = "${params.results}/${tumor_meta.id}/caveman"
        [ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir ]
    }

    CAVEMAN_SETUP(ch_pair_ctx, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    // one task per chromosome in the fai (mirrors NCHROM in the manual script)
    // NOTE: filtered_fai is a .collect()-ed channel, so it emits a List
    // (even with a single file inside) rather than a bare Path - index in.
    ch_nchrom = filtered_fai.map { fai -> (fai instanceof List ? fai[0] : fai).readLines().findAll { it.trim() }.size() }

    ch_split_in = CAVEMAN_SETUP.out
        .combine(ch_nchrom)
        .flatMap { tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir, nchrom ->
            (1..nchrom).collect { idx -> [ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir, idx ] }
        }

    CAVEMAN_SPLIT(ch_split_in, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    // wait for ALL split tasks of a given pair before concatenating
    ch_split_grouped = CAVEMAN_SPLIT.out
        .groupTuple(by: 0)
        .map { tm, tb, tbi, nm, nb, nbi, od -> [ tm, tb[0], tbi[0], nm[0], nb[0], nbi[0], od[0] ] }

    CAVEMAN_SPLIT_CONCAT(ch_split_grouped, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    // caveman's split step subdivides each chromosome further - read the
    // actual splitList it produced (dynamic, unknown until now) and fan out
    // one mstep task per chunk, capped at 200 concurrent via maxForks
    ch_mstep_in = CAVEMAN_SPLIT_CONCAT.out
        .flatMap { tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir ->
            def nchunks = file("${outdir}/tmpCaveman/splitList").readLines().findAll { it.trim() }.size()
            (1..nchunks).collect { idx -> [ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir, idx ] }
        }

    CAVEMAN_MSTEP(ch_mstep_in, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    ch_mstep_grouped = CAVEMAN_MSTEP.out
        .groupTuple(by: 0)
        .map { tm, tb, tbi, nm, nb, nbi, od -> [ tm, tb[0], tbi[0], nm[0], nb[0], nbi[0], od[0] ] }

    CAVEMAN_MERGE(ch_mstep_grouped, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    ch_estep_in = CAVEMAN_MERGE.out
        .flatMap { tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir ->
            def nchunks = file("${outdir}/tmpCaveman/splitList").readLines().findAll { it.trim() }.size()
            (1..nchunks).collect { idx -> [ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir, idx ] }
        }

    CAVEMAN_ESTEP(ch_estep_in, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    ch_estep_grouped = CAVEMAN_ESTEP.out
        .groupTuple(by: 0)
        .map { tm, tb, tbi, nm, nb, nbi, od -> [ tm, tb[0], tbi[0], nm[0], nb[0], nbi[0], od[0] ] }

    CAVEMAN_MERGE_RESULTS(ch_estep_grouped, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)
    CAVEMAN_ADD_IDS(CAVEMAN_MERGE_RESULTS.out, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)

    emit:
    CAVEMAN_ADD_IDS.out
}

process CAVEMAN_SETUP {
    tag "${tumor_meta.id}_setup"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    mkdir -p ${outdir}
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process setup -index 1 2>&1
    """
}

process CAVEMAN_SPLIT {
    tag "${tumor_meta.id}_split_${idx}"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir), val(idx)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process split -index ${idx} 2>&1
    """
}

process CAVEMAN_SPLIT_CONCAT {
    tag "${tumor_meta.id}_split_concat"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process split_concat -index 1 2>&1
    """
}

process CAVEMAN_MSTEP {
    tag "${tumor_meta.id}_mstep_${idx}"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '64 GB'
    maxForks 200

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir), val(idx)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process mstep -index ${idx} 2>&1
    """
}

process CAVEMAN_MERGE {
    tag "${tumor_meta.id}_merge"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '48 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process merge -index 1 2>&1
    """
}

process CAVEMAN_ESTEP {
    tag "${tumor_meta.id}_estep_${idx}"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '64 GB'
    maxForks 200

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir), val(idx)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process estep -index ${idx} 2>&1
    """
}

process CAVEMAN_MERGE_RESULTS {
    tag "${tumor_meta.id}_merge_results"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '48 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process merge_results -index 1 2>&1
    """
}

process CAVEMAN_ADD_IDS {
    tag "${tumor_meta.id}_add_ids"
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy', pattern: 'caveman_done.flag'
    cpus 1
    memory '48 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    path "caveman_done.flag"

    script:
    """
    touch empty.txt
    mkdir -p empty_dir
    ln -fs ${genome_fasta} local_genome.fa
    ln -fs ${filtered_fai} local_genome.fa.fai

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb ${tumor_bam} -nb ${normal_bam} \\
        -ig ${caveman_blacklist} -tc empty.txt -td 3 -nc empty.txt -nd 3 \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in ${caveman_indels} \\
        -st genome -u empty_dir -noflag \\
        -process add_ids -index 1 2>&1

    echo "done" > caveman_done.flag
    """
}

// --- NEW GATHER PROCESSES ---
process MERGE_GATK_VCFS {
    tag "${tumor_meta.id}"
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy'

    input:
    tuple val(tumor_meta), path(vcf_list)

    output:
    path "gatk/${tumor_meta.id}_merged.vcf"

    script:
    def input_list = vcf_list.collect { "-I ${it}" }.join(' ')
    """
    module load GATK
    mkdir -p gatk
    gatk MergeVcfs ${input_list} -O "gatk/${tumor_meta.id}_merged.vcf"
    """
}

process RUN_ASCAT {
    tag "${tumor_meta.id}_vs_${normal_meta.id}"
    container "${params.cgpwgs_sif}"
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai)
    path ascat_gc_correction

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
        -sg ${ascat_gc_correction} -c 8 2>&1
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
