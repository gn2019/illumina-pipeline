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

    ch_lumpy = RUN_LUMPY(ch_all_samples, PREPARE_BEDS.out.exclude_bed.collect())
    RUN_STATS(ch_all_samples)
    RUN_AMPLICONARCHITECT(ch_tumor_samples)
    ch_coverage = COVERAGE(ch_tumor_samples.map { meta, bam, bai -> meta })

    // 9. Somatic Paired Analyses
    ch_paired_somatic = ch_tumor_samples.combine(ch_normal_sample)

    // --- SCATTER PREPARATION (GATK only - Mutect2 benefits from BED-region chunking) ---
    SPLIT_BED_INTO_CHUNKS(PREPARE_BEDS.out.include_bed.collect())
    chunk_beds_ch = SPLIT_BED_INTO_CHUNKS.out.flatten()

    ch_scattered_somatic = ch_paired_somatic.combine(chunk_beds_ch)

    // --- RUN SCATTERED TOOLS ---
    RUN_GATK(ch_scattered_somatic)

    // --- GATHER RESULTS ---
    // merge GATK results
    gatk_per_sample_ch = RUN_GATK.out.groupTuple(by: 0)
    MERGE_GATK_VCFS(gatk_per_sample_ch)

    // ASCAT now has to run BEFORE CaVEMan: CaVEMan's copy-number input
    // (-tc/-nc) and normal contamination (-k) are derived directly from
    // ASCAT's per-tumor copynumber.caveman.csv and samplestatistics.txt.
    RUN_ASCAT(ch_paired_somatic, DOWNLOAD_REFS.out.ascat_gc_correction.collect())
    ASCAT_TO_CAVEMAN(RUN_ASCAT.out.ascat_out)

    // Build the caveman.pl copy-number flags ONCE per tumor/normal pair here,
    // so every caveman.pl step downstream (setup, split, mstep, estep, ...)
    // is invoked with identical, consistent flags. Default is to follow the
    // ASCAT->CaVEMan flow (params.follow_ascat_caveman_flow = true); set it
    // to false to always use the flat -td/-nd defaults instead. Even when
    // set to follow, if ASCAT didn't produce usable cn.bed files or a normal
    // contamination value for a sample, this still falls back to -td/-nd
    // rather than failing the pair.
    ch_cn_args = ASCAT_TO_CAVEMAN.out.cn_data.map { tumor_meta, normal_meta, tumor_cn_bed, normal_cn_bed, normal_contamination ->
        def has_cn_data = params.follow_ascat_caveman_flow &&
            tumor_cn_bed.size() > 0 && normal_cn_bed.size() > 0 && normal_contamination?.trim()
        def cn_args = has_cn_data
            ? "-tc tumor_cn.bed -nc normal_cn.bed -k ${normal_contamination.trim()}"
            : "-td 5 -nd 2"
        [ tumor_meta.id, tumor_meta, tumor_cn_bed, normal_cn_bed, cn_args ]
    }

    // CAVEMAN is NOT BED-chunked here. caveman.pl does its own internal
    // chromosome-level split/mstep/estep parallelization, so we drive that
    // explicitly as one Nextflow/LSF task per split index (mirrors the
    // manual bsub job-array script) instead of letting a single caveman.pl
    // invocation loop over everything internally with -t.
    ch_caveman_done = CAVEMAN(
        ch_paired_somatic,
        ch_genome_fa,
        PREPARE_BEDS.out.filtered_fai.collect(),
        DOWNLOAD_REFS.out.caveman_blacklist.collect(),
        DOWNLOAD_REFS.out.caveman_indels.collect(),
        DOWNLOAD_REFS.out.caveman_indels_tbi.collect(),
        ch_cn_args
    )

    // --- FINAL VISUALIZATION ARCHIVE ---
    // There's exactly one normal sample for the whole run (params.normal),
    // so it's referenced by id directly below rather than carried through
    // every channel as its own value.
    ch_tumor_lumpy_dir = ch_lumpy.lumpy_out
        .filter { meta, dir -> meta.type == 'tumor' }
        .map { meta, dir -> [ meta.id, dir ] }

    ch_normal_lumpy_dir = ch_lumpy.lumpy_out
        .filter { meta, dir -> meta.type == 'normal' }
        .map { meta, dir -> dir }
        .collect()

    ch_coverage_bg = ch_coverage.coverage
        .map { meta, bg -> [ meta.id, bg ] }

    ch_ascat_dir = RUN_ASCAT.out.ascat_out
        .map { tumor_meta, normal_meta, dir -> [ tumor_meta.id, dir ] }

    ch_caveman_flag = ch_caveman_done
        .map { tumor_meta, flag -> [ tumor_meta.id, flag ] }

    ch_archive_in = ch_paired_somatic
        .map { tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai -> [ tumor_meta.id, tumor_meta ] }
        .join(ch_tumor_lumpy_dir)
        .join(ch_coverage_bg)
        .join(ch_ascat_dir)
        .join(ch_caveman_flag)
        .combine(ch_normal_lumpy_dir)
        .map { id, tumor_meta, tumor_lumpy_dir, coverage_bg, ascat_dir, caveman_flag, normal_lumpy_dir ->
            [ tumor_meta, tumor_lumpy_dir, normal_lumpy_dir, coverage_bg, ascat_dir, caveman_flag ]
        }

    ARCHIVE(ch_archive_in)
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
    publishDir "${params.results}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path exclude_bed

    output:
    tuple val(meta), path("lumpy"), emit: lumpy_out

    script:
    """
    bash ${params.scripts}/sv.sh ${meta.id} ${bam} ${params.genome} ${exclude_bed} 2>&1
    """
}

// Locates the per-lane *_sorted.dedup.bw coverage tracks for a tumor sample
// (they're written into the fastq directory tree by preprocess.sh, next to
// the source fastqs, not into the Nextflow work dir) and merges them into a
// single sorted bedGraph.
process COVERAGE {
    tag "${tumor_meta.id}_coverage"
    conda "${params.bw_env}"
    publishDir "${params.results}/${tumor_meta.id}/coverage", mode: 'copy'

    input:
    val tumor_meta

    output:
    tuple val(tumor_meta), path("coverage.sorted.bedGraph"), emit: coverage

    script:
    """
    find "${params.results}/${tumor_meta.id}" -type f -name "*_sorted.dedup.bw" | sort > bw_list.txt

    bigWigMerge -inList bw_list.txt coverage.bedGraph

    sort -k1,1 -k2,2n coverage.bedGraph > coverage.sorted.bedGraph
    """
}

process RUN_STATS {
    tag "${meta.id}_${meta.type}"
    conda "${params.lumpy_env}"
    publishDir "${params.results}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "stats"

    script:
    """
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
    ch_cn_args             // tuple(tumor_meta.id, tumor_meta, tumor_cn_bed, normal_cn_bed, cn_args) - see ASCAT_TO_CAVEMAN

    main:
    ch_pair_ctx = ch_pairs
        .map { tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
            def outdir = "${params.results}/${tumor_meta.id}/caveman"
            [ tumor_meta.id, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir ]
        }
        .join(ch_cn_args.map { id, tumor_meta, tumor_cn_bed, normal_cn_bed, cn_args -> [ id, tumor_cn_bed, normal_cn_bed, cn_args ] })
        .map { id, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir, tumor_cn_bed, normal_cn_bed, cn_args ->
            [ tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, outdir, tumor_cn_bed, normal_cn_bed, cn_args ]
        }

    // SETUP is the only step that touches the bam/reference/cn.bed channels
    // - it symlinks everything caveman needs into the shared outdir ONCE,
    // using fixed relative names, and carries the resolved cn_args string
    // forward. Every step after this only needs (tumor_meta, outdir,
    // cn_args[, idx]) - it cds into outdir and the relative names already
    // resolve, so there's no repeated readlink/relinking per task.
    CAVEMAN_SETUP(ch_pair_ctx, genome_fasta, filtered_fai, caveman_blacklist, caveman_indels, caveman_indels_tbi)
    // CAVEMAN_SETUP.out: tuple(tumor_meta, outdir, cn_args)

    // one task per chromosome in the fai (mirrors NCHROM in the manual script)
    // NOTE: filtered_fai is a .collect()-ed channel, so it emits a List
    // (even with a single file inside) rather than a bare Path - index in.
    ch_nchrom = filtered_fai.map { fai -> (fai instanceof List ? fai[0] : fai).readLines().findAll { it.trim() }.size() }

    ch_split_in = CAVEMAN_SETUP.out
        .combine(ch_nchrom)
        .flatMap { tumor_meta, outdir, cn_args, nchrom ->
            (1..nchrom).collect { idx -> [ tumor_meta, outdir, cn_args, idx ] }
        }

    CAVEMAN_SPLIT(ch_split_in)

    // wait for ALL split tasks of a given pair before concatenating
    ch_split_grouped = CAVEMAN_SPLIT.out
        .groupTuple(by: 0)
        .map { tm, od, ca -> [ tm, od[0], ca[0] ] }

    CAVEMAN_SPLIT_CONCAT(ch_split_grouped)

    // caveman's split step subdivides each chromosome further - read the
    // actual splitList it produced (dynamic, unknown until now) and fan out
    // one mstep task per chunk, capped at 200 concurrent via maxForks
    ch_mstep_in = CAVEMAN_SPLIT_CONCAT.out
        .flatMap { tumor_meta, outdir, cn_args ->
            def nchunks = file("${outdir}/tmpCaveman/splitList").readLines().findAll { it.trim() }.size()
            (1..nchunks).collect { idx -> [ tumor_meta, outdir, cn_args, idx ] }
        }

    CAVEMAN_MSTEP(ch_mstep_in)

    ch_mstep_grouped = CAVEMAN_MSTEP.out
        .groupTuple(by: 0)
        .map { tm, od, ca -> [ tm, od[0], ca[0] ] }

    CAVEMAN_MERGE(ch_mstep_grouped)

    ch_estep_in = CAVEMAN_MERGE.out
        .flatMap { tumor_meta, outdir, cn_args ->
            def nchunks = file("${outdir}/tmpCaveman/splitList").readLines().findAll { it.trim() }.size()
            (1..nchunks).collect { idx -> [ tumor_meta, outdir, cn_args, idx ] }
        }

    CAVEMAN_ESTEP(ch_estep_in)

    ch_estep_grouped = CAVEMAN_ESTEP.out
        .groupTuple(by: 0)
        .map { tm, od, ca -> [ tm, od[0], ca[0] ] }

    CAVEMAN_MERGE_RESULTS(ch_estep_grouped)
    CAVEMAN_ADD_IDS(CAVEMAN_MERGE_RESULTS.out)

    emit:
    CAVEMAN_ADD_IDS.out
}

process CAVEMAN_SETUP {
    tag "${tumor_meta.id}_setup"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(tumor_meta), path(tumor_bam), path(tumor_bai), val(normal_meta), path(normal_bam), path(normal_bai), val(outdir), path(tumor_cn_bed), path(normal_cn_bed), val(cn_args)
    path genome_fasta
    path filtered_fai
    path caveman_blacklist
    path caveman_indels
    path caveman_indels_tbi

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    mkdir -p ${outdir}

    # Symlink everything caveman needs into the shared outdir under FIXED
    # relative names, once. Every later step just cds into outdir and
    # references these names directly - no re-linking, no path resolution.
    # tumor_cn.bed/normal_cn.bed are always linked (even when cn_args falls
    # back to -td/-nd and doesn't reference them) so the outdir layout is
    # identical either way.
    ln -fs \$(readlink -f ${genome_fasta}) ${outdir}/local_genome.fa
    ln -fs \$(readlink -f ${filtered_fai}) ${outdir}/local_genome.fa.fai
    ln -fs \$(readlink -f ${tumor_bam}) ${outdir}/tumor.bam
    ln -fs \$(readlink -f ${tumor_bai}) ${outdir}/tumor.bam.bai
    ln -fs \$(readlink -f ${normal_bam}) ${outdir}/normal.bam
    ln -fs \$(readlink -f ${normal_bai}) ${outdir}/normal.bam.bai
    ln -fs \$(readlink -f ${tumor_cn_bed}) ${outdir}/tumor_cn.bed
    ln -fs \$(readlink -f ${normal_cn_bed}) ${outdir}/normal_cn.bed
    ln -fs \$(readlink -f ${caveman_blacklist}) ${outdir}/blacklist.bed
    ln -fs \$(readlink -f ${caveman_indels}) ${outdir}/indels.vcf.gz
    ln -fs \$(readlink -f ${caveman_indels_tbi}) ${outdir}/indels.vcf.gz.tbi
    touch ${outdir}/empty.txt
    mkdir -p ${outdir}/empty_dir

    # caveman.pl records the CWD it was invoked from as part of its setup
    # state and refuses to run again from a different CWD, so every step
    # must cd into this same fixed, shared outdir before invoking it.
    cd ${outdir}

    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
        -process setup -index 1 2>&1
    """
}

process CAVEMAN_SPLIT {
    tag "${tumor_meta.id}_split_${idx}"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(tumor_meta), val(outdir), val(cn_args), val(idx)

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
        -process split -index ${idx} 2>&1
    """
}

process CAVEMAN_SPLIT_CONCAT {
    tag "${tumor_meta.id}_split_concat"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
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
    tuple val(tumor_meta), val(outdir), val(cn_args), val(idx)

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
        -process mstep -index ${idx} 2>&1
    """
}

process CAVEMAN_MERGE {
    tag "${tumor_meta.id}_merge"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '48 GB'

    input:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
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
    tuple val(tumor_meta), val(outdir), val(cn_args), val(idx)

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
        -process estep -index ${idx} 2>&1
    """
}

process CAVEMAN_MERGE_RESULTS {
    tag "${tumor_meta.id}_merge_results"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '48 GB'

    input:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    output:
    tuple val(tumor_meta), val(outdir), val(cn_args)

    script:
    """
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
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
    tuple val(tumor_meta), val(outdir), val(cn_args)

    output:
    tuple val(tumor_meta), path("caveman_done.flag")

    script:
    """
    ORIG_DIR=\$PWD
    cd ${outdir}
    caveman.pl -o ${outdir} \\
        -r local_genome.fa.fai \\
        -tb tumor.bam -nb normal.bam \\
        -ig blacklist.bed ${cn_args} \\
        -s Human -sa GRCh38 -b empty.txt \\
        -in indels.vcf.gz \\
        -st genome -u empty_dir -noflag \\
        -e 3500000 \\
        -process add_ids -index 1 2>&1

    echo "done" > ${outdir}/caveman_done.flag
    ln -fs ${outdir}/caveman_done.flag \$ORIG_DIR/caveman_done.flag
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
    tuple val(tumor_meta), val(normal_meta), path("ascat"), emit: ascat_out

    script:
    // -pu/-pi are only passed when BOTH ascat_purity and ascat_ploidy are
    // set - otherwise ascat.pl estimates them itself.
    def purity_ploidy_args = (params.ascat_purity && params.ascat_ploidy) ? "-pu ${params.ascat_purity} -pi ${params.ascat_ploidy}" : ''
    """
    mkdir -p ascat
    ascat.pl -o ascat \\
        -t ${tumor_bam} \\
        -n ${normal_bam} \\
        -r ${params.genome} -pr WGS -g XY -gc chrY \\
        -rs ${params.species} \\
        -ra ${params.assembly} \\
        -sg ${ascat_gc_correction} -c 8 ${purity_ploidy_args} 2>&1
    """
}

// ASCAT's per-tumour copynumber.caveman.csv already carries both the tumour
// and normal copy-number segments (columns: chr,start,stop,normal_total_cn,
// normal_minor_cn,tumour_total_cn,tumour_minor_cn - 1-indexed after
// splitting on comma), and its samplestatistics.txt carries the normal
// contamination fraction CaVEMan expects via -k. This process turns those
// into the two BED files + contamination value CaVEMan needs.
process ASCAT_TO_CAVEMAN {
    tag "${tumor_meta.id}_ascat2caveman"
    container "${params.cgpwgs_sif}"
    cpus 1
    memory '1 GB'

    input:
    tuple val(tumor_meta), val(normal_meta), path(ascat_dir)

    output:
    tuple val(tumor_meta), val(normal_meta), path("${tumor_meta.id}.cn.bed"), path("${normal_meta.id}.cn.bed"), env(NORMAL_CONTAMINATION), emit: cn_data

    script:
    """
    set +e

    perl -ne '@F=(split q{,}, \$_)[1,2,3,6]; \$F[1]--; print join("\\t",@F)."\\n";' \\
        < ${ascat_dir}/${tumor_meta.id}.copynumber.caveman.csv > ${tumor_meta.id}.cn.bed

    perl -ne '@F=(split q{,}, \$_)[1,2,3,4]; \$F[1]--; print join("\\t",@F)."\\n";' \\
        < ${ascat_dir}/${tumor_meta.id}.copynumber.caveman.csv > ${normal_meta.id}.cn.bed

    NORMAL_CONTAMINATION=\$(awk '(\$1=="NormalContamination") {print \$2}' ${ascat_dir}/${tumor_meta.id}.samplestatistics.txt)

    # if ASCAT's outputs weren't there or were malformed, still land empty
    # files/an empty value here rather than failing the pair - ch_cn_args
    # in the main workflow falls back to flat -td/-nd when this happens.
    touch ${tumor_meta.id}.cn.bed ${normal_meta.id}.cn.bed
    true
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
    AmpliconSuite-pipeline.py -s ${meta.id} -t 16 --bam ${bam} --run_AA --run_AC 2>&1
    """
}

// Bundles the per-tumor visualization inputs into a single tar.gz. sample.vcf
// (the tumor's lumpy SV calls) is the only mandatory member - everything
// else is included if present and silently skipped otherwise. The caveman
// SNV vcf is read from its well-known published path (same convention the
// CAVEMAN subworkflow above uses for its shared outdir) rather than staged
// as a formal Nextflow input, since it's one of several optional members.
process ARCHIVE {
    tag "${tumor_meta.id}_vs_${params.normal}_archive"
    conda "${params.bw_env}"
    publishDir "${params.results}/${tumor_meta.id}", mode: 'copy'

    input:
    tuple val(tumor_meta), path(tumor_lumpy_dir, stageAs: 'tumor_lumpy'), path(normal_lumpy_dir, stageAs: 'normal_lumpy'), path(coverage_bg), path(ascat_dir), path(caveman_flag)

    output:
    path "${tumor_meta.id}_vs_${params.normal}.visualize.tar.gz"

    script:
    def caveman_vcf_gz = "${params.results}/${tumor_meta.id}/caveman/${tumor_meta.id}_vs_${params.normal}.muts.ids.vcf.gz"
    def archive_name    = "${tumor_meta.id}_vs_${params.normal}.visualize.tar.gz"
    """
    set -euo pipefail
    workdir=\$(mktemp -d)

    # sample.vcf - the tumor's lumpy SV calls. Required.
    if [ ! -s tumor_lumpy/lumpy_sv.vcf ]; then
        echo "ERROR: required file missing: tumor lumpy SV vcf (tumor_lumpy/lumpy_sv.vcf)" >&2
        exit 1
    fi
    cp tumor_lumpy/lumpy_sv.vcf \$workdir/sample.vcf

    # background.vcf - the normal's lumpy SV calls. Optional.
    if [ -s normal_lumpy/lumpy_sv.vcf ]; then
        cp normal_lumpy/lumpy_sv.vcf \$workdir/background.vcf
    fi

    # snv.vcf - caveman's id'd somatic SNVs, ungzipped. Optional.
    if [ -s "${caveman_vcf_gz}" ]; then
        gunzip -c "${caveman_vcf_gz}" > \$workdir/snv.vcf
    fi

    # coverage.bedGraph - from the COVERAGE step above. Optional.
    if [ -s ${coverage_bg} ]; then
        cp ${coverage_bg} \$workdir/coverage.bedGraph
    fi

    # baf.bedGraph - from ASCAT's per-tumor BAF bigwig. Optional.
    baf_bw="${ascat_dir}/${tumor_meta.id}.copynumber.baf.bw"
    if [ -s "\$baf_bw" ]; then
        bigWigToBedGraph "\$baf_bw" \$workdir/baf.bedGraph
    fi

    tar -czf ${archive_name} -C \$workdir .
    rm -rf \$workdir
    """
}
