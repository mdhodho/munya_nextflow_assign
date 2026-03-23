nextflow.enable.dsl=2

params.reads                 = "data/*_{R1,R2}.fastq.gz"
params.reference             = "ref/reference.fa"
params.adapters              = "ref/adapters.fa"
params.outdir                = "results"
params.trimmomatic_container = "${baseDir}/containers/trimmomatic.sif"
params.threads               = 4
params.minlen                = 32

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")

    script:
    """
    fastqc -t ${params.threads} ${reads[0]} ${reads[1]}
    """
}

process TRIMMOMATIC {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    container "${params.trimmomatic_container}"

    input:
    tuple val(sample_id), path(reads)
    path adapters

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz")

    script:
    """
    java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads ${params.threads} \
        -phred33 \
        ${reads[0]} ${reads[1]} \
        ${sample_id}_R1.trim.fastq.gz ${sample_id}_R1.unpaired.fastq.gz \
        ${sample_id}_R2.trim.fastq.gz ${sample_id}_R2.unpaired.fastq.gz \
        ILLUMINACLIP:${adapters}:2:30:10 \
        LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:${params.minlen}
    """
}

process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")

    script:
    """
    fastqc -t ${params.threads} ${r1} ${r2}
    """
}

process PREP_REF {
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path ref

    output:
    tuple path("reference.fa"),
          path("reference.fa.amb"),
          path("reference.fa.ann"),
          path("reference.fa.bwt"),
          path("reference.fa.pac"),
          path("reference.fa.sa"),
          path("reference.fa.fai")

    script:
    """
    bwa index reference.fa
    samtools faidx reference.fa
    """
}

process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/sam", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    tuple path(ref),
          path(amb),
          path(ann),
          path(bwt),
          path(pac),
          path(sa),
          path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bwa mem -t ${params.threads} ${ref} ${r1} ${r2} > ${sample_id}.sam
    """
}

process SAM_TO_SORTED_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools view -@ ${params.threads} -bS ${sam_file} > ${sample_id}.bam
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam ${sample_id}.bam
    samtools index ${sample_id}.sorted.bam
    """
}

process VARCALL {
    tag "$sample_id"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    tuple path(ref),
          path(amb),
          path(ann),
          path(bwt),
          path(pac),
          path(sa),
          path(fai)

    output:
    path("${sample_id}.vcf")

    script:
    """
    bcftools mpileup -Ou -f ${ref} ${bam} | \
    bcftools call -mv -Ov -o ${sample_id}.vcf
    """
}

workflow {
    read_pairs = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { sample_id, reads -> tuple(sample_id, reads) }

    adapters_ch = Channel.fromPath(params.adapters, checkIfExists: true)
    ref_ch      = Channel.fromPath(params.reference, checkIfExists: true)

    FASTQC_RAW(read_pairs)

    trimmed_reads = TRIMMOMATIC(read_pairs, adapters_ch)

    FASTQC_TRIMMED(trimmed_reads)

    indexed_ref = PREP_REF(ref_ch)

    sam_files = BWA_MEM(trimmed_reads, indexed_ref)

    bam_files = SAM_TO_SORTED_BAM(sam_files)

    VARCALL(bam_files, indexed_ref)
}
