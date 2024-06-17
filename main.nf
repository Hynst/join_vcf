// According GATK best practices for ACGT:
// https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
// plus VEP annotation

process HAPLOTYPECALLER_GVCF {

    publishDir "${params.pubdir}/results/HC", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    cpus 2
    memory 22.GB

    input:
      tuple val(sample), val(bam), val(bai)

    output:
      path '*'

    script:
    """
    gatk --java-options "-Xmx22g -Xms22g -XX:ParallelGCThreads=2" HaplotypeCaller  \
    -R $params.ref \
    -I $bam \
    -L $params.interval \
    -O ${sample}.g.vcf.gz \
    --min-base-quality-score 20 \
    -ERC GVCF
    """    
}

process COMBINE_GVCF {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport

    publishDir "${params.pubdir}/results/combine_vcf/", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    cpus 2
    memory 16.GB

    input:
      tuple val(reg), val(bed)

    output:
      tuple val(reg), file("./${reg}_acgt_database")

    script:

      """
      for i in `ls ${params.pubdir}/results/HC/*gz | xargs -n1 basename`
      do
        awk -v id="\${i%.g.vcf.gz}" -v vcf="\$i" -v dir="${params.pubdir}" 'BEGIN{print id, dir "results/HC/" vcf}' | sed 's/ /\t/g' >> samples.names
      done

      gatk --java-options "-Xmx4g -Xms4g -XX:ParallelGCThreads=2" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ./${reg}_acgt_database \
        --batch-size 80 \
        --genomicsdb-shared-posixfs-optimizations \
        --merge-input-intervals \
        --sample-name-map samples.names \
        --L ${bed} \
        --tmp-dir . \
        --reader-threads 2
      """
}

process JOIN_GVCF {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs

    publishDir "${params.pubdir}/results/join_vcf/", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    
    input:
      tuple val(reg), file(db)

    output:
      path "${reg}_ACGT_joint*"

    script:
      """
      gatk --java-options "-Xms4g -Xmx4g -XX:ParallelGCThreads=2" GenotypeGVCFs \
      -R $params.ref \
      -V gendb://${db} \
      -O ${reg}_ACGT_joint.vcf.gz
      """
}

process MERGE_VCFS {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036713331-MergeVcfs-Picard

    publishDir "${params.pubdir}/results/merged_vcf/", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'

    input:
      file joint_vcf

    output:
      path 'ACGT_joint_merged.vcf*'

    script:
      """
      touch jointVCFS_files.list
      for file in `ls ${params.pubdir}/results/join_vcf/*vcf.gz`
      do
        gatk SortVcf -I \${file} -O \${file%.vcf.gz}_sorted.vcf.gz -SD ${params.dict}
        echo \${file%.vcf.gz}_sorted.vcf.gz >> jointVCFS_files.list
      done

      gatk MergeVcfs \
      -I jointVCFS_files.list \
      -O ACGT_joint_merged.vcf.gz \
      -D ${params.dict}
      """
}

process VAR_RECALL {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator

    publishDir "${params.pubdir}/results/model/", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    
    input:
      file merged_vcf

    output:
      path 'ACGT_*'

    script:
      """
      gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
        -R $params.ref \
        -V ACGT_joint_merged.vcf.gz \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.annot}/hapmap_3.3.hg38.vcf.gz \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.annot}/1000G_omni2.5.hg38.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.annot}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.annot}/dbsnp_146.hg38.vcf.gz \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        --max-gaussians 4 \
        -mode SNP \
        -O ACGT_variants.recal \
        --tranches-file ACGT_variants.tranches \
        --rscript-file ACGT_variants.plots.R
      """
}

process APPLY_RECALL {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR

    publishDir "${params.pubdir}/results/recal/", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    
    input:
      file recal

    output:
      path 'ACGT_variants_recall.vcf.*'

    script:
      """
       gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" ApplyVQSR \
        -R $params.ref \
        -V ${params.pubdir}/results/merged_vcf/ACGT_joint_merged.vcf.gz \
        -O ACGT_variants_recall.vcf.gz \
        --tranches-file ACGT_variants.tranches \
        --recal-file ACGT_variants.recal \
        -mode SNP
      """
}

workflow {
    // create channel
    //input_ch = Channel.empty()
    //tsv = file(params.input)
    //input_ch = extractFastq(tsv)

    // execute workflow
    //vcf = HAPLOTYPECALLER_GVCF(input_ch)
    //batch_vcf = vcf.collect()
    bed_ch = Channel.empty()
    int_tsv = file(params.interval)
    bed_ch = extractBeds(int_tsv)
    //com_ch = Channel.empty()
    //com_tsv = file(params.interval)
    //com_ch = extractComDB(com_tsv)


    comb_gvcf = COMBINE_GVCF(bed_ch)
    genotypegvcf = JOIN_GVCF(comb_gvcf)
    merged_vcfs = MERGE_VCFS(genotypegvcf.collect())
    var_recall_model = VAR_RECALL(merged_vcfs)
    APPLY_RECALL(var_recall_model)
}

def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->           
            def sample    = row[0]
            def bam       = returnFile(row[1])
            def bai       = returnFile(row[2])

            [sample, bam, bai]
        }
}

def extractBeds(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->           
            def reg    = row[0]
            def bed       = returnFile(row[1])

            [reg, bed]
        }
}


def extractComDB(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def db    = row[0]
            def reg   = row[1]

            [db, reg]
        }
}
