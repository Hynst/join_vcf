// According GATK best practices for ACGT:
// https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
// plus VEP annotation

process HAPLOTYPECALLER_GVCF {

    publishDir "${params.pubdir}/results/HC", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
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

    publishDir "${params.pubdir}/results/combine_vcf", mode: 'copy'
    container 'broadinstitute/gatk:4.2.3.0'
    cpus 16
    memory 10.GB

    //input:
    //  path vcf_list

    output:
      path './acgt_database'

    script:

      """
      for i in `ls ${params.pubdir}/results/HC/*gz | xargs -n1 basename`
      do
        awk -v id="\${i%.g.vcf.gz}" -v vcf="\$i" -v dir="${params.pubdir}" 'BEGIN{print id, dir "results/HC/" vcf}' | sed 's/ /\t/g' >> samples.names
      done

      gatk --java-options "-Xmx4g -Xms4g -XX:ParallelGCThreads=2" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ./acgt_database \
        --batch-size 3 \
        --genomicsdb-shared-posixfs-optimizations \
        --merge-input-intervals \
        --sample-name-map samples.names \
        --L $params.interval \
        --tmp-dir . \
        --reader-threads 8
      """
}

process JOIN_GVCF {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs

    publishDir "${params.pubdir}/results/join_vcf/", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    
    input:
      val combined_vcf

    output:
      path '*'

    script:
      """
      gatk --java-options "-Xms4g -Xmx4g -XX:ParallelGCThreads=2" GenotypeGVCFs \
      -R $params.ref \
      -V gendb://${params.pubdir}/results/combine_vcf/acgt_database \
      -O ACGT_joint.vcf.gz
      """
}

process VAR_RECALL {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator

    publishDir "${params.pubdir}/results/filtered/", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    
    input:
      file genotype_gvcf

    output:
      path '*'

    script:
      """
      gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
        -R $params.ref \
        -V ACGT_joint.vcf.gz \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${params.annot}/hapmap_3.3.hg38.vcf.gz \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${params.annot}/1000G_omni2.5.hg38.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${params.annot}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.annot}/dbsnp_146.hg38.vcf.gz \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O ACGT_variants.recal \
        --tranches-file ACGT_variants.tranches \
        --rscript-file ACGT_variants.plots.R
      """
}

process APPLY_RECALL {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR

    publishDir "${params.pubdir}/results/filtered/", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    
    input:
      file recall_gvcf

    output:
      path '*'

    script:
      """
       gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" ApplyVQSR \
        -R $params.ref \
        -V ${params.pubdir}/results/join_vcf/ACGT_joint.vcf.gz \
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
    //input_ch = Channel.fromPath('/mnt2/shared/ACGT/join_VC/input/HaplotypeCaller_ACGT00{1,2}.vcf.gz').view()

    // execute workflow
    //vcf = HAPLOTYPECALLER_GVCF(input_ch)
    //batch_vcf = vcf.collect()
    comb_gvcf = COMBINE_GVCF()
    genotypegvcf = JOIN_GVCF(comb_gvcf)
    var_recall_model = VAR_RECALL(genotypegvcf)
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
