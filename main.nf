// According GATK best practices for ACGT:
// https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
// plus VEP annotation
process HAPLOTYPECALLER_GVCF {

    publishDir "${params.pubdir}/results/HC", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    cpus 2

    input:
      tuple val(sample), val(bam), val(bai)

    output:
      path '*'

    script:
    """
    gatk --java-options "-Xmx16g -Xms16g" HaplotypeCaller  \
    -R $params.ref \
    -I $bam \
    -L chr10:18000-45500 \
    -O ${sample}.g.vcf.gz \
    -ERC GVCF

    """    
}

process COMBINE_GVCF {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport

    publishDir "${params.pubdir}/results/combine_vcf", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    cpus 32
    memory 256.GB

    input:
      path vcf_list

    output:
      path './acgt_database'

    script:

      """
      for i in `ls ${params.pubdir}/results/HC/*gz | xargs -n1 basename`
      do
        awk -v id="\${i%.g.vcf.gz}" -v vcf="\$i" -v dir="${params.pubdir}" 'BEGIN{print id, dir "results/HC/" vcf}' | sed 's/ /\t/g' >> samples.names
      done

      gatk --java-options "-Xmx4g -Xms4g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ./acgt_database \
        --batch-size 100 \
        --sample-name-map samples.names \
        --L $params.interval \
        --tmp-dir=./ \
        --reader-threads 32
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
      gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -R $params.ref \
      -V gendb://./results/combine_vcf/acgt_database \
      -O ACGT_joint.vcf.gz
      """
}

process VAR_RECALL {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator

    publishDir "${launchDir}/results/filtered/", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    
    input:
      file joint_gvcf

    output:
      path '*'

    script:
      """
      


      """
}

process APPLY_RECALL {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR

    publishDir "${launchDir}/results/filtered/", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    
    input:
      file recall_gvcf

    output:
      path '*'

    script:
      """
      


      """
}

workflow {
    // create channel
    input_ch = Channel.empty()
    tsv = file(params.input)
    input_ch = extractFastq(tsv)
    //input_ch = Channel.fromPath('/mnt2/shared/ACGT/join_VC/input/HaplotypeCaller_ACGT00{1,2}.vcf.gz').view()

    // execute workflow
    vcf = HAPLOTYPECALLER_GVCF(input_ch)
    batch_vcf = vcf.collect()
    comb_gvcf = COMBINE_GVCF(batch_vcf)
    JOIN_GVCF(comb_gvcf)
    //var_recall_model = VAR_RECALL(join_vcf)
    //APPLY_RECALL(var_recall_model)
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
