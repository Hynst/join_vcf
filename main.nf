// According GATK best practices for ACGT:
// https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
// plus VEP annotation
process HAPLOTYPECALLER_GVCF{

    publishDir "${launchDir}/results/HC", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    cpus 2

    input:
      tuple val(sample), val(bam), val(bai)

    output:
      path '*.g.vcf.gz'

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller  \
    -R $params.ref \
    -I $bam \
    -O ${sample}.g.vcf.gz \
    -ERC GVCF
    """    

}

process COMBINE_GVCF {
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport
    
    publishDir "${launchDir}/results/combine_vcf", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    cpus 32
    memory 256.GB
    
    input:
      path gvcfs

    output:
      path './acgt_database'

    script:
      """
      # create vcf files with sample id
      dir=`echo ${launchDir}/results/HC/`
      for i in $(ls ${launchDir}/results/HC)
      do
        awk -v id="${i%.g.vcf.gz}" -v vcf="$i" -v dir="$dir" 'BEGIN{print id, dir vcf}' | \
        sed 's/ /\t/g' > samples.names
      done

      # run gatk DB for GenotypeGVCFs
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

    publishDir "${launchDir}/results/join_vcf/", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    
    input:
      val combined_vcf

    output:
      path '*'

    script:
      """
      gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -R $params.ref \
      -V gendb://${launchDir}/results/combine_vcf/acgt_database \
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
    comb_gvcf = COMBINE_GVCF(vcf)
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