// According GATK best practices:
// https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
// plus VEP annotation

process COMBINE_GVCF {
    // !!! CombineGVCFs will be probably replaced with:
    // https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport
    
    publishDir "${launchDir}/results/combine_vcf", mode: 'copy'
    container 'broadinstitute/gatk:4.1.3.0'
    cpus 32
    memory 256.GB
    
    input:
      file gvcfs_list

    output:
      path './acgt_database'

    script:
      """
   
       gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path ./acgt_database \
       --batch-size 100 \
       --sample-name-map $gvcfs_list \
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
    // create channels
    input_ch = Channel.fromPath(params.input)
    //input_ch = Channel.fromPath('/mnt2/shared/ACGT/join_VC/input/HaplotypeCaller_ACGT00{1,2}.vcf.gz').view()

    // execute workflow
    //COMBINE_GVCF(input_ch)
    comb_gvcf = COMBINE_GVCF(input_ch)
    JOIN_GVCF(comb_gvcf)
    //var_recall_model = VAR_RECALL(join_vcf)
    //APPLY_RECALL(var_recall_model)
}
