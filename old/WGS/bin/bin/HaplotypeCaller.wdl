
workflow HaplotypeCaller {
  call haplotypeCaller
}

task haplotypeCaller {
  String sampleName
  File inputBAM
  File BAMindex
  String workdir
  File RefFasta
  File Refindex
  File Refdict
  command {
    java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R ${RefFasta} \
        -I ${inputBAM} \
        -o ${workdir}/${sampleName}.gvcfs \
		-nct 8 \
		--genotyping_mode DISCOVERY \
		--emitRefConfidence GVCF \
		-stand_call_conf 30 \
		-variant_index_type LINEAR \
		-variant_index_parameter 128000 \
		-filterNoBases \
		-filterMBQ \
		-filterRNC \
		-log ${workdir}/${sampleName}.gvcf.log
  }
  output {
    File rawVCF = "${workdir}/${sampleName}.gvcfs"
  }
}
