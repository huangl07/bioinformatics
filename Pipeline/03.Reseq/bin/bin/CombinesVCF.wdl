workflow Combines {
  call combines
}

task combines {
  File VCFlist
  String workdir
  File RefFasta
  File Refindex
  File Refdict
  command {
    java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar \
        -T CombineVariants \
        -V ${VCFlist} \
        -R ${RefFasta} \
        -o ${workdir}/pop.final.vcf \
	--genotypemergeoption UNSORTED \
	-log $out/pop.merge.log
  }
  output {
    File rawVCF = "${workdir}/pop.final.vcf"
  }
}
