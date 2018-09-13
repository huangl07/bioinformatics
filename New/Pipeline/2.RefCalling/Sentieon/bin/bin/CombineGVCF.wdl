workflow CombinesGVCF {
  call CombineVCF
}


task CombineVCF {
  File inputVCFs
  String workdir
  String Output
  File RefFasta
  File Refindex
  File Refdict
  command {
    java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar \
        -T  CombineGVCFs \
        -R ${RefFasta} \
        -V ${inputVCFs} \
        -o ${workdir}/${Output}.g.vcf \
        -log ${workdir}/${Output}.combine.log
  }
  output {
    File rawVCF = "${workdir}/${Output}.g.vcf"
  }
}
