workflow Gvcftyping {
  call gvcftyping
}

task gvcftyping {
  File inputVCFs
  String workdir
  File RefFasta
  File Refindex
  File Refdict
  command {
    java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar \
        -T  GenotypeGVCFs\
        -R ${RefFasta} \
        -V ${inputVCFs} \
        -o ${workdir}/pop.noid.vcf \
	-nt 16 \
	-log ${workdir}/pop.vcf-typing.log
  }
  output {
    File rawVCF = "${workdir}/pop.noid.vcf"
  }
}

