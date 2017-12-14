workflow Gvcftyping {
  call gvcftyping
}

task gvcftyping {
  File inputVCFs
  String workdir
  File RefFasta
  File Refindex
  File Refdict
  File Internal
  String Filename
  command {
    java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar \
        -T  GenotypeGVCFs\
        -R ${RefFasta} \
        -V ${inputVCFs} \
        -o ${workdir}/${Filename}.noid.vcf \
        --never_trim_vcf_format_field \
		-L ${Internal}
	-jdk_inflater \
	-jdk_deflater \
	-nt 32 \
	-log ${workdir}/${Filename}.vcf-typing.log
  }
  output {
    File rawVCF = "${workdir}/${Filename}.noid.vcf"
  }
}

