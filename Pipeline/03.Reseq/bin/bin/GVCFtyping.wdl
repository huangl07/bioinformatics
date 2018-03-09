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
  String NT
  String Xmx
  command {
    java -XX:+UseSerialGC -Xmx${Xmx} -Djava.io.tmpdir=${workdir}/tmp/ -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar \
        -T  GenotypeGVCFs\
        -R ${RefFasta} \
        -V ${inputVCFs} \
        -o ${workdir}/${Filename}.noid.vcf \
        --never_trim_vcf_format_field \
	--max_alternate_alleles 3 \
	-L ${Internal} \
	-jdk_inflater \
	-jdk_deflater \
	-nt ${NT} \
	-log ${workdir}/${Filename}.vcf-typing.log \
	 --disable_auto_index_creation_and_locking_when_reading_rods 
  }
  output {
    File rawVCF = "${workdir}/${Filename}.noid.vcf"
  }
}

