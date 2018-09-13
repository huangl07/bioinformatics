workflow SelectVariant {
	File refFasta
	File refIndex
	File refDict
	File rawvcf
	String workdir
	call select {
	input:
		RefFasta=refFasta,
		RefIndex=refIndex,
		RefDict=refDict,
		rawVCF=rawvcf,
		type="SNP",
		workdir=workdir
	}
	call filterSNP {
	input:
		RefFasta=refFasta, 
		RefIndex=refIndex, 
		RefDict=refDict, 
		type="SNP",
		rawVCF=select.VCF,
		workdir=workdir
	}
	call recodeSNP {
	input:
		RefFasta=refFasta, 
		RefIndex=refIndex, 
		RefDict=refDict, 
		type="SNP",
		rawVCF=filterSNP.VCF,
		workdir=workdir
	}

}
task select {
	File RefFasta
	File RefIndex
	File RefDict
	String type
	File rawVCF
	String workdir
	command {
		java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar  \
			-T SelectVariants \
			-R ${RefFasta} \
			-V ${rawVCF} \
			-nt 8 \
			-selectType ${type} \
			-o ${workdir}/pop.${type}.vcf
	}
	output {
		File VCF = "${workdir}/pop.${type}.vcf"
	}
}
task filterSNP {
	File RefFasta
	File RefIndex
	File RefDict
	File rawVCF
	String workdir
	String type
	command {
		java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar  \
		-T VariantFiltration \
		-R ${RefFasta} \
		-V ${rawVCF} \
		-o ${workdir}/pop.${type}.filtered.vcf \
		--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SQR > 3.0" \
		--filterName Failer
	}
	output {
		File VCF = "${workdir}/pop.${type}.filtered.vcf"
	}
}
task recodeSNP {
	File RefFasta
	File RefIndex
	File RefDict
	File rawVCF
	String workdir
	String type
	command {
		java -jar /mnt/ilustre/users/dna/.env//bin//GenomeAnalysisTK.jar  \
		-T SelectVariants \
		-R ${RefFasta} \
		-V ${rawVCF} \
		-o ${workdir}/pop.${type}.final.vcf \
		-nt 8 \
		--setFilteredGtToNocall \
		--excludeFiltered \
		--excludeNonVariants	
	}
	output {
		File VCF = "${workdir}/pop.${type}.final.vcf"
	}
}
