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
		type="INDEL",
		workdir=workdir
	}
	call filterINDEL {
	input:
		RefFasta=refFasta, 
		RefIndex=refIndex, 
		RefDict=refDict, 
		type="INDEL",
		rawVCF=select.VCF,
		workdir=workdir
	}
	call recodeINDEL {
	input:
		RefFasta=refFasta, 
		RefIndex=refIndex, 
		RefDict=refDict, 
		type="INDEL",
		rawVCF=filterINDEL.VCF,
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
			-selectType ${type} \
			-nt 8 \
			-o ${workdir}/pop.${type}.vcf
	}
	output {
		File VCF = "${workdir}/pop.${type}.vcf"
	}
}
task filterINDEL {
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
		--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoef < -0.5 ||SQR > 10.0" \
		--filterName Failer
	}
	output {
		File VCF = "${workdir}/pop.${type}.filtered.vcf"
	}
}
task recodeINDEL {
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
		--setFilteredGtToNocall \
		--excludeFiltered \
		-nt 8 \
		--excludeNonVariants	
	}
	output {
		File VCF = "${workdir}/pop.${type}.final.vcf"
	}
}
