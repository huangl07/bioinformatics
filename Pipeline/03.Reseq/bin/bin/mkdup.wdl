workflow MarkDuplicate {
	String sample
	File bam
	String workdir
	call markduplicate{
	input:
		sample=sample,
		bam=bam,
		workdir=workdir
	}
}
task markduplicate {
	String sample
	File bam
	String workdir
	command {
		java -jar /mnt/ilustre/users/dna/.env/bin/picard.jar MarkDuplicates \
			TMP_DIR=${workdir}/MKDUP/ \
			MAX_FILE_HANDLES=100 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true \
			I=${bam} \
			O=${workdir}/${sample}.mkdup.bam \
			M=${workdir}/${sample}.metric \
			CREATE_INDEX=TRUE
	}
	output {
		File outbam = "${workdir}/${sample}.mkdup.bam"
		File metric = "${workdir}/${sample}.metric"
	}
}
