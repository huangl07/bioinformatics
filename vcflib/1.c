#include <htslib/vcf.h>
#include <config.h>
#include <stdio.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
void main(int argc,char **argv){
	char *fname = argc>1 ? argv[1] : "rmme.bcf";
	htsFile *fp    = hts_open(fname,"rb");
//	int nsample=bcf_hdr_nsample(hdr)
	bcf_hdr_t *hdr=bcf_hdr_read(fp);
	bcf_hdr_t=bcf_hdr_dup(hdr);
	while(bcf_read1(fp,hdr,rec)>0){
	`		
	}
}


