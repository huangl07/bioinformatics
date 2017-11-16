#include <ngsuniform.h>
#include <limits.h>

int usage(){
	printf("#--------------------------------------------- usage ------------------------------------------------------\n");
	printf("Example: ./ngsqc -1 read1.fq -2 read2.fq -q 33 -k demo -o ./\n");
	printf("Option:\topt\ttype\tdescription\n");
	printf("\t-1\tfile\tfastq read1 file,gz or not gz\n");
	printf("\t-2\tfile\tfastq read2 file,gz or not gz\n");
	printf("\t-l\tnum\tlength to uniform\n");
	printf("\t-a\tfile\toutput read1 file\n");
	printf("\t-b\tfile\toutput read2 file\n");
//	printf("\t-k\tstr\toutput keys of filename,always samples ID\n");
	printf("\t-h\tusage\n");
	printf("#--------------------------------------------- usage ------------------------------------------------------\n");
	// over
	exit(0);
}
_Opt_ *Get_opt(int argc, char *argv[]){
	_Opt_ *opt;
	opt=malloc(sizeof(_Opt_));
	opt->lenth=-1;
	int opti;
	if( argc <= 1 ) {usage();};
	while((opti = getopt(argc, argv, "1:2:a:b:l:")) != -1){
		switch (opti){
		case '1': opt->read1 = malloc(strlen(optarg)+1);strcpy(opt->read1,optarg);break;
		case '2': opt->read2 = malloc(strlen(optarg)+1);strcpy(opt->read2,optarg);break;
		case 'a': opt->fastq1 = malloc(strlen(optarg)+1);strcpy(opt->fastq1,optarg);break;
		case 'b': opt->fastq2 = malloc(strlen(optarg)+1);strcpy(opt->fastq2,optarg);break;
		case 'l': opt->lenth = atoi(optarg);break;
		case 'h': usage(); break;
		};
	};
	if (opt->read1 == NULL){usage();};
	if (opt->read2 == NULL){usage();};
	if (opt->fastq1 == NULL){usage();};
	if (opt->fastq2 == NULL){usage();};
	if (opt->lenth == -1){usage();};
	return opt;
}
void *Get_read(gzFile FP,_Fastq_ *read){
	if (gzgets(FP,read->id, 256)!=NULL)
	{
		gzgets(FP,read->seq,256);
		gzgets(FP,read->qual,256);
		gzgets(FP,read->qual,256);
	}else{
		return NULL;
	}
	if(read->seq[strlen(read->seq)-1] == '\n'){read->seq[strlen(read->seq)-1] = '\0';}
	if(read->id[strlen(read->id)-1] == '\n'){read->id[strlen(read->id)-1] = '\0';}
	if(read->qual[strlen(read->qual)-1] == '\n'){read->qual[strlen(read->qual)-1] = '\0';}
	return read;
}

int main(int argc, char *argv[]){
	int time=clock();
	_Opt_ *opt=Get_opt(argc,argv);
	gzFile FP1,FP2;
	if( (FP1=gzopen(opt->read1, "rb")) == NULL ){printf("open file %s is failed!", opt->read1);exit(0);};
    if( (FP2=gzopen(opt->read2, "rb")) == NULL ){printf("open file %s is failed!", opt->read2);exit(0);};
	gzFile OFP1,OFP2;
	if( (OFP1=gzopen(opt->fastq1, "w")) == NULL ){printf("open file %s is failed!", opt->fastq2);exit(0);};
	if( (OFP2=gzopen(opt->fastq2, "w")) == NULL ){printf("open file %s is failed!", opt->fastq1);exit(0);};
	_Fastq_	*read1;
	_Fastq_ *read2;
	read1=malloc(sizeof(_Fastq_));
	read1->id=malloc(sizeof(char)*256);
	read1->qual=malloc(sizeof(char)*256);
	read1->seq=malloc(sizeof(char)*256);
	read2=malloc(sizeof(_Fastq_));
	read2->id=malloc(sizeof(char)*256);
	read2->qual=malloc(sizeof(char)*256);
	read2->seq=malloc(sizeof(char)*256);
	while ((Get_read(FP1,read1))!=NULL && (Get_read(FP2,read2))!=NULL)
	{
		if(strlen(read1->seq) < opt->lenth){continue;}
		if(strlen(read2->seq) < opt->lenth){continue;}
		read1->seq[opt->lenth]='\0';
		read1->qual[opt->lenth]='\0';
		read2->seq[opt->lenth]='\0';
		read2->qual[opt->lenth]='\0';
		gzprintf(OFP1,"%s\n",read1->id);
		gzprintf(OFP1,"%s\n",read1->seq);
		gzprintf(OFP1,"+\n");
		gzprintf(OFP1,"%s\n",read1->qual);
		gzprintf(OFP2,"%s\n",read2->id);
		gzprintf(OFP2,"%s\n",read2->seq);
		gzprintf(OFP2,"+\n");
		gzprintf(OFP2,"%s\n",read2->qual);
	}
	gzclose(FP1);
	gzclose(FP2);
	gzclose(OFP1);
	gzclose(OFP2);
	time=(clock()-time);
	printf("Our program is Done! It takes %d s\n",time/1000000);
	return(0);
	
}
