#include <ngsduplicate.h>
int usage(){
	printf("#--------------------------------------------- usage ------------------------------------------------------\n");
	printf("Example: ./ngsduplicate -1 read1.fq -2 read2.fq -q 33 -k demo -s a.stat\n");
	printf("Option:\topt\ttype\tdescription\n");
	printf("\t-1\tfile\tfastq read1 file,gz or not gz\n");
	printf("\t-2\tfile\tfastq read2 file,gz or not gz\n");
	printf("\t-q\tnum\tquality base,default 33\n");
	printf("\t-s\tfile\toutput duplication stat file [forced]\n");
	printf("\t-a\tfile\toutput uniq fragment 1.fq.gz\n");
	printf("\t-b\tfile\toutput uniq fragment 2.fq.gz\n");
	printf("\t-c\tfile\toutput duplication fragment 1.fq.gz\n");
	printf("\t-d\tfile\toutput duplication fragment 2.fq.gz\n");
	printf("\t-h\tusage\n");
	printf("#--------------------------------------------- usage ------------------------------------------------------\n");
	// over
	exit(0);
}
_Opt_ *Get_opt(int argc, char *argv[]){
	_Opt_ *opt;
	opt=malloc(sizeof(_Opt_));
	int opti;
	opt->qual=33;
	if( argc <= 1 ) {usage();};
	while((opti = getopt(argc, argv, "1:2:q:a:b:s:c:d:h:")) != -1){
		switch (opti){
		case '1': opt->read1 = malloc(strlen(optarg)+1);strcpy(opt->read1,optarg);break;
		case '2': opt->read2 = malloc(strlen(optarg)+1);strcpy(opt->read2,optarg);break;
		case 's': opt->ostat = malloc(strlen(optarg)+1);strcpy(opt->ostat,optarg);break;
		case 'a': opt->ofq1 = malloc(strlen(optarg)+1);strcpy(opt->ofq1,optarg);break;
		case 'b': opt->ofq2 = malloc(strlen(optarg)+1);strcpy(opt->ofq2,optarg);break;
		case 'c': opt->ofq3 = malloc(strlen(optarg)+1);strcpy(opt->ofq3,optarg);break;
		case 'd': opt->ofq4 = malloc(strlen(optarg)+1);strcpy(opt->ofq4,optarg);break;
		case 'q': opt->qual = atoi(optarg);break;
		case 'h': usage(); break;
		};
	};
	if (opt->read1 == NULL){usage();};
	if (opt->read2 == NULL){usage();};
	if (opt->ostat == NULL){usage();};
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
	gzFile FP1;
	gzFile FP2;
	if( (FP1=gzopen(opt->read1, "rb")) == NULL ){printf("open file %s is failed!", opt->read1);exit(0);};
    if( (FP2=gzopen(opt->read2, "rb")) == NULL ){printf("open file %s is failed!", opt->read2);exit(0);};
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
	char *readmerge;
	char *qualmerge;
	readmerge=malloc(sizeof(char)*500);
	qualmerge=malloc(sizeof(char)*500);
	_Dupstat_ *Dupstat=NULL;
	_Dupstat_ *tmp;
	int sum_reads=0;
	while ((Get_read(FP1,read1))!=NULL && (Get_read(FP2,read2))!=NULL)
	{
		sum_reads++;
		strcat(readmerge,read1->seq);
		strcat(readmerge," ");
		strcat(readmerge,read2->seq);
		HASH_FIND_STR(Dupstat,readmerge,tmp);
		strcat(qualmerge,read1->qual);
		strcat(qualmerge," ");
		strcat(qualmerge,read2->qual);
		if (tmp == NULL)
		{
			tmp=malloc(sizeof(_Dupstat_));
			tmp->seq=malloc(strlen(readmerge)+1);
			strcpy(tmp->seq,readmerge);
			tmp->qual=malloc(strlen(qualmerge)+1);
			strcpy(tmp->qual,qualmerge);
			tmp->num=1;
			HASH_ADD_STR(Dupstat,seq,tmp);
		}else{
			tmp->num++;
		}
		readmerge[0]='\0';
		qualmerge[0]='\0';
	//printf("%d\n",sum_reads);
	}
	
	gzclose(FP1);
	gzclose(FP2);
	FILE *FS;
	if( (FS=fopen(opt->ostat, "w")) == NULL ){printf("open file %s is failed!", opt->ostat);exit(0);};
	gzFile FQ1,FQ2,FQ3,FQ4;
	if (opt->ofq1 != NULL && opt->ofq2 !=NULL )
	{
		if( (FQ1=gzopen(opt->ofq1, "w")) == NULL ){printf("open file %s is failed!", opt->ofq1);exit(0);};
		if( (FQ2=gzopen(opt->ofq2, "w")) == NULL ){printf("open file %s is failed!", opt->ofq2);exit(0);};
	}
	if (opt->ofq3 != NULL && opt->ofq4 !=NULL )
	{
		if( (FQ3=gzopen(opt->ofq3, "w")) == NULL ){printf("open file %s is failed!", opt->ofq3);exit(0);};
		if( (FQ4=gzopen(opt->ofq4, "w")) == NULL ){printf("open file %s is failed!", opt->ofq4);exit(0);};
	}

	int dup_frag=0;	
	int total_frag=HASH_COUNT(Dupstat);
	char *seq1,*qual1;
	char *seq2,*qual2;
	int seqid=0;
	int dupseqid=0;
	int dup_read=0;
	int uniq_read=0;
	int uniq_frag=0;
	for(tmp=Dupstat; tmp != NULL; tmp=tmp->hh.next) {  
		if (tmp->num > 1)
		{
			dup_frag++;
			dup_read+=tmp->num;
			if (opt->ofq3 != NULL && opt->ofq4 != NULL)
			{
				dupseqid++;
				seq1 = strtok (tmp->seq," ");
				seq2 = strtok (tmp->seq," ");
				qual1 = strtok (tmp->qual," ");
				qual2 = strtok (tmp->qual," ");
				gzprintf(FQ3,"@duplicate:%d:fragment:%d\n",tmp->num,dupseqid);
				gzprintf(FQ3,"%s\n",seq1);
				gzprintf(FQ3,"+\n");
				gzprintf(FQ3,"%s\n",qual1);
				gzprintf(FQ4,"@duplicate:%d:fragment:%d\n",tmp->num,dupseqid);
				gzprintf(FQ4,"%s\n",seq2);
				gzprintf(FQ4,"+\n");
				gzprintf(FQ4,"%s\n",qual2);
			}
		}else{
			seqid++;
			uniq_frag++;
			uniq_read++;
			if (opt->ofq1 != NULL && opt->ofq2 != NULL)
			{
				seq1 = strtok (tmp->seq," ");
				seq2 = strtok (tmp->seq," ");
				qual1 = strtok (tmp->qual," ");
				qual2 = strtok (tmp->qual," ");
				gzprintf(FQ1,"@uniq:fragment:%d\n",seqid);
				gzprintf(FQ1,"%s\n",seq1);
				gzprintf(FQ1,"+\n");
				gzprintf(FQ1,"%s\n",qual1);
				gzprintf(FQ2,"@uniq:fragment:%d\n",seqid);
				gzprintf(FQ2,"%s\n",seq2);
				gzprintf(FQ2,"+\n");
				gzprintf(FQ2,"%s\n",qual2);
			}
		}
	}  
	fprintf(FS,"#Type\tTotal Reads\tDump Reads\tUniq Read\tTotal Fragment\tDump Fragment\tUniq Fragment\n");
	fprintf(FS,"Number\t%d\t%d\t%d\t%d\t%d\t%d\n",sum_reads,dup_read,uniq_read,total_frag,dup_frag,uniq_frag);
	fprintf(FS,"Percentage(%)\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",1.0*sum_reads/sum_reads*100,1.0*dup_read/sum_reads*100,1.0*uniq_read/sum_reads*100,1.0*total_frag/total_frag*100,1.0*dup_frag/total_frag*100,1.0*uniq_frag/total_frag*100);
	fclose(FS);
	if (opt->ofq1 != NULL && opt->ofq2 != NULL)
	{
		gzclose(FQ1);
		gzclose(FQ2);
	}
	if (opt->ofq3 != NULL && opt->ofq4 != NULL){
		gzclose(FQ3);
		gzclose(FQ4);
	}
	time=(clock()-time);
	printf("Our program is Done! It takes %d s\n",time/1000000);
	return(0);
	
}
