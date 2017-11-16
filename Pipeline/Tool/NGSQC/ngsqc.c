#include <ngsqc.h>
#include <limits.h>

int usage(){
	printf("#--------------------------------------------- usage ------------------------------------------------------\n");
	printf("Example: ./ngsqc -1 read1.fq -2 read2.fq -q 33 -k demo -o ./\n");
	printf("Option:\topt\ttype\tdescription\n");
	printf("\t-1\tfile\tfastq read1 file,gz or not gz\n");
	printf("\t-2\tfile\tfastq read2 file,gz or not gz\n");
	printf("\t-q\tnum\tquality base,default 33\n");
	printf("\t-o\tdir\toutput dir\n");
	printf("\t-k\tstr\toutput keys of filename,always samples ID\n");
	printf("\t-e\tstr\tenzyme 1\n");
	printf("\t-p\tstr\tenzyme 2\n");
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
	while((opti = getopt(argc, argv, "1:2:o:k:q:e:p:")) != -1){
		switch (opti){
		case '1': opt->read1 = malloc(strlen(optarg)+1);strcpy(opt->read1,optarg);break;
		case '2': opt->read2 = malloc(strlen(optarg)+1);strcpy(opt->read2,optarg);break;
		case 'o': opt->dir = malloc(strlen(optarg)+1);strcpy(opt->dir,optarg);break;
		case 'k': opt->key = malloc(strlen(optarg)+1);strcpy(opt->key,optarg);break;
		case 'q': opt->qual = atoi(optarg);break;
		case 'h': usage(); break;
		};
	};
	if (opt->read1 == NULL){usage();};
	if (opt->read2 == NULL){usage();};
	if (opt->key == NULL){usage();};
	if (opt->seq1 == NULL){opt->seq1="HHHH";};
	if (opt->seq2 == NULL){opt->seq2="HHHH";};
	if (opt->dir == NULL){opt->dir=malloc(sizeof("./"));strcpy(opt->dir,",.");};
	return opt;
}
int main(int argc, char *argv[]){
	int time=clock();
	_Opt_ *opt=Get_opt(argc,argv);
	gzFile FP1,FP2;
//	printf("%c\t%c\t%c\t\%c\n",opt->read1[strlen(opt->read1)],opt->read1[strlen(opt->read1)-1],opt->read1[strlen(opt->read1)-2],opt->read1[strlen(opt->read1)-3]);
//	opt->read1[strlen(opt->read1)]='\0';
//	printf("%c\t%c\t%c\t\%c\n",opt->read1[strlen(opt->read1)],opt->read1[strlen(opt->read1)-1],opt->read1[strlen(opt->read1)-2],opt->read1[strlen(opt->read1)-3]);
//	opt->read2[strlen(opt->read2)]='\0';
	if( (FP1=gzopen(opt->read1, "rb")) == NULL ){printf("open file %s is failed!", opt->read1);exit(0);};
    if( (FP2=gzopen(opt->read2, "rb")) == NULL ){printf("open file %s is failed!", opt->read2);exit(0);};
	if (access(opt->dir,0)==-1){ if((mkdir(opt->dir,0755))==-1){printf("error make dir!pleasecheack!");exit(0);} }
	char outstat[1000];sprintf(outstat, "%s/%s.stat",opt->dir,opt->key);
	char outqual[1000];sprintf(outqual, "%s/%s.qual",opt->dir,opt->key);
	char outatgc[1000];sprintf(outatgc, "%s/%s.atgc",opt->dir,opt->key);
	char line1[256];
	char line2[256];
	char qual1[256];
	char qual2[256];
	char seq1[256];
	char seq2[256];
	char *id1;
	char *id2;
	long long  ATGCN1[100]={0};
	long long  QUAL1[100]={0};
	long long  Pos_ATGCN1[500][100]={0};
	long long  Pos_QUAL1[500][100]={0};
	long long  ATGCN2[100]={0};
	long long  QUAL2[100]={0};
	long long Pos_ATGCN2[500][100]={0};
	long long Pos_QUAL2[500][100]={0};
	long long readnum=0;
	long long enzy1=0;
	long long enzy2=0;
	long long enzy1_read1=0;
	long long enzy1_read2=0;
	long long enzy2_read1=0;
	long long enzy2_read2=0;
	int lenth;
	int maxqua=0;
	int maxlen=0;
	int i,j;
	while ((gzgets(FP1,line1, sizeof(line1))!=NULL) && (gzgets(FP2, line2,sizeof(line2)) !=NULL))
	{
		readnum++;
		id1 = strtok (line1," ");
		id2 = strtok (line2," ");
		if (strcmp(id1,id2) != 0){printf("not a pair of fastqfiles! please check!\n%s\n%s\n",opt->read1,opt->read2);exit(0);}
		gzgets(FP1,seq1, sizeof(seq1));lenth =strlen(seq1);
		gzgets(FP2,seq2, sizeof(seq2));lenth =strlen(seq2);
		gzgets(FP1,qual1, sizeof(qual1));
		gzgets(FP2,qual2, sizeof(qual1));
		gzgets(FP1,qual1, sizeof(qual1));lenth =strlen(qual1);
		gzgets(FP2,qual2, sizeof(qual2));lenth =strlen(qual2);
		if (strlen(qual1) != strlen(seq1)){printf("quality length diff read length! error! please check\n%s\n%s\n%s\n!",line1,seq1,qual1);exit(0);}
		if (strlen(qual2) != strlen(seq2)){printf("quality length diff read length! error! please check\n%s\n%s\n%s\n!",line2,seq2,qual2);exit(0);}
		if (strlen(seq1) > maxlen){maxlen=strlen(seq1);}
		if (strlen(seq2) > maxlen){maxlen=strlen(seq2);}
		//	printf("seq1:%s\n",seq1)
		for ( i =0;i<strlen(seq1);i++ )
		{
			ATGCN1[seq1[i]]++;
			Pos_ATGCN1[i][seq1[i]]++;
			int qua=qual1[i]-opt->qual;
			if (maxqua < qua){maxqua=qua;}
			QUAL1[qua]++;
			Pos_QUAL1[i][qua]++;
		}
		for ( i =0;i<strlen(seq2);i++ )
		{
			ATGCN2[seq2[i]]++;
			Pos_ATGCN2[i][seq2[i]]++;
			int qua=qual2[i]-opt->qual;
			if (maxqua < qua){maxqua=qua;}
			QUAL2[qua]++;
			Pos_QUAL2[i][qua]++;
		}
	}
	gzclose(FP1);
	gzclose(FP2);
	FILE *Out1,*Out2;
    if( (Out1=fopen(outstat, "w")) == NULL ){printf("open file %s is failed!", outstat);exit(0);};
	long long sum=ATGCN1['A']+ATGCN1['T']+ATGCN1['G']+ATGCN1['C']+ATGCN1['N']+ATGCN2['A']+ATGCN2['T']+ATGCN2['G']+ATGCN2['C']+ATGCN2['N'];
	long long  GC=ATGCN1['G']+ATGCN1['C']+ATGCN2['G']+ATGCN2['C'];
	long long Q30=0;
	long long Q20=0;
	long long Q=0;
	long long Q1=0;
	long long Q2=0;
	for (i=0;i<=maxqua;i++)
	{
		Q+=QUAL1[i]*i+QUAL2[i]*i;
		Q1+=QUAL1[i]*i;
		Q2+=QUAL2[i]*i;
	}
	for ( i=21;i<=maxqua;i++)
	{
		if (i >30)
		{
			Q30+=QUAL1[i]+QUAL2[i];
		}
		Q20+=QUAL1[i]+QUAL2[i];
	}
	long long sum1=ATGCN1['A']+ATGCN1['T']+ATGCN1['G']+ATGCN1['C']+ATGCN1['N'];
	long long GC1=ATGCN1['G']+ATGCN1['C'];
	long long Q30_1=0;
	long long Q20_1=0;
	for ( i=21;i<=maxqua;i++)
	{
		if (i >30)
		{
			Q30_1+=QUAL1[i];
		}
		Q20_1+=QUAL1[i];
	}
	long long sum2=ATGCN2['A']+ATGCN2['T']+ATGCN2['G']+ATGCN2['C']+ATGCN2['N'];
	long long GC2=ATGCN2['G']+ATGCN2['C'];
	long long Q30_2=0;
	long long Q20_2=0;
	for ( i=21;i<=maxqua;i++)
	{
		if (i >30)
		{
			Q30_2+=QUAL2[i];
		}
		Q20_2+=QUAL2[i];
	}

	fprintf(Out1,"#Flag\tstatics Flag\n");
	fprintf(Out1,"#Sample ID\tsampleID\n");
	fprintf(Out1,"#ReadNumber\treads number\n");
	fprintf(Out1,"#BaseNumber\tbase number\n");
	fprintf(Out1,"#GC(%)\tGC percentage(%)\n");
	fprintf(Out1,"#Q30(%)\tQ30 percentage(%)\n");
	fprintf(Out1,"#Q20(%)\tQ20 percentage(%)\n");
	fprintf(Out1,"#maxlen:%d\n",maxlen);
	fprintf(Out1,"#maxqua:%d\n",maxqua);
	fprintf(Out1,"#Flag\tSample ID\tRead Number\tBaseNumber\tA(%)\tT(%)\tG(%)\tC(%)\tN(%)\tGC(%)\tQ30(%)\tQ20(%)\tAverageQ\tenzyme1(%)\tenzyme2(%)\n");
	fprintf(Out1,"Total\t%s\t%lld\t%lld\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",opt->key,readnum,sum,1.0*(ATGCN1['A']+ATGCN2['A'])/sum*100,1.0*(ATGCN1['T']+ATGCN2['T'])/sum*100,1.0*(ATGCN1['G']+ATGCN2['G'])/sum*100,1.0*(ATGCN1['C']+ATGCN2['C'])/sum*100,1.0*(ATGCN1['N']+ATGCN2['N'])/sum*100,1.0*GC/sum*100,1.0*Q30/sum*100,1.0*Q20/sum*100,1.0*Q/sum);
	fprintf(Out1,"Read1\t%s\t%lld\t%lld\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",opt->key,readnum,sum1,1.0*ATGCN1['A']/sum1*100,1.0*ATGCN1['T']/sum1*100,1.0*ATGCN1['G']/sum1*100,1.0*ATGCN1['C']/sum1*100,1.0*ATGCN1['N']/sum1*100,1.0*GC1/sum1*100,1.0*Q30_1/sum1*100,1.0*Q20_1/sum1*100,1.0*Q1/sum1);
	fprintf(Out1,"Read2\t%s\t%lld\t%lld\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",opt->key,readnum,sum2,1.0*ATGCN1['A']/sum1*100,1.0*ATGCN1['T']/sum1*100,1.0*ATGCN1['G']/sum1*100,1.0*ATGCN1['C']/sum1*100,1.0*ATGCN1['N']/sum1*100,1.0*GC2/sum2*100,1.0*Q30_2/sum2*100,1.0*Q20_2/sum2*100,1.0*Q2/sum2);
	fclose(Out1);
    if( (Out1=fopen(outatgc, "w")) == NULL ){printf("open file %s is failed!", outatgc);exit(0);};
    if( (Out2=fopen(outqual, "w")) == NULL ){printf("open file %s is failed!", outqual);exit(0);};

	fprintf(Out1,"#pos\tA\tT\tG\tC\tN\n");
	fprintf(Out2,"#pos\t");
	for ( i=0;i<=maxqua;i++)
	{
		fprintf(Out2,"%d\t",i);
	}
	fprintf(Out2,"Aver\n");
	for ( i=0;i<maxlen;i++)
	{
		long long	sum=Pos_ATGCN1[i]['A']+Pos_ATGCN1[i]['T']+Pos_ATGCN1[i]['G']+Pos_ATGCN1[i]['C']+Pos_ATGCN1[i]['N'];
		long long Q=0;
		fprintf(Out2,"%d\t",i+1);
		for ( j=0;j<=maxqua;j++ )
		{
			Q+=Pos_QUAL1[i][j]*j;
			fprintf(Out2,"%.2f\t",1.0*Pos_QUAL1[i][j]/sum*100);
		}
		fprintf(Out2,"%.2f\n",1.0*Q/sum);
		fprintf(Out1,"%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",i+1,1.0*Pos_ATGCN1[i]['A']/sum*100,1.0*Pos_ATGCN1[i]['T']/sum*100,1.0*Pos_ATGCN1[i]['G']/sum*100,1.0*Pos_ATGCN1[i]['C']/sum*100,1.0*Pos_ATGCN1[i]['N']/sum*100);
	}
	for ( i=0;i<maxlen;i++)
	{
		long long	sum=Pos_ATGCN2[i]['A']+Pos_ATGCN2[i]['T']+Pos_ATGCN2[i]['G']+Pos_ATGCN2[i]['C']+Pos_ATGCN2[i]['N'];
		long long Q=0;
		fprintf(Out2,"%d\t",i+151);
		for ( j=0;j<=maxqua;j++ )
		{
			Q+=Pos_QUAL2[i][j]*j;
			fprintf(Out2,"%.2f\t",1.0*Pos_QUAL2[i][j]/sum*100);
		}
		fprintf(Out2,"%.2f\n",1.0*Q/sum);
		fprintf(Out1,"%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",i+151,1.0*Pos_ATGCN2[i]['A']/sum*100,1.0*Pos_ATGCN2[i]['T']/sum*100,1.0*Pos_ATGCN2[i]['G']/sum*100,1.0*Pos_ATGCN2[i]['C']/sum*100,1.0*Pos_ATGCN2[i]['N']/sum*100);
	}
	fclose(Out1);
	fclose(Out2);

	time=(clock()-time);
	printf("Our program is Done! It takes %d s\n",time/1000000);
	return(0);
	
}
