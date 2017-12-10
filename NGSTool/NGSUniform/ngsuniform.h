#ifndef ngsqc_H
#define ngsqc_H
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <zlib.h>
#include <sys/stat.h>
typedef struct {
	char *read1;
	char *read2;
	char *fastq1;
	char *fastq2;
	int lenth;
}_Opt_;
typedef struct{
	char *seq;
	char *id;
	char *qual;
}_Fastq_;


// extern function
// processing line




#endif	// end BIOIO

