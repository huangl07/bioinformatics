#ifndef ngsqc_H
#define ngsqc_H
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "/mnt/ilustre/users/long.huang/bin/C/lib/zlib/zlib.h"  
#include <sys/stat.h>
typedef struct {
	char *read1;
	char *read2;
	char *dir;
	char *key;
	char *seq1;
	char *seq2;
	int qual;
}_Opt_;
typedef struct{
	char *seq;
	char *id;
	char *qual;
}_Fastq_;


// extern function
// processing line




#endif	// end BIOIO

