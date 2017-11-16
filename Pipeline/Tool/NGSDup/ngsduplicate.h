#ifndef ngsqc_H
#define ngsqc_H
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "/mnt/ilustre/users/long.huang/bin/C/lib/zlib/zlib.h"  
#include <sys/stat.h>
#include "/mnt/ilustre/users/long.huang/Software/uthash-master/include/uthash.h"
typedef struct {
	char *read1;
	char *read2;
	char *ofq1;
	char *ofq2;
	char *ofq3;
	char *ofq4;
	int qual;
	char *ostat;
}_Opt_;
typedef struct{
	char *seq;
	char *id;
	char *qual;
}_Fastq_;

typedef struct{  
	char  *seq; /* key */  
	int num;/*duplication numbers*/
	char *qual;
	UT_hash_handle hh; /* makes this structure hashable */  
}_Dupstat_;  
// extern function
// processing line




#endif	// end BIOIO

