/*
 *
 * main of mapping using maximum likelihood algorithm 
 *
 */

// Copyright 2012, Ma Chouxian <macx@biomarker.com.cn>

#include "sgsMap.h"

#define MAPTERM_TEST 
#ifdef MAPTERM_TEST

int main(int argc, char *argv[]) {
	// stat the time of map sort
	double cpubegin = clock();
	long runbegin = time(NULL);
	
	// parse command line option 
	MapOpt *opt = parse_command_line(argc, argv);
	
	// read loc and pwd info 
	LocsINFO *locsINFO = read_loc_pwd_info(opt);

	// get linear arrangement of loci in a linkage group 
	map_building(locsINFO, opt);

	// over and calculate elapsed time 
	double cputime = (clock() - cpubegin) / CLOCKS_PER_SEC;
	long runtime = time(NULL) - runbegin;
	fprintf(opt->logfp, "CPU time using: %.0lf seconds or %.3lf mins or %.3lf hours\n", cputime, cputime/60, cputime/3600);
	fprintf(opt->logfp, "Run time using: %ld seconds or %.3lf mins or %.3lf hours\n", runtime, runtime/60.0, runtime/3600.0);
	fclose(opt->logfp); // close log file 

	printf("Done.Total elapsed time :%ld seconds or %.2lf mins or %.2lf hours\n",runtime, runtime/60.0, runtime/3600.0);
	return 0;
}

#endif
