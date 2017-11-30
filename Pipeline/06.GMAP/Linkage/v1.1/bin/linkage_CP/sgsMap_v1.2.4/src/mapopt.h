#ifndef MAPOPT_H
#define MAPOPT_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include "bioio.h"
#include "biomemory.h"
#include "biostring.h"
#include "bioerror.h"


/* Map option instruction : employed from Joinmap 
 * 
 * === Map building ===
 *   -- Spatial sampling threshold (five): default {[0.1], [0.05], [0.03], [0.02], [0.01]}
 *   -- Nr. of map optimization rounds per sample: default [3]
 *
 * === Map order optimization === 
 *   -- Chain length(with constant acc.prob): default [1000]
 *   -- Initial acceptance probability: default [0.25] 
 *   -- Cooling control parameter: default [0.001] 
 *   -- Stop after # chains without improvement: default [1000] 
 *
 * === Multipoint estimation of recombination frequencies
 *   -- Length of burn-in chain: default [1000]
 *   -- Nr. of Monte Carlo EM cycles: default [4]
 *   -- Chain length per Monte Carlo EM cycle: default [1000]
 *   -- sampling period for rec.freq matrix samples: default [5] 
 *
 * === Plausible map positions
 *   -- Acceptance control parameter: default [1]
 *   -- Burn-in chain length: default [10000]
 *   -- Number of map samples to draw: default [1000]
 *   -- Sampling period for map samples: default [1000]
 *
 */

// define MapOpt
typedef struct {
	// map building: spatial sampling 
	double *SS_rec_thresholds;
	int MB_nopt;
	
	// map optimization: simulating annealing
	int SA_nstep;
	double SA_initial_acceptance;
	double SA_cooling_control;
	int SA_unimproved_chain_max;
	
	// multipoint estimation of recombination frequencies: Gibbs sampling 
	int Gibbs_burn_in;
	int Gibbs_ncycle_MCEM;
	int Gibbs_chain_length;
	int Gibbs_sample_period;
	int Gibbs_global_sampling;

	// generate plausible position matrix: Metropolis Hasting algorithm
	double PP_acceptance_control;
	int PP_burn_in;
	int PP_nsample;
	int PP_sample_period;
	
	// map function 
	int mapFunc; // 0 for Handane, 1 for Kosambi 

	// map IO option 
	char *locFile;			// locus genotype file
	char *pwdFile;			// pwd file 
	char *startOrderFile;   // start order file  
	char *fixOrderFile;		// fix order file 
	char *fKey;		        // key of output files 
	char *sexAverMapFile;   // integrited map file 
	char *maleMapFile;      // male map file 
	char *femaleMapFile;    // female map file 
	char *sexAverPwdFile;   // sexAver pwd file: final multi-point estimation 
	char *malePwdFile;      // male pwd file 
	char *femalePwdFile;    // female pwd file 
	char *malePlausiblePosFile;
	char *femalePlausiblePosFile;
	char *singletonFile;
	char *haplotypeFile;
	char *logFile;			// log file of mapping process 
	FILE *logfp;			// file pointer of log file

}MapOpt;

// extern methods
extern MapOpt *get_map_opt(int argc, char *argv[]);
extern void write_map_opt_into_logfile(MapOpt *opt);

#endif


