#ifndef SIMAN_H
#define SIMAN_H

#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_machine.h>
#include <string.h>
#include <math.h>
#include "bioerror.h"
#include "biomemory.h"
#include "map.h"
#include "mapopt.h"

typedef struct {
	int *mapOrder;   // map order 
	int maplen;      // map length 	
	int *rLocs;     // random loci space, {non_fix, fix_have_walking_space}, random loci picked in sa process, because of the restriction of fixorder, this list may be empty 
	int rLocsLen;    // length of rLocs, if rLocsLen equals to 0, means all loci fixed 
	
	double p_sarf; // optimization criterion of paternal map
	double m_sarf; // optimization criterion of maternal map
	double E;      // optimization criterion : p_sarf + m_sarf
}Map;

typedef struct {
	Map *curMap; // current map 
	Map *newMap; // new map generated from curMap 
	Map *bestmap;  // the best map found  so far 
	
	LocsINFO *locsINFO; // locus information 
	/* the fix order info : if loci in fix order, then isfixed[loci] == 1;  */
	int *isfixed;   // a link to last fix order --> the mapping process ensure that the last fixed order contains ohters    

	/* stat simulating annealing process */
	int n_evals;    // Nr. of evaluate Ef
	int n_iters;    // # iterations performed before system stops
	int n_rejects;  // Nr. of refusing new orders  
	int n_accepts;  // Nr. of accepting new orders 
	int n_eless;    // Nr. of downhill steps 
	double t_final; // final temperature
	
	/* gsl random number generator */
	const gsl_rng_type *T;
	gsl_rng *r;

	// simulating annealing parameters 
	int nstep;      // Chain length(with constant acc.prob)
	double t_initial; // Initial acceptance probability, or initial temperature
	double cooling_factor;  // cooling factor : cooling_factor = 1/(1 + cooling control parameter)
	int unimproved_chain_max; // Stop after # chains without improvement
}simAN;

/* extern functions */
extern inline int map_optimization (LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);
extern inline double boltzmann(double E, double new_E, double T);
extern inline void repalce_random_locus(int *newMap, int rLocus, int rPos);

#endif
