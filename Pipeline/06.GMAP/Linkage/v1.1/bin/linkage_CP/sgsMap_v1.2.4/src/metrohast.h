#ifndef METROHAST_H
#define METROHAST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "map.h"
#include "mapopt.h"
#include "siman.h"

//===================================
// metropolis hasting struct 
//===================================

typedef struct{
	GenoMatrix **gMatrix;  // genotype matrix 
	RecLOD **recLOD;       // rec.freq matrix 
	
	// indi 
	double nind;

	// map
	int *map;             // current sex-averaged map 
	int maplen;              // s map len
	double init_E;          // initial E: best map 

	int **pmMapPos;     // initial p and m map position of loci 
	int **p_plausible_position;
	int plen;
	int **m_plausible_position;
	int mlen;

	/* gsl random number generator */
	const gsl_rng_type *T;
	gsl_rng *r;

	// calculateion parameters 
	double control_probability;  // acceptance control parameters
	double burn_in;            // length of burn-in chain 
	double n_samples;          // number of map samples to draw 
	double sample_period;      // sampling period for map samples
	
}MetroHast;

// extern functions 
extern inline void generate_plausible_position_matrix(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);

#endif
