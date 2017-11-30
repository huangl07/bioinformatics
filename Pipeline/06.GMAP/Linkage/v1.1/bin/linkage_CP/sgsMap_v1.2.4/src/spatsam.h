#ifndef SPATSAM_H
#define SPATSAM_H


#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "bioerror.h"
#include "biomemory.h"
#include "map.h"
#include "mapopt.h"

/* define spatial sampling struct */
typedef struct {
	int round;             // sample stage 
	double radius;         // sampling radius 
	int *samSpace;         // locus set of spatial sampling, always tspm->T
	int ssl;               // space len
	int *isAnchor;         // indictive array of anchor loci 
	
	GenoMatrix **gMatrix;  // genotype matrix 
	RecLOD **recLOD;       // pairwise rec.matrix info 

	// gsl random number generator 
	const gsl_rng_type *T;
	gsl_rng *r;

	// sampling result 
	int *sample;              // a spatial sample 
	int samLen;               // sample size 
	int *unSampled;           // unsampled 
	int unSamLen;             // unsampled size 
}SpatSAM;

// extern function 
extern int spatial_sampling(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);

#endif
