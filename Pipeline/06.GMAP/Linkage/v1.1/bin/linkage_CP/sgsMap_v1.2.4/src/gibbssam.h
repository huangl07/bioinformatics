#ifndef GIBBSSAM_H
#define GIBBSSAM_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "bioerror.h"
#include "biomemory.h"
#include "biostring.h"
#include "map.h"
#include "mapopt.h"

/* define Coordinate struct */
typedef struct {
	char type;       // type: 0 for abxcd & efxeg, 1 for hkxhk, 2 for lmxll, 3 for nnxnp
	char isMulti;    // is multiple haplotypes: 0 for miss observation, 1 for unique haplo, 2 for multiple haplos(hk)
	int indi;        // indi ID 
	int loci;        // loci index 
	int pPreLoc;     // prefix position in sMap of paternal map's loci 
	int pNextLoc;    // next position in sMap of paternal map's loci  
	int mPreLoc;     // prefix position in sMap of maternal map's loci 
	int mNextLoc;    // next position in sMap of maternal map's loci
	int sMapPos;     // position in sMap of loci
	int pMapPos;     // position in pMap 
	int mMapPos;     // position in mMap 

	int *pBinLoc;    // loci in a bin on paternal map and genotypes in indi are missing observations or partial informative
	int pBinLen;     // bin marker with loci on paternal map 
	int *mBinLoc;    // loci in a bin on maternal map and genotypes in indi are missing observations or partial informative 
	int mBinLen;     // bin marker with loci in maternal map 
	
	char pSamHap;    // sampled haplotype of P 
	char mSamHap;    // sampled haplotype of M 
}Coordinate;

typedef struct {
	double p_sarf_start;   // sum of paternal rec.freq's of adjacent segments at initial state 
	double m_sarf_start;   // sum of maternal rec.freq's of adjacent segments at initial state 
	
	double p_mrn_start;
	double m_mrn_start;
	
	double p_singleton_start;
	double m_singleton_start;
	
	double p_mrn_burn_in;
	double m_mrn_burn_in;
	
	double p_singleton_burn_in;
	double m_singleton_burn_in;
	
	double p_sarf_end;     // sum of paternal rec.freq's of adjacent segments at end of one MCEM cycle 
	double m_sarf_end;     // sum of maternal rec.freq's of adjacent segments at end of one MCEM cycle 

	double p_mrn_end;
	double m_mrn_end;

	double p_singleton_end;
	double m_singleton_end;
	
}GibbsStat;

typedef struct {
	
	// haplo matrix  
	HapElem ***hapMatrix;   // haplotype matrix 

	// recombination matrix 
	RecLOD **recMatrix;     // pairwise recombination info of sMap, pMap and mMap 
	RecLOD **recMatrixSum;  // sum of sampled recMatrix

	// sarf 
	double p_sarf;      // sum of adjacent paternal rec.freq 
	double m_sarf;      // sum of adjacent maternal rec.freq

	// mrn 
	double p_mrn;       
	double m_mrn;
	
	// singleton
	double pSingleton;  // sum of paternal singleton 
	double mSingleton;  // sum of maternal singleton 

}MCEM;

/* define GibbsSAM for Gibbs sampling */
typedef struct {
	// genotype and haplotype matrix
	int nind;
	GenoMatrix **gMatrix;     // sampling gMatrix 
	HapElem ***hapMatrix;       // initial haplo Matrix 
	
	// pairwise rec matrix 
	RecLOD **recMatrix;       // recombination frequency matrix: a link to locsINFO->recLOD

	// sample space 
	Coordinate **sampleSpace;  // sample from this list 
	int sampleLen;            // length of sample list  
	
	// mcem cycles 
	MCEM *mcem;               // MCEM cycles 

	// statistic 
	GibbsStat **mcemMornitor;        // statistic for mornitoring gibbs sampling 
	
	// gsl random number generator 
	const gsl_rng_type * T;
	gsl_rng *r;
	
	// map info 
	int **pmMapPos;           // loci position info in sex-specific maps, 1D smap order, 2d - p and m pos  
	int *sMap;			      // map order used in multi-point estimation of recombination frequency
	int sLen;                 // sMap len 
	double *SD;               // sMap distance 
	double s_singleton;       // singleton of sMap
	double s_sarf;            // sum of adjacent sexaver rec.freq 
	int *pMap;                // paternal map 
	int pLen;                 // paternal map len 
	double *PD;                // paternal map distance 
	double p_singleton;       // singleton of pMap
	double p_sarf;
	int *mMap;                // maternal map 
	int mLen;                 // maternal map len
	double *MD;                // maternal map distance 
	double m_singleton;       // singleton of mMap 
	double m_sarf;
	
	// calculate options 
	int burn_in;              // burn in chain 
	int ncycle_MCEM;          // cycle number of MCEM
	int chain_length;         // Markov chain length of each cycle 
	int sample_period;        // sampling period, recombination frequncy will be calculated, other chains discarded 
	int isGlobal;             // Gibbs sampling in overall haploMatrix, not restrict to missing obseration and partial informative genotypes  

}GibbsSAM;

// extern functions 
extern inline void multi_estimation_rf(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);

#endif
