#ifndef MLMAP_H
#define MLMAP_H

#include <pthread.h>
#include "biohash.h"

#define MAXIDLENGTH 128
#define MAXFIXORDER 128
#define DEBUG 0


/* define HapElem*/
// type:  abxcd & efxeg 0
//        hkxhk: 1
//        lmxll: 2
//        nnxnp: 3
// multiple: 0: only one haps; 1: two haps
// hap:	00, 01, 10, 11; -: 2
typedef struct {
    char type;                                  // maker type: for gibbs sampling use
    char isMulti;                               // if there are multiple haps
    char hapP;                                  // hap P, when isMulti equals to 2, the antithesis of hapP is 1-hapP, similarly for hapM
    char hapM;                                  // hap M
}HapElem;

/* Define genotype matrix */ 
typedef struct {
    char *marker;                               // Marker ID 
    char *type;                                 // segregation genotype 
    char *phase;                                // linkage phase 
    char **genotype;                            // individual genotypes 
}GenoMatrix;

/* define RecLOD for pairwise data */
typedef struct {
    double pRec;                                // paternal recombination frequency 
    double mRec;                                // maternal recombination frequency 
    double sRec;                                // average recombination frequency of P and M
    int pnRec;                                  // the number of paternal recombination events
    int mnRec;                                  // the number of maternal recombination events 
    double iLOD;                                // independence LOD score 
}RecLOD;

/* define LocsINFO for single linakge group */
typedef struct {
	// locus info 
    BioStrIntHASH *locs;                        // Locs => int; Hash
    GenoMatrix **gMatrix;                       // genotype matrix 
    HapElem ***hapMatrix;                       // haplotype matrix (1D loci, 2D indi, 3D HapElem pointer)
    int nloc;                                   // number of locus 
    int nind;                                   // number of individuals 
    int lmxllSum;                               // Sum number of lmxll loci 
    int nnxnpSum;                               // Sum number of nnxnp loci
    int *isAnchor;                              // indictive anchor loci array 

	// pairwise linkage information  
    RecLOD ***recLOD;                           // 2D recFreq array; recFreq for p.rec, m.rec, s.rec 
	
	// start order 
    int startOrderLen;                          // the number of loci in Start order
    int *startOrder;                            // start order

	// fix order 
    int *fixOrderLen;                           // the number of loci in fix order 
    int **fixOrder;                             // loci array of fix order 
    int fixOrderSum;                            // the number of fix order(s)
    int **isfix;                                // is a loci in fix order or not: 0->no, 1->yes

}LocsINFO;

/* define TSPMSET for map building */ 
typedef struct { 
    int round;                                  // map building round: specified by threshold of rec.freq in spatial sampling, five or six round 
    double radius;                              // threshold of rec.freq in spatial sampling, init by calOpt
    int cycle;                                  // map optimization cycle 
    int isStartSample;                          // is a sample start sample or not, 1 for start sample, 0 not 
	// currently unmapped 
    int *T;                                     // loci unmapped 
    int Tlen;                                   // T length 
	// s map 
    int *S;                                     // sexAver Map 
    int Slen;                                   // sexAver Map len 
    double *SD;                                 // sexAver map distance 
    double S_sarf;                              // not sum of P_sarf and M_sarf;
	double S_singleton;
	// p map 
    int *P;                                     // paternal map
    int Plen;                                   // paternal map len 
    double *PD;                                 // paternal map distance   
    double P_sarf;                              // sum of adjacent paternal recombination frequency
	double P_singleton;
	// m map 
    int *M;                                     // maternal map 
    int Mlen;                                   // maternal map len 
    double *MD;                                 // maternal map distance 
    double M_sarf;                              // sum of adjacent maternal recombination frequency
	double M_singleton;

}TSPMSET;


#endif
