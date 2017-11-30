/*
 *
 *  ripple map obtained from simulated annealing to avoid getting trapped in a local optimum 
 *
 */

// Copyright 2012, Ma Chouxian <macx@biomarker.com.cn>

/**************************************************************
	This modual is used to optimize map order obtained 
from simulated annealing.Two heuristic methods are used here:
	(1) The Initial map is digested into "blocks" by 
large gaps of genetic map. These "blocks" are rearranged 
and tested for any further improvement of objective function:
either the maximum likelihood or the sum of adjacent fractions.
	(2) Each locus is relocated to other positions to gain any 
improvement of map order. Time complexity of this operation is 
O(n^2).

****************************************************************/

#include "ripple.h"

//===================================
// declaration of static functions
//===================================


//===================================
// interface
//===================================
void ripple(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	
	// get original map and map length 
	int *map = tspm->S;
	int maplen = tspm->Slen;
	int p_sarf;
	int m_sarf;

	// relocate each locus 
	relocate_each_locus();

	// rearrange blocks
	rarrange_blocks(LocsINFO *locsINFO, int *map, int *maplen);
}

//===================================
// static functions 
//===================================
typedef struct{
	int *bestmap;  // hold the global best map, initial tspm->S
	int maplen;    // current map length

	int *isvisited;  // is a loci visited in ripple process  
//	int *isanchor;   // is a anchor loci, initialized by locsINFO->isAnchor

	int isimproved;  // whether the map changed during the whole process 

	LocsINFO *locsINFO; // info of loci 
//	RecLOD **recLOD;   // pairwise rec.freq matrix 

}Ripple;

typedef struct BLOCK{
	int pos; // block position
	int id;  // block ID 
	int nblock;  // how many blocks
	
	int *locus; // locus ID in block 
	int blen;  // block length 
	int orient; // direction of block, 1 for "+", -1 for "-" 
	
	struct BLOCK *pre;  // previous block  
	struct BLOCK *next; // next block 

}Block;



static inline void relocate_each_locus(Ripple *r){
	// get current map 
	int *curmap = r->bestmap;
	int maplen = r->maplen;
	// isvisited 
	int *isvisited = r->isvisited;

	// variables 
	int *globbest = (int *) BioMalloc(sizeof(int) * maplen); // hold global best map in mapping process 
	int cur_p_sarf = r->p_sarf;
	int cur_m_sarf = r->p_sarf;
	
	// iterately ripple 
	int i, vi, j;
	for (i=0; i<maplen; i++){
		vi = curmap[i];
		if (isvisited[vi] == 1) continue;  // if this loci has been visited in history ripple process, skip this one
		// try all other positions
		for (j=0; j<maplen; j++){
			if (j == i) continue;  // skip current position 
			/* judge whether the objective function improved, if improved, update current map and relocate next loci */
			// 
			
		}

	}

}







