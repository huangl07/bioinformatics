/*
 *	
 * metrohast.c : evaluate uncertainty of final maps due to viaration  of data
 *               using metropolis-hasting algorithm
 *
 */

// Copyright 2012 , Ma Chouxian <macx@biomarker.com.cn>

#include "metrohast.h"

//===================================
// static functions 
//===================================

static inline MetroHast *init_metro_hasting(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);
static inline void _metropolis_hasting(MetroHast *mh);
static inline void init_p_m_map_pos(GenoMatrix **gm, int **mappos, int *map, int maplen);
static inline double calculate_energy(GenoMatrix **gMatrix, RecLOD **redLOD, int *sOrder, int maplen);
static inline void get_neighbor_map(int *map, int maplen, gsl_rng *r);
static inline void update_plausible_position_CP(GenoMatrix **gm, int *map, int maplen, int **pmPos, int **ppp, int **mpp);
static inline void output_pp_matrix(MetroHast *mh, TSPMSET *tspm, MapOpt *opt);
static inline void free_MetroHast(MetroHast **mh);

//===================================
// interface 
//===================================

inline void generate_plausible_position_matrix(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){

	// init metropolis-hasting struct 
	MetroHast *mh = init_metro_hasting(locsINFO, tspm, opt);

	// draw map samples and generate plausible position matrix 
	_metropolis_hasting(mh);

	// output pp matrix 
	output_pp_matrix(mh, tspm, opt);

	// free 
	free_MetroHast(&mh);
}

//===================================
// utilities 
//===================================

static inline void _metropolis_hasting(MetroHast *mh){
	// get genotype matrix 
	GenoMatrix **gm = mh->gMatrix;
	// recLOD 
	RecLOD **recLOD = mh->recLOD;
	// int maplen 
	int maplen = mh->maplen;
	// current map 
	int *curMap = mh->map;
	// new Map 
	int *newMap = (int *) BioMalloc(sizeof(int) * maplen);
	// get smaple period 
	int sam_interval = mh->sample_period;

	double cur_E, new_E;
	int i, nsam = 0, total_sam = 0;

	// metropolis burn-in
	cur_E = calculate_energy(gm, recLOD, curMap, maplen);
	for (i=0; i<mh->burn_in; i++){
		memcpy(newMap, curMap, sizeof(int)*maplen);
		get_neighbor_map(newMap, maplen, mh->r);
		new_E = calculate_energy(gm, recLOD, newMap, maplen);

		if (new_E < cur_E){
			/* yay! take a step */
			memcpy(curMap, newMap, sizeof(int) * maplen);
			cur_E = new_E;
		}else if (gsl_rng_uniform(mh->r) < boltzmann(cur_E, new_E, mh->control_probability/mh->nind)){
			/* yay! take a step */
			memcpy(curMap, newMap, sizeof(int) * maplen);
			cur_E = new_E;
		}
	}
	// sampling maps and stat
	cur_E = calculate_energy(gm, recLOD, curMap, maplen);
	while (nsam < mh->n_samples){
		memcpy(newMap, curMap, sizeof(int)*maplen);
		get_neighbor_map(newMap, maplen, mh->r);
		new_E = calculate_energy(gm, recLOD, newMap, maplen);

		if (new_E < cur_E){
			/* yay! take a step */
			memcpy(curMap, newMap, sizeof(int) * maplen);
			cur_E = new_E;
		}else if (gsl_rng_uniform(mh->r) < boltzmann(cur_E, new_E, mh->control_probability/mh->nind)){
			/* yay! take a step */
			memcpy(curMap, newMap, sizeof(int) * maplen);
			cur_E = new_E;
		}

		if (total_sam % sam_interval == 0){
			update_plausible_position_CP(gm, curMap, maplen, mh->pmMapPos, mh->p_plausible_position, mh->m_plausible_position);
			nsam++;
		}
		total_sam++;
	}
	// free
	free(newMap);newMap = NULL;
}

static inline void update_plausible_position_CP(GenoMatrix **gm, int *map, int maplen, int **pmPos, int **ppp, int **mpp){
	int i, a, origin_pPos, pPos = 0, origin_mPos, mPos = 0;
	for (i=0; i<maplen; i++){
		a = map[i];
		if (strcmp(gm[a]->type, "lmxll") != 0){
			origin_mPos = pmPos[a][1];
			if (origin_mPos == -1) BioDie("%s not found in female map\n", gm[a]->marker);
			mpp[origin_mPos][mPos++]++;
		}
		if (strcmp(gm[a]->type, "nnxnp") != 0){
			origin_pPos = pmPos[a][0];
			if (origin_pPos == -1) BioDie("%s not found in male map\n", gm[a]->marker);
			ppp[origin_pPos][pPos++]++;
		}
	}
}

static inline void get_neighbor_map(int *map, int maplen, gsl_rng *r){
	int rLocs, rPos;
	while (1){
		rLocs = gsl_rng_uniform_int(r, maplen);
		rPos = gsl_rng_uniform_int(r, maplen);
		if (rLocs == rPos) continue;
		repalce_random_locus(map, rLocs, rPos);
		break;
	}
}

static inline double calculate_energy(GenoMatrix **gMatrix, RecLOD **recLOD, int *sOrder, int maplen){	
	double E = 0;
	int i, isPstart = 0, isMstart = 0, pLast = -1, mLast = -1, loci;
	for (i=0;i<maplen ;i++ ){
		loci = sOrder[i];
		if (isPstart == 0){
			if (strcmp(gMatrix[loci]->type,"nnxnp") != 0){
				pLast = loci;
				isPstart = 1;
			}
		}
		if (isMstart == 0){
			if (strcmp(gMatrix[loci]->type,"lmxll") != 0){
				mLast = loci;
				isMstart = 1;
			}
		}
		if (isPstart == 1 && pLast != loci){
			if (strcmp(gMatrix[loci]->type,"nnxnp") != 0){
				E += recLOD[loci][pLast].pRec;
				pLast = loci;
			}
		}
		if (isMstart == 1 && mLast != loci){
			if (strcmp(gMatrix[loci]->type,"lmxll") != 0){
				E += recLOD[loci][mLast].mRec;
				mLast = loci;
			}
		}
	}
	return E;
}

static inline MetroHast *init_metro_hasting(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	// metropolis hasting struct
	MetroHast *mh = (MetroHast *) BioMalloc(sizeof(MetroHast));
	// genotype matrix 
	mh->gMatrix = locsINFO->gMatrix;
	// rec.freq matrix 
	mh->recLOD = locsINFO->recLOD[0];
	// nloc 
	int nloc = locsINFO->nloc;
	// indi 
	mh->nind = (double) locsINFO->nind;

	// map 
	mh->maplen = tspm->Slen;
	mh->map = (int *) BioMalloc(sizeof(int) * mh->maplen);
	memcpy(mh->map, tspm->S, sizeof(int) * mh->maplen);
	// E
	mh->init_E = tspm->P_sarf + tspm->M_sarf;

	
	// p and m Map pos info 
	mh->pmMapPos = (int **) BioMalloc(sizeof(int *) * nloc);
	int i;
	for (i=0; i<nloc; i++){
		mh->pmMapPos[i] = (int *) BioMalloc(sizeof(int) * 2);
	}
	init_p_m_map_pos(locsINFO->gMatrix, mh->pmMapPos, mh->map, nloc);

	// gsl random number generator 
	gsl_rng_env_setup();
	mh->T = gsl_rng_default;
	mh->r = gsl_rng_alloc(mh->T);
	gsl_rng_set(mh->r,(unsigned long int)time(NULL));

	// plausible position 
	mh->plen = tspm->Plen;
	mh->mlen = tspm->Mlen;

	mh->p_plausible_position = (int **) BioMalloc(sizeof(int *) * mh->plen);
	for (i=0; i<mh->plen; i++){
		mh->p_plausible_position[i] = (int *) BioMalloc(sizeof(int) * mh->plen);
		memset(mh->p_plausible_position[i], 0, sizeof(int) * mh->plen);
	}
	mh->m_plausible_position = (int **) BioMalloc(sizeof(int *) * mh->mlen);
	for (i=0; i<mh->mlen; i++){
		mh->m_plausible_position[i] = (int *) BioMalloc(sizeof(int) * mh->mlen);
		memset(mh->m_plausible_position[i], 0, sizeof(int) * mh->mlen);
	}

	// calculate option 
	mh->control_probability = opt->PP_acceptance_control;
	mh->burn_in = opt->PP_burn_in;
	mh->n_samples = opt->PP_nsample;
	mh->sample_period = opt->PP_sample_period;

	return mh;
}

static inline void free_MetroHast(MetroHast **mh){
	free((*mh)->map); (*mh)->map = NULL;
	int i;
	for (i=0; i<(*mh)->maplen; i++) free((*mh)->pmMapPos[i]);
	for (i=0; i<(*mh)->plen; i++) free((*mh)->p_plausible_position[i]);
	for (i=0; i<(*mh)->mlen; i++) free((*mh)->m_plausible_position[i]);
	free((*mh)->pmMapPos); (*mh)->pmMapPos = NULL;
	free((*mh)->p_plausible_position); (*mh)->p_plausible_position = NULL;
	free((*mh)->m_plausible_position); (*mh)->m_plausible_position = NULL;
	gsl_rng_free((*mh)->r);

	free(*mh); *mh = NULL;
}

static inline void init_p_m_map_pos(GenoMatrix **gm, int **mappos, int *map, int maplen){
	int i, a, pPos = 0, mPos = 0;
	for (i=0; i<maplen; i++){
		a = map[i];
		if (strcmp(gm[a]->type, "abxcd") == 0 || strcmp(gm[a]->type, "efxeg") == 0 || strcmp(gm[a]->type, "hkxhk") == 0){
			mappos[a][0] = pPos;
			mappos[a][1] = mPos;
			pPos++;
			mPos++;
		}else if (strcmp(gm[a]->type, "lmxll") == 0){
			mappos[a][0] = pPos;
			mappos[a][1] = -1;
			pPos++;
		}else if (strcmp(gm[a]->type, "nnxnp") == 0){
			mappos[a][0] = -1;
			mappos[a][1] = mPos;
			mPos++;
		}
	}
}

static inline void output_pp_matrix(MetroHast *mh, TSPMSET *tspm, MapOpt *opt){
	// open file
	FILE *pfp, *mfp; 
	if( (pfp = fopen(opt->malePlausiblePosFile, "w")) == NULL )
		BioDie("open or create ppm file %s is failed!", opt->malePlausiblePosFile);
	if( (mfp = fopen(opt->femalePlausiblePosFile, "w")) == NULL )
		BioDie("open or create ppm file %s is failed!", opt->femalePlausiblePosFile);
	// get plen and mlen 
	int plen = mh->plen;
	int mlen = mh->mlen;
	// origianl p and m map 
	int *P = tspm->P;
	int *M = tspm->M;
	// genotype matrix 
	GenoMatrix **gm = mh->gMatrix;

	// get ppp and mpp
	int **ppp = mh->p_plausible_position;
	int **mpp = mh->m_plausible_position;

	// output
	int i, j;
	for (i=0; i<plen; i++){
		fprintf(pfp, "%s\t%s\t%.3f", gm[P[i]]->marker, gm[P[i]]->type, tspm->PD[i]);
		for (j=0; j<plen; j++){
			fprintf(pfp, "\t%d",  ppp[i][j]);
		}
		fprintf(pfp, "\n");
	}
	// Maternal 
	for (i=0; i<mlen; i++){
		fprintf(mfp, "%s\t%s\t%.3f", gm[M[i]]->marker, gm[M[i]]->type, tspm->MD[i]);
		for (j=0; j<mlen; j++){
			fprintf(mfp, "\t%d", mpp[i][j]);
		}
		fprintf(mfp, "\n");
	}
}


//#define METROHAST_TEST
#ifdef METROHAST_TEST

int main (void){

	return 0;
}


#endif
