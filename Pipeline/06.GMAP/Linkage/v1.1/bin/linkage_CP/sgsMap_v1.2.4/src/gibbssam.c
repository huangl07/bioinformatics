/*
 *
 * gibbssam.c: Gibbs sampling library 
 *
 */

// Copyright 2012, Ma Chouxian <macx@biomarker.com.cn>

#include "gibbssam.h"
#include "gibbsutil.h"

//===================================
// static functions 
//===================================

static inline void gibbs_init_state(GibbsSAM *gSam, int n);
static inline void _gibbs_sampling(GibbsSAM *gSam);
static inline void gibbs_burn_in(GibbsSAM *gSam, int nCycle);
static inline void gibbs_MCEM(GibbsSAM *gSam, int nCycle);
static inline void calculate_gibbs_sarf(GibbsSAM *gSam, double *pSarf, double *mSarf);
static inline GibbsSAM *init_gibbsSAM(LocsINFO *locsINFO,TSPMSET *tspm, MapOpt *opt );
static inline GibbsSAM *allocate_memory_gibbsSAM(LocsINFO *locsINFO,TSPMSET *tspm, MapOpt *opt );
static inline void free_gSam(GibbsSAM *gSam);
static inline void update_map(GibbsSAM *gSam, TSPMSET *tspm);
static inline void reorgnize_map_CP(GenoMatrix **gMatrix, int **pmMapPos, int *S, double *SD, double *PD, double *MD, int maplen);

//===================================
// Gibbs sampling interface 
//===================================
/*
 *  Function: multi_estimation_rf()
 *
 *  Purpose: MCEM multipoint estimation of rec.freq using gibbs sampling, interface with other moduals 
 *
 *  Args:    TSPMSET *tspm  -- map info   
 *           LocsINFO *locsINFO -- loc info 
 *           MapOpt *opt -- map calculate option 
 *
 *  Return:  None
 */

inline void multi_estimation_rf(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	// create a GibbsSAM object and init it 
	GibbsSAM *gSam = init_gibbsSAM(locsINFO, tspm, opt);
	
	// no miss observation and partial informative data 
	if (gSam->sampleLen == 0){
		fprintf(opt->logfp, "rec.freq not re-estimated due to complete data\n\n");
		// update map 
		update_map(gSam, tspm);
		// print best map info into log 
		print_best_map(gSam, opt->logfp);
	}else{ 
		// gibbs sampling process, recMatrix info will be updated 
		_gibbs_sampling(gSam);

		// debug: dump recombination frequencies matrix 
//		int i, j, a, b;
//		int *S = tspm->S;
//		int Slen = tspm->Slen;
//		for (i=0; i<Slen - 1; i++){
//			a = S[i];
//			for (j=i+1; j<Slen; j++){
//				b = S[j];
//				printf("%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%lf\t%lf\t%lf\n", i, j, a, b, locsINFO->gMatrix[a]->marker, locsINFO->gMatrix[b]->marker, locsINFO->gMatrix[a]->type, locsINFO->gMatrix[b]->type, locsINFO->recLOD[0][a][b].pRec, locsINFO->recLOD[0][a][b].mRec, locsINFO->recLOD[0][a][b].sRec);
//			}
//		}

		// update map, calculate paternal and maternal map distance and integrate map 
		update_map(gSam, tspm);

		// print log file 
		print_gibbs_info_into_log(gSam, opt->logfp);
	}

	// free 
	free_gSam(gSam);
}

//===================================
//  Gibbs sampling 
//===================================
/*
 *  Function: _gibbs_sampling()
 *
 *  Purpose: bulit-in static function  of GibbsSAM modual, multipoint estimation of rec.freq  
 *
 *  Args:    GibbsSAM *gSam  -- GibbsSAM struct   
 *
 *  Return: NONE
 */
static inline void _gibbs_sampling(GibbsSAM *gSam){
	// get ncycle
	int n = gSam->ncycle_MCEM;
	int m = n < 3 ? n : 3;
	
	// the first 3 cycles(or less than three) is used to mornitor markov chain 
	int i;
	for (i=0;i<m ;i++ ){
		// init stat 
		gibbs_init_state(gSam, i);
		// gibbs burn-in period 
		gibbs_burn_in(gSam, i);
		// mcem cycle: sampling process
		gibbs_MCEM(gSam, i);
	}
	// following cycles have no burn-in period 
	for (i=m;i<n;i++ ){
		// reset sum of rec matrix to 0
//		reset_sum_recMatrix(gSam->recMatrixSum, gSam->sMap, gSam->sLen);
		// mcem cycle: sampling process
		gibbs_MCEM(gSam, i);	
	}

	// update singleton info
	GibbsStat *gStat = gSam->mcemMornitor[n-1];
	gSam->p_singleton = gStat->p_singleton_end;
	gSam->m_singleton = gStat->m_singleton_end;
	gSam->s_singleton = (gSam->p_singleton + gSam->m_singleton) / 2;
	gSam->p_sarf = gStat->p_sarf_end;
	gSam->m_sarf = gStat->m_sarf_end;

}

//===================================
// gibbs sampling utilities  
//===================================

static inline void gibbs_init_state(GibbsSAM *gSam, int n){
	// impute missing observations
	if (n == 0){
		init_miss_observation(gSam);
	}
	// calculate total recombination counts and singleton 
	calculate_singleton_and_rn(gSam, &(gSam->mcemMornitor[n]->p_singleton_start), &(gSam->mcemMornitor[n]->m_singleton_start), &(gSam->mcemMornitor[n]->p_mrn_start), &(gSam->mcemMornitor[n]->m_mrn_start));
	// calculate sarf 
	calculate_gibbs_sarf(gSam, &(gSam->mcemMornitor[n]->p_sarf_start), &(gSam->mcemMornitor[n]->m_sarf_start));
}

static inline void gibbs_burn_in(GibbsSAM *gSam, int n){
	// get length of burn-in chain 
	int burn_in_len = gSam->burn_in;

	// iterate sampling missing data and partial informtive genotypes 
	// these samples are discarded and not used to estimate rec.freq 
	int i;
	for (i=0; i<burn_in_len; i++){
		generate_next_sample(gSam);
	}
	//  update GibbsStat info 
	calculate_singleton_and_rn(gSam, &(gSam->mcemMornitor[n]->p_singleton_burn_in), &(gSam->mcemMornitor[n]->m_singleton_burn_in), &(gSam->mcemMornitor[n]->p_mrn_burn_in), &(gSam->mcemMornitor[n]->m_mrn_burn_in));
}

static inline void gibbs_MCEM(GibbsSAM *gSam, int n){
	// get chain length 
	int chain_length = gSam->chain_length;
	// get sample interval 
	int interval = gSam->sample_period;
	// mcem 
	MCEM *mcem = gSam->mcem;

	int i, nsam = 0;
	double ps, ms, prn, mrn, ps_sum = 0, ms_sum = 0, prn_sum = 0, mrn_sum = 0;
	for (i=1; i<=chain_length; i++){
		// generate next haplotype matrix 
		generate_next_sample(gSam); 
		// sampling
		if (i % interval == 0){ 
			calculate_pairwise_rec_info(gSam); // calculate rec matrix 
			calculate_singleton_and_rn(gSam, &ps, &ms, &prn, &mrn); // calculate singleton and rn 
			ps_sum += ps; // paternal singleton ratio of sampled haplotype matrix 
			ms_sum += ms; // maternal singleton ratio of sampled haplotype matrix
			prn_sum += prn;  // paternal recombination events 
			mrn_sum += mrn;  // maternal recombination events
			nsam++;
		}
	}
	// update recMatrix 
	update_rec_matrix(mcem, gSam->sMap, gSam->sLen, nsam);
	// update sarf, singleton and rn 
	calculate_gibbs_sarf(gSam, &(gSam->mcemMornitor[n]->p_sarf_end), &(gSam->mcemMornitor[n]->m_sarf_end));
	gSam->mcemMornitor[n]->p_mrn_end = prn_sum / nsam;
	gSam->mcemMornitor[n]->m_mrn_end = mrn_sum / nsam;
	gSam->mcemMornitor[n]->p_singleton_end = ps_sum / nsam;
	gSam->mcemMornitor[n]->m_singleton_end = ms_sum / nsam;
}

/*
 *  Function: calculate_gibbs_sarf()
 *
 *  Purpose:  calculate sum of adj.rec.freq, update gSam->mcem->p_sarf and gSam->mcem->m_sarf
 *
 *  Args:     GibbsSAM *gSam  - Gibbs sampling struct
 *
 *  Return:   None
 *
 */
static inline void calculate_gibbs_sarf(GibbsSAM *gSam, double *pSarf, double *mSarf){
	// get original genotype matrix 
	GenoMatrix **gMatrix = gSam->gMatrix;
	// get mcem 
	MCEM *mcem = gSam->mcem;
	// get sexAverm map info 
	int *sMap = gSam->sMap;
	int maplen = gSam->sLen;
	// get recMatrix
	RecLOD **recMatrix = mcem->recMatrix;

	int i, isPstart = 0, isMstart = 0, pLast = -1, mLast = -1, loci;
	for (i=0;i<maplen ;i++ ){
		loci = sMap[i];
		if (isPstart == 0){
			if (strcmp(gMatrix[loci]->type, "nnxnp") != 0){
				pLast = loci;
				isPstart = 1;
			}
		}
		if (isMstart == 0){
			if (strcmp(gMatrix[loci]->type, "lmxll") != 0){
				mLast = loci;
				isMstart = 1;
			}
		}
		if (isPstart == 1 && pLast != loci){
			if (strcmp(gMatrix[loci]->type, "nnxnp") != 0){
				*pSarf += recMatrix[loci][pLast].pRec;
				pLast = loci;
			}
		}
		if (isMstart == 1 && mLast != loci){
			if (strcmp(gMatrix[loci]->type, "lmxll") != 0){
				*mSarf += recMatrix[loci][mLast].mRec;
				mLast = loci;
			}
		}
	}
}

static inline GibbsSAM *init_gibbsSAM(LocsINFO *locsINFO,TSPMSET *tspm, MapOpt *opt ){
	GibbsSAM *gSam = allocate_memory_gibbsSAM(locsINFO, tspm, opt);

	// map info
	gSam->sMap = tspm->S;  // sexAver map 
	gSam->sLen = tspm->Slen;  // sexAver map len
	gSam->SD = tspm->SD;    // sexAver distance (accumutive)
	gSam->s_singleton = 0;   // sexAver noise ratio 
	gSam->pMap = tspm->P;     // paternal map 
	gSam->pLen = tspm->Plen;  // paternal map length 
	gSam->PD = tspm->PD;    // paternal map distance 
	gSam->p_singleton = 0;    // paternal noise ratio 
	gSam->mMap = tspm->M;    // maternal map 
	gSam->mLen = tspm->Mlen;  // maternal map length 
	gSam->MD = tspm->MD;      // maternal map distance 
	gSam->m_singleton = 0;    // maternal noise ratio 

	// init P and M map info and sample space 
	init_map_and_sample_space(gSam);

	// calculate options 
	gSam->burn_in = opt->Gibbs_burn_in;
	gSam->ncycle_MCEM  = opt->Gibbs_ncycle_MCEM;
	gSam->chain_length = opt->Gibbs_chain_length;
	gSam->sample_period = opt->Gibbs_sample_period;

	return gSam;
}

static inline GibbsSAM *allocate_memory_gibbsSAM(LocsINFO *locsINFO,TSPMSET *tspm, MapOpt *opt ){
	// nind 
	int nind = locsINFO->nind;
	// get current map length 
	int Slen = tspm->Slen;
	// get n_cycle 
	int nCycle = opt->Gibbs_ncycle_MCEM;

	int i;
	GibbsSAM *gSam = (GibbsSAM *) BioMalloc(sizeof(GibbsSAM) * 1);
	// nind 
	gSam->nind = nind;
	gSam->gMatrix = locsINFO->gMatrix; // Note: don't free this 
	gSam->hapMatrix = locsINFO->hapMatrix; // haplotype matrix, do not free  	
	// sample space 
	gSam->sampleSpace = (Coordinate **) BioMalloc (sizeof(Coordinate*) * (Slen * nind));
	for (i=0; i<Slen * nind; i++){
		gSam->sampleSpace[i] = (Coordinate *) BioMalloc (sizeof(Coordinate));
		gSam->sampleSpace[i]->pBinLoc = (int *) BioMalloc(sizeof(int) * Slen);
		gSam->sampleSpace[i]->mBinLoc = (int *) BioMalloc(sizeof(int) * Slen);
		memset(gSam->sampleSpace[i]->pBinLoc, -1, sizeof(int) * Slen);
		memset(gSam->sampleSpace[i]->mBinLoc, -1, sizeof(int) * Slen);
		gSam->sampleSpace[i]->pBinLen = 0;
		gSam->sampleSpace[i]->mBinLen = 0;
	}
	gSam->sampleLen = 0;
	//mcem cycle 
	gSam->mcem = (MCEM *) BioMalloc(sizeof(MCEM) * 1); 
	gSam->mcem->hapMatrix = gSam->hapMatrix;  // haplotype matrix 
	gSam->mcem->recMatrix = locsINFO->recLOD[0]; // recombination frequency matrix 
	gSam->mcem->recMatrixSum = locsINFO->recLOD[1];  // sum of sampled rec matrix 
	// gibbs statistic: to morniter the convergence of  markov chain
	gSam->mcemMornitor = (GibbsStat **) BioMalloc (sizeof(GibbsStat *) * nCycle); 
	for (i=0; i<nCycle; i++){
		gSam->mcemMornitor[i] = (GibbsStat *) BioMalloc(sizeof(GibbsStat) * 1);
	}
	// init gibbs random number generator 
	gsl_rng_env_setup();
	gSam->T = gsl_rng_default;
	gSam->r = gsl_rng_alloc(gSam->T);
	gsl_rng_set(gSam->r,(unsigned long int)time(NULL));
	// map pos 
	gSam->pmMapPos = (int **) BioMalloc(sizeof(int *) * Slen); 
	for (i=0;i<Slen ;i++ ){
		gSam->pmMapPos[i] = (int *) BioMalloc(sizeof(int) * 2);
	}

	return gSam;
}

static inline void free_gSam(GibbsSAM *gSam){
	gSam->gMatrix = NULL;
	gSam->hapMatrix = NULL;
	int i;
	for (i=0; i<gSam->sampleLen; i++){
		free(gSam->sampleSpace[i]->pBinLoc);gSam->sampleSpace[i]->pBinLoc = NULL;
		free(gSam->sampleSpace[i]->mBinLoc);gSam->sampleSpace[i]->mBinLoc = NULL;
		free(gSam->sampleSpace[i]); gSam->sampleSpace[i] = NULL;
	}
	free(gSam->sampleSpace); gSam->sampleSpace = NULL;
	// mcem 
	gSam->mcem->hapMatrix = NULL;
	gSam->mcem->recMatrix = NULL;
	gSam->mcem->recMatrixSum = NULL;
	free(gSam->mcem); gSam->mcem = NULL;
	// GibbsStat
	for (i=0; i<gSam->ncycle_MCEM; i++){
		free(gSam->mcemMornitor[i]); gSam->mcemMornitor[i] = NULL;
	}
	free(gSam->mcemMornitor); gSam->mcemMornitor = NULL;
	// r 
	gsl_rng_free(gSam->r);

	// map pos info 
	for (i=0; i<gSam->sLen; i++){
		free(gSam->pmMapPos[i]);gSam->pmMapPos[i] = NULL;
	}
	free(gSam->pmMapPos);gSam->pmMapPos = NULL;
	// map and map dist 
	gSam->sMap = NULL;
	gSam->pMap = NULL;
	gSam->mMap = NULL;
	gSam->SD = NULL;
	gSam->PD = NULL;
	gSam->MD = NULL;
	// gSam 
	free(gSam); gSam = NULL;
}

static inline void update_map(GibbsSAM *gSam, TSPMSET *tspm){
	// get S and sLen 
	int *S = gSam->sMap;
	double *SD = gSam->SD;
	int sLen = gSam->sLen;
	int *P = gSam->pMap;
	int pLen = gSam->pLen;
	double *PD = gSam->PD;
	int *M = gSam->mMap;
	int mLen = gSam->mLen;
	double *MD = gSam->MD;
	// get gMatrix and recLOD 
	GenoMatrix **gMatrix = gSam->gMatrix;
	RecLOD **recLOD = gSam->mcem->recMatrix;

	// calculate paternal and maternal map distance
	int i;
	double d;
	PD[0] = 0;
	for (i = 1;i < pLen ; i++){
		d = haldane(recLOD[P[i]][P[i-1]].pRec);
		PD[i] = PD[i-1] + d;
	}

	MD[0] = 0;
	for (i = 1;i < mLen ; i++){
		d = haldane(recLOD[M[i]][M[i-1]].mRec);
		MD[i] = MD[i-1] + d;
	}

	// reorgnize sMap 
	reorgnize_map_CP(gMatrix, gSam->pmMapPos, S, SD, PD, MD, sLen);
	// update singleton
	tspm->Slen = sLen;
	tspm->Plen = pLen;
	tspm->Mlen = mLen;
	tspm->P_singleton = gSam->s_singleton;
	tspm->M_singleton = gSam->m_singleton;
	tspm->S_singleton = (tspm->P_singleton + tspm->M_singleton) / 2;
	tspm->P_sarf = gSam->p_sarf;
	tspm->M_sarf = gSam->m_sarf;
	tspm->S_sarf = 0;
	for (i = 1;i < sLen ; i++){
		tspm->S_sarf += inverseHaldane(SD[i] - SD[i-1]);
	}
}
/*
 *  Function: reorgnize_map_CP()
 *
 *  Purpose:  integrate Paternal and maternal maps 
 *
 *  Args:     GenoMatrix **gMatrix -- genotype matrix 
 *            RecLOD **recLOD -- pairwise rec.freq 
 *            int *P -- paternal map 
 *            int *M -- maternal map 
 *            int *S -- integrated map 
 *            int meplen -- integrated map length
 *            double *SD -- integrated map distance 
 *
 *  Return:   None
 *
 */
static inline void reorgnize_map_CP(GenoMatrix **gMatrix, int **pmMapPos, int *S, double *SD, double *PD, double *MD, int maplen){
	// anchoring loci
	int *order = (int *) BioMalloc (sizeof(int) * maplen);
	// get anchor loci
	int i, j, anchorNum = 0, lmxllNum = 0, nnxnpNum = 0;
	for (i = 0;i < maplen ; i++){
		if (strcmp(gMatrix[S[i]]->type, "lmxll") == 0 ){
			lmxllNum++;
		}else if (strcmp(gMatrix[S[i]]->type, "nnxnp") == 0){
			nnxnpNum++;
		}else{
			order[anchorNum] = i;
			anchorNum++;
		}
//		if (strcmp(gMatrix[S[i]]->type, "lmxll") != 0 && strcmp(gMatrix[S[i]]->type, "nnxnp") != 0) 
//			order[anchorNum++] = i;
	}
//	if (anchorNum == 0) BioDie("can't find anchor loci during integrition of p and m maps");
	
	if (anchorNum == 0){ // no anchor loci
		if (lmxllNum != 0 && nnxnpNum != 0){
			BioDie("can't find anchor loci during integrition of p and m maps");
		}else if (nnxnpNum != 0){
			memcpy(SD, MD, sizeof(double) * nnxnpNum);
		}else if (lmxllNum != 0){
			memcpy(SD, PD, sizeof(double) * lmxllNum);
		}

		free(order);order = NULL;
		return;
	}

	/* estimate map distance using interpolation and extrapolation algorithm */
	// the outermost anchor markers: head --- extrapolation
	int pa, pb, ppa, mpa, ppb, mpb;
	double d1, d2, d;
	pa = order[0]; // pos of first anchor loci
	ppa = pmMapPos[pa][0];
	mpa = pmMapPos[pa][1];
	d1 = PD[ppa];
	d2 = MD[mpa];
	
	SD[pa] = (d1 > d2) ? d1: d2;
	for (i = 0;i <pa ; i++){ // lmxll and nnxnp loci position 
		if (strcmp(gMatrix[S[i]]->type, "lmxll") == 0){ // p loci
			ppa = pmMapPos[i][0];
			SD[i] = SD[pa] - (d1 - PD[ppa]);
		}else if (strcmp(gMatrix[S[i]]->type, "nnxnp") == 0){ // m loci 
			mpa = pmMapPos[i][1];
			SD[i] = SD[pa] - (d2 - MD[mpa]);
		}else{
			BioDie("error in map integration");
		}
	}
	
	// middle
	double delta_dist;
	for (i = 1;i <anchorNum  ; i++){
		pa = order[i-1];
		pb = order[i];
		ppa = pmMapPos[pa][0];
		mpa = pmMapPos[pa][1];
		ppb = pmMapPos[pb][0];
		mpb = pmMapPos[pb][1];
		// dist of anchor
		d1 = PD[ppb] - PD[ppa];
		d2 = MD[mpb] - MD[mpa];
		delta_dist = (d1 + d2) / 2; // averaged p and m distance 
		SD[pb] = SD[pa] + delta_dist;
		// dist of p and m loci 
		for (j = pa+1;j < pb ; j++){
			if (strcmp(gMatrix[S[j]]->type, "lmxll") == 0){ // p loci
				ppb = pmMapPos[j][0];
				d = PD[ppb] - PD[ppa];
				SD[j] = d1 == 0 ? SD[pa] : SD[pa] + delta_dist * (d / d1); // interpolation
			}else if (strcmp(gMatrix[S[j]]->type, "nnxnp") == 0){ // m loci 
				mpb = pmMapPos[j][1];
				d = MD[mpb] - MD[mpa];
				SD[j] = d2 == 0 ? SD[pa] : SD[pa] + delta_dist * (d / d2); // interpolation
			}else{
				BioDie("error in map integration");
			}
		}
	}
	// the outermost anchor markers: tail --- extrapolation
	pb = order[anchorNum-1];
	ppb = pmMapPos[pb][0];
	mpb = pmMapPos[pb][1];
	d1 = PD[ppb];
	d2 = MD[mpb];
	if (pb != maplen - 1){ // not the last marker 
		for (i = pb + 1;i < maplen ; i++){
			if (strcmp(gMatrix[S[i]]->type, "lmxll") == 0){ // p loci
				ppb = pmMapPos[i][0];
				SD[i] = SD[pb] + (PD[ppb] - d1); // extrapolation
			}else if (strcmp(gMatrix[S[i]]->type, "nnxnp") == 0){ // m loci 
				mpb = pmMapPos[i][1];
				SD[i] = SD[pb] + (MD[mpb] - d2); // extrapolation
			}else{
				BioDie("error in map integration");
			}
		}
	}
	// order loci according to map dist obtained above, actually determine order between lmxll and nnxnp
	int swap;
	double swapDist;
//	printf(">\n");
	for (i = 0;i < maplen - 1 ; i++){
		for (j = i+1;j < maplen ; j++){
//			printf("%s\t%s\t%lf\n", gMatrix[S[i]]->type, gMatrix[S[j]]->type, SD[i] - SD[j]);
			if (SD[i] > SD[j]){ // swap loci and dist 	
				swap = S[i]; S[i] = S[j]; S[j] = swap;
				swapDist = SD[i]; SD[i] = SD[j]; SD[j] = swapDist;
			}
		}
	}
//	printf("== after changed\n");
//	for (i = 0;i < maplen - 1 ; i++){
//		for (j = i+1;j < maplen ; j++){
//			printf("%s\t%s\t%lf\n", gMatrix[S[i]]->type, gMatrix[S[j]]->type, SD[i] - SD[j]);
//		}
//	}
	// free 
	free(order);order = NULL;
}

//#define GIBBSSAM_TEST
#ifdef GIBBSSAM_TEST

int main(int argc, char *argv[]){

	return 0;
}

#endif
