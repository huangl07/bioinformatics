/*
 *
 * gibbsutil.c: utility for gibbs sampling 
 *
 */
// Copyright 2012, Ma Chouxian <macx@biomarker.com.cn>

#include "gibbsutil.h"

inline void print_gibbs_info_into_log(GibbsSAM *gSam, FILE *fp){
	// printing gibbs stat 
	print_gibbs_stat(gSam, fp);
	// printing best map 
	print_best_map(gSam, fp);
}

inline void print_gibbs_stat(GibbsSAM *gSam, FILE *fp){
	// get gibbs stat 
	GibbsStat **gStat = gSam->mcemMornitor;
	int n = gSam->ncycle_MCEM;
	int m = n > 3 ? 3 : n;

	fprintf(fp, "multipoint estimation of rec. frequencies:\n");
	fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%-10s  %-10s\n", "posotion", "p.sarf", "m.sarf", "t.sarf", "p.mrn", "m.mrn", "t.mrn", "p.sigton", "m.sigton");
	fprintf(fp, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%-10s  %-10s\n", "--------", "------", "------", "------", "-----", "-----", "-----", "--------", "--------");
	// at start 
	int i;
	for (i=0; i<m; i++){
		fprintf(fp, "at start %d\t%.3f\t%.3f\t%.3f\t%.0f\t%.0f\t%.0f\t%-10.3f  %-10.3f\n", i+1, gStat[i]->p_sarf_start, gStat[i]->m_sarf_start, gStat[i]->p_sarf_start+gStat[i]->m_sarf_start, gStat[i]->p_mrn_start, gStat[i]->m_mrn_start, gStat[i]->p_mrn_start+gStat[i]->m_mrn_start, gStat[i]->p_singleton_start, gStat[i]->m_singleton_start);
	}
	// after burn in 
	for (i=0; i<m; i++){
		fprintf(fp, "after burn-in %d\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%-10.3f  %-10.3f\n", i+1, "", "", "", gStat[i]->p_mrn_burn_in, gStat[i]->m_mrn_burn_in, gStat[i]->p_mrn_burn_in+gStat[i]->m_mrn_burn_in, gStat[i]->p_singleton_burn_in, gStat[i]->m_singleton_burn_in);
	}
	// at end 
	for (i=0; i<n; i++){
		fprintf(fp, "after cycle %d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%-10.3f  %-10.3f\n", i+1, gStat[i]->p_sarf_end, gStat[i]->m_sarf_end, gStat[i]->p_sarf_end+gStat[i]->m_sarf_end, gStat[i]->p_mrn_end, gStat[i]->m_mrn_end, gStat[i]->p_mrn_end+gStat[i]->m_mrn_end, gStat[i]->p_singleton_end, gStat[i]->m_singleton_end);
	}

	fprintf(fp,"\n");
}
inline void print_best_map(GibbsSAM *gSam, FILE *fp){
	// get gMatrix 
	GenoMatrix **gm = gSam->gMatrix;
	// get mapPos 
	int *sMap = gSam->sMap;
	int slen = gSam->sLen;
	
	double *SD = gSam->SD;
	double *PD = gSam->PD;
	double *MD = gSam->MD;
	
	fprintf(fp, "best map order:\n");
	fprintf(fp, "%4s/%-4s\t%-15s\t%-8s   %-10s\t%-10s\t%-10s\n", "Nr.", "pos.", "locus", "seg.type", "paternal", "maternal", "sex-average");
	fprintf(fp, "%4s/%-4s\t%-15s\t%-8s   %-10s\t%-10s\t%-10s\n", "---", "----", "-----", "--------", "--------", "--------", "-----------");
	int i, a, pPos = 0, mPos = 0;
	for (i=0; i<slen; i++){
		a = sMap[i];
		if (strcmp(gm[a]->type, "lmxll") == 0){
			fprintf(fp, "%4d/%-4d\t%-15s\t<%s>    %-10.3f\t%-10s\t%-10.3f\n", a, i, gm[a]->marker, gm[a]->type, PD[pPos]*100, " ", SD[i]*100);
			pPos++;
		}else if (strcmp(gm[a]->type, "nnxnp") == 0){
			fprintf(fp, "%4d/%-4d\t%-15s\t<%s>    %-10s\t%-10.3f\t%-10.3f\n", a, i, gm[a]->marker, gm[a]->type, " ", MD[mPos]*100, SD[i]*100);
			mPos++;
		}else if (strcmp(gm[a]->type, "abxcd") == 0 || strcmp(gm[a]->type, "efxeg") == 0 || strcmp(gm[a]->type, "hkxhk") == 0){
			fprintf(fp, "%4d/%-4d\t%-15s\t<%s>+   %-10.3f\t%-10.3f\t%-10.3f\n", a, i, gm[a]->marker, gm[a]->type, PD[pPos]*100, MD[mPos]*100, SD[i]*100);
			pPos++;
			mPos++;
		}else{
			BioDie("type error: %s\t%s", gm[a]->marker, gm[a]->type);
		}
	}
	fprintf(fp, "\n");
}

/*
 *  Function: calculate_singleton_and_rn()
 *
 *  Purpose:  calculate recombination counts and singlton ratio 
 *
 *  Args:	  GibbsSAM *gSam, double *ps, double *ms, int *pn, int *mn
 *
 *  Return:   None
 */
inline void calculate_singleton_and_rn(GibbsSAM *gSam, double *ps, double *ms, double *pn, double *mn){
	// get p and m Map 
	int *P = gSam->pMap;
	int *M = gSam->mMap;
	int plen = gSam->pLen;
	int mlen = gSam->mLen;

	int n = gSam->nind;

	// get mcem and hap matrix 
	MCEM *mcem = gSam->mcem;
	HapElem ***hm = mcem->hapMatrix;
	
	int i, j, pSingleton = 0, mSingleton = 0;
	int fa, a, na;
	double psr = 0, msr = 0, prn = 0, mrn = 0; 
	// P 
	for (i=1; i<plen - 1; i++){
		a = P[i];
		fa = P[i-1];
		na = P[i+1];	
		for (j=0; j<n; j++){
			if (hm[a][j]->hapP ^ hm[fa][j]->hapP) prn++;
			if ((hm[a][j]->hapP ^ hm[fa][j]->hapP) && (hm[a][j]->hapP ^ hm[na][j]->hapP)) pSingleton++;
		}
	}
	// M 
	for (i=1; i<mlen - 1; i++){
		a = M[i];
		fa = M[i-1];
		na = M[i+1];	
		for (j=0; j<n; j++){
			if (hm[a][j]->hapM ^ hm[fa][j]->hapM) mrn++;
			if ((hm[a][j]->hapM ^ hm[fa][j]->hapM) && (hm[a][j]->hapM ^ hm[na][j]->hapM)) mSingleton++;
		}
	}
	
	psr = pSingleton / (double)(plen * n);
	msr = mSingleton / (double)(mlen * n);
	// result 
	*ps = psr;
	*ms = msr;
	*pn = prn;
	*mn = mrn;
}

/*
 *  Function: update_rec_matrix()
 *
 *  Purpose:  at the and of markov chain, calculate average rec info and update recMatrix 
 *
 *  Args:     
 *
 *  Return:   None
 *
 */
inline void update_rec_matrix(MCEM *mcem, int *map, int maplen, int n){
	// get recmatrix and recMatrix  
	RecLOD **recMatrix = mcem->recMatrix; 
	RecLOD **recMatrixSum = mcem->recMatrixSum; 

	// iterate update 
	int i, j, a, b;
	for (i=0; i<maplen-1; i++){
		a = map[i];
		for (j=i+1; j<maplen; j++){
			b = map[j];
			recMatrix[a][b].pRec = recMatrixSum[a][b].pRec / (double)n;
			recMatrix[b][a].pRec = recMatrix[a][b].pRec;
			recMatrix[a][b].mRec = recMatrixSum[a][b].mRec / (double)n;
			recMatrix[b][a].mRec = recMatrix[a][b].mRec;

			recMatrixSum[a][b].pRec = 0;
			recMatrixSum[b][a].pRec = 0;
			recMatrixSum[a][b].mRec = 0;
			recMatrixSum[b][a].mRec = 0;
		}
	}
}
/*
 *  Function: calculate_pairwise_rec_info()
 *
 *  Purpose:  calculate pairwise rec.freq at sample interval and update recMatrixSum
 *
 *  Args:     
 *
 *  Return:   None
 *
 */

inline void calculate_pairwise_rec_info(GibbsSAM *gSam){
	// genotype matrix 
	GenoMatrix **gMatrix = gSam->gMatrix;
	// recMatrix and recMatrixSum
	RecLOD **recMatrix = gSam->mcem->recMatrix;
	RecLOD **recMatrixSum = gSam->mcem->recMatrixSum;
	// haplo Matrix 
	HapElem ***hm = gSam->mcem->hapMatrix;
	// get sMap and sMap len 
	int *S = gSam->sMap;
	int maplen = gSam->sLen;
	// nind 
	double n = (double) gSam->nind;
	double n2 = n / 2.0;

	// iter calculate 
	int i, j, k, a, b;
	double pRecSum, mRecSum;
	for (i=0; i<maplen - 1; i++){
		a = S[i];
		for (j=i+1; j<maplen; j++){
			b = S[j];
			pRecSum = 0;
			mRecSum = 0;
			if ((strcmp(gMatrix[a]->type,"lmxll") == 0 && strcmp(gMatrix[b]->type,"nnxnp") != 0)
				|| (strcmp(gMatrix[b]->type,"lmxll") == 0 && strcmp(gMatrix[a]->type,"nnxnp") != 0)){ // lmxll series 
				// stat paternal recombinants 
				for (k=0; k<n; k++){
					pRecSum += hm[a][k]->hapP ^ hm[b][k]->hapP;
				}
				// update recMatrix and recMatrix sum
				recMatrix[a][b].pnRec = pRecSum;
				recMatrix[b][a].pnRec = pRecSum;
				recMatrix[a][b].pRec = pRecSum < n2 ? pRecSum / n : 0.4999;
				recMatrix[b][a].pRec = pRecSum < n2 ? pRecSum / n : 0.4999;
				
				recMatrixSum[a][b].pRec += recMatrix[a][b].pRec;
				recMatrixSum[b][a].pRec = recMatrixSum[a][b].pRec;
				recMatrixSum[a][b].mRec += 1; // means record not exists 
				recMatrixSum[b][a].mRec = recMatrixSum[a][b].mRec;
			}else if ((strcmp(gMatrix[a]->type,"nnxnp") == 0 && strcmp(gMatrix[b]->type,"lmxll") != 0)
				|| (strcmp(gMatrix[b]->type,"nnxnp") == 0 && strcmp(gMatrix[a]->type,"lmxll") != 0)){ // nnxnp series 
				// stat maternal recombinants 
				for (k=0; k<n; k++){
					mRecSum += hm[a][k]->hapM ^ hm[b][k]->hapM;
				}
				// update recMatrix and recMatrix sum
				recMatrix[a][b].mnRec = mRecSum;
				recMatrix[b][a].mnRec = mRecSum;
				recMatrix[a][b].mRec = mRecSum < n2 ? mRecSum / n : 0.4999;
				recMatrix[b][a].mRec = mRecSum < n2 ? mRecSum / n : 0.4999;
				
				recMatrixSum[a][b].mRec += recMatrix[a][b].mRec;
				recMatrixSum[b][a].mRec = recMatrixSum[a][b].mRec;
				recMatrixSum[a][b].pRec += 1;
				recMatrixSum[b][a].pRec = recMatrixSum[a][b].pRec;
			}else if (strcmp(gMatrix[a]->type,"lmxll") != 0 
				&& strcmp(gMatrix[a]->type,"nnxnp") != 0
				&& strcmp(gMatrix[b]->type,"lmxll") != 0
				&& strcmp(gMatrix[b]->type,"nnxnp") != 0){ // between anchor loci and anchor loci
				// stat paternal and maternal recombinants 
				for (k=0; k<n; k++){
					pRecSum += hm[a][k]->hapP ^ hm[b][k]->hapP;
					mRecSum += hm[a][k]->hapM ^ hm[b][k]->hapM;
				}
				// update recMatrix and recMatrix sum
				recMatrix[a][b].pnRec = pRecSum;
				recMatrix[b][a].pnRec = pRecSum;
				recMatrix[a][b].pRec = pRecSum < n2 ? pRecSum / n : 0.4999;
				recMatrix[b][a].pRec = pRecSum < n2 ? pRecSum / n : 0.4999;
				recMatrix[a][b].mnRec = mRecSum;
				recMatrix[b][a].mnRec = mRecSum;
				recMatrix[a][b].mRec = mRecSum < n2 ? mRecSum / n : 0.4999;
				recMatrix[b][a].mRec = mRecSum < n2 ? mRecSum / n : 0.4999;
				
				recMatrixSum[a][b].pRec += recMatrix[a][b].pRec;
				recMatrixSum[b][a].pRec = recMatrixSum[a][b].pRec;
				recMatrixSum[a][b].mRec += recMatrix[a][b].mRec;
				recMatrixSum[b][a].mRec = recMatrixSum[a][b].mRec;
			}
		}
	}	
}
/*
 *  Function: generate_next_sample()
 *
 *  Purpose:  generate next haplotype matrix using gibbs sampler, conditioned on map order, map dist and haplotype of adjacent loci 
 *
 *  Args:     GibbsSAM *gSam
 *
 *  Return:   None
 *
 */

inline void generate_next_sample(GibbsSAM *gSam){
	// get sample space info 
	Coordinate **ss = gSam->sampleSpace;
	int sl = gSam->sampleLen;
	if (sl == 0) return;  // sample space is empty, means no missing observation and partial informative loci 

	// shuffling sample space  
	int *sIndex = (int *) BioMalloc(sizeof(int) * sl);
	int i;
	for (i=0;i<sl ;i++ ){
		sIndex[i] = i;
	}
	gsl_ran_shuffle(gSam->r, sIndex, sl, sizeof(int)); // shufling 
	// sampling missing observation and partial informative genotypes conditioned on paternal and maternal maps and haplotypes of adjacent loci   
	for (i=0;i<sl ;i++ ){
		sampling_one_observation(gSam,ss[sIndex[i]]);
	}
	// free 
	free(sIndex);sIndex=NULL;
}
/*
 *  Function: sampling_one_observation()
 *
 *  Purpose:  sampling one element of haplotype matrix conditioned on paternal and maternal maps and haplotypes of adjacent loci
 *
 *  Args:     GibbsSAM *gSam, Coordinate *ss
 *
 *  Return:   None
 *
 */

inline void sampling_one_observation(GibbsSAM *gSam, Coordinate *ss){
	// sampling missing observation or partial informative 
	if (ss->isMulti == 0){ // missing
		sampling_miss_observation(gSam, ss);
	}else if (ss->isMulti == 2){ // partial informative 
		sampling_partial_informative(gSam, ss);
	}else{
		BioDie("Gibbs sampling error: not missing observation or partial informative");
	}
}
/*
 *  Function: sampling_partial_informative()
 *
 *  Purpose:  sampling hk obervation of hkxhk loci to eliminate ambiguity of parental allele contribution   
 *
 *  Args:    GibbsSAM *gSam, Coordinate *ss 
 *
 *  Return:   None
 *
 */
inline void sampling_partial_informative(GibbsSAM *gSam, Coordinate *ss){
	// MCEM cycle 
	MCEM *mcem = gSam->mcem;
	// haplo matrix 
	HapElem ***hm = mcem->hapMatrix;
	// get recMatrix 
	RecLOD **rm = mcem->recMatrix;
	// get p and m map len 
	int pLen = gSam->pLen;
	int mLen = gSam->mLen;
	// get indi and loci 
	int indi = ss->indi;
	int loci = ss->loci;
	// sample space for hk
	char sp[2][2];
	double p[2] = {0,0};
	sp[0][0] = hm[loci][indi]->hapP;
	sp[0][1] = hm[loci][indi]->hapM;
	sp[1][0] = 1 - hm[loci][indi]->hapP;
	sp[1][1] = 1 - hm[loci][indi]->hapM;
	// calculate probability 
	double pProb, mProb;
	int i;
	for (i=0; i<2; i++){
		pProb = calculate_paternal_condition_prob(hm, rm, ss, sp[i][0], pLen, gSam->r);
		mProb = calculate_maternal_condition_prob(hm, rm, ss, sp[i][1], mLen, gSam->r);
		p[i] = pProb * mProb;
//		printf("p1 = %lf, m1 = %lf, p[%d] = %lf\n",pProb, mProb, i, p[i]);
	}
	// standarlization: conditonal on observation of hk
	p[0] = p[0]/(p[0] + p[1]);
	// sampling 
	if (gsl_rng_uniform(gSam->r) >= p[0]){ // set to antithesis haplo of current haplo
		ss->pSamHap = sp[1][0];
		ss->mSamHap = sp[1][1];
		hm[loci][indi]->hapP = sp[1][0];
		hm[loci][indi]->hapM = sp[1][1];
	}else{
		ss->pSamHap = sp[0][0];
		ss->mSamHap = sp[0][1];
	}
	// update bin marker
	int *pBin = ss->pBinLoc;
	int *mBin = ss->mBinLoc;
	int a, n = ss->pBinLen, m = ss->mBinLen;
	for (i=0; i<n; i++){
		a = pBin[i];
		if (hm[a][indi]->isMulti != 2){
			hm[a][indi]->hapP = ss->pSamHap;
		}else if (hm[a][indi]->hapP != ss->pSamHap){
			hm[a][indi]->hapP = ss->pSamHap;
			hm[a][indi]->hapM = 1 - hm[a][indi]->hapM; // change to 
		}
	}
	for (i=0; i<m; i++){
		a = mBin[i];
		if (hm[a][indi]->isMulti != 2){
			hm[a][indi]->hapM = ss->mSamHap;
		}else if (hm[a][indi]->hapM != ss->mSamHap){
			hm[a][indi]->hapM = ss->mSamHap;
			hm[a][indi]->hapP = 1 - hm[a][indi]->hapP;
		}
	}

}

/*
 *  Function: sampling_miss_observation()
 *
 *  Purpose:  sampling missing observations 
 *
 *  Args:    GibbsSAM *gSam, Coordinate *ss 
 *
 *  Return:   None
 *
 */

inline void sampling_miss_observation(GibbsSAM *gSam, Coordinate *ss){
	// MCEM cycle 
	MCEM *mcem = gSam->mcem;
	// haplo matrix 
	HapElem ***hm = mcem->hapMatrix;
	// get recMatrix 
	RecLOD **rm = mcem->recMatrix;
	// get p and m map len 
	int pLen = gSam->pLen;
	int mLen = gSam->mLen;
	// get type
	int type = ss->type;

	// sampling: different cases 
	if (type == 2){ // lmxll loci 
		sample_paternal_haplo(hm, rm, ss, pLen, gSam->r);
	}else if (type == 3){ // nnxnp loci
		sample_maternal_haplo(hm, rm, ss, mLen, gSam->r);
	}else{ // other case 
		sample_anchor_haplo(hm, rm, ss, pLen, mLen, gSam->r); 
	}
}
/*
 *  Function: sample_anchor_haplo()
 *
 *  Purpose:  sampling missing observation of abxcd,efxeg and hkxhk 
 *
 *  Args:     ...
 *
 *  Return:   None
 *
 */
inline void sample_anchor_haplo(HapElem ***hm, RecLOD **rm, Coordinate *ss, int pLen, int mLen, gsl_rng *r){
	// get indi id and loci
	int indi = ss->indi;
	int loci = ss->loci;
	
	// sample space 
	char sp[4][2] = {{0,0},{0,1},{1,0},{1,1}};
	double p[] = {0,0,0,0};
	// calculate joint probability  
	int i;
	double pProb, mProb, u;
	for (i=0; i<4; i++){
		pProb = calculate_paternal_condition_prob(hm, rm, ss, sp[i][0], pLen, r);
		mProb = calculate_maternal_condition_prob(hm, rm, ss, sp[i][1], mLen, r);
		p[i] = pProb * mProb;
	}
	// sample from {00,01,10,11} 
	for (i=1; i<4; i++) p[i] += p[i-1];  // accumutive probability 
	u = gsl_rng_uniform(r);
	if (u < p[0]){
		ss->pSamHap = 0;
		ss->mSamHap = 0;
	}else if (u >= p[0] && u < p[1]){
		ss->pSamHap = 0;
		ss->mSamHap = 1;
	
	}else if (u >= p[1] && u < p[2]){
		ss->pSamHap = 1;
		ss->mSamHap = 0;
	}else if (u >= p[2] && u < p[3]){
		ss->pSamHap = 1;
		ss->mSamHap = 1;
	}
	// update haplo matrix
	hm[loci][indi]->hapP = ss->pSamHap;
	hm[loci][indi]->hapM = ss->mSamHap ;
	// update bin markers 
	int *pBin = ss->pBinLoc;
	int *mBin = ss->mBinLoc;
	int a, n = ss->pBinLen, m = ss->mBinLen;
	for (i=0; i<n; i++){
		a = pBin[i];
		if (hm[a][indi]->isMulti != 2){
			hm[a][indi]->hapP = ss->pSamHap;
		}else if (hm[a][indi]->hapP != ss->pSamHap){
			hm[a][indi]->hapP = ss->pSamHap;
			hm[a][indi]->hapM = 1 - hm[a][indi]->hapM;
		}	
	}
	for (i=0; i<m; i++){
		a = mBin[i];
		if (hm[a][indi]->isMulti != 2){
			hm[a][indi]->hapM = ss->mSamHap;
		}else if (hm[a][indi]->hapM != ss->mSamHap){
			hm[a][indi]->hapM = ss->mSamHap;
			hm[a][indi]->hapP = 1 - hm[a][indi]->hapP;
		}
	}
}
/*
 *  Function: sample_anchor_haplo()
 *
 *  Purpose:  sampling missing observation of lmxll 
 *
 *  Args:     ...
 *
 *  Return:   None
 *
 */
inline void sample_paternal_haplo(HapElem ***hm, RecLOD **rm, Coordinate *ss, int pLen, gsl_rng *r){
	// get indi and loci 
	int indi = ss->indi;
	int loci = ss->loci;
	// calculate probability  
	double p = calculate_paternal_condition_prob(hm, rm, ss, 0, pLen, r);
	// update haplo matrix
	ss->pSamHap = (gsl_rng_uniform(r) < p) ? 0 : 1 ;
	hm[loci][indi]->hapP = ss->pSamHap;
	// update bin marker
	int *pBin = ss->pBinLoc;
	int i, a, n = ss->pBinLen;
	for (i=0; i<n; i++){
		a = pBin[i];
		if (hm[a][indi]->isMulti != 2){
			hm[a][indi]->hapP = ss->pSamHap;
		}else if (hm[a][indi]->hapP != ss->pSamHap){
			hm[a][indi]->hapP = ss->pSamHap;
			hm[a][indi]->hapM = 1 - hm[a][indi]->hapM;
		}
	}
}
/*
 *  Function: sample_anchor_haplo()
 *
 *  Purpose:  sampling missing observation of nnxnp 
 *
 *  Args:     ...
 *
 *  Return:   None
 *
 */
inline void sample_maternal_haplo(HapElem ***hm, RecLOD **rm, Coordinate *ss, int mLen, gsl_rng *r){
	// get indi and loci 
	int indi = ss->indi;
	int loci = ss->loci;
	// calculate probability 
	double p = calculate_maternal_condition_prob(hm, rm, ss, 0, mLen, r);
	// update haplo matrix
	ss->mSamHap = (gsl_rng_uniform(r) < p) ? 0 : 1 ;
	hm[loci][indi]->hapM = ss->mSamHap ;
	// update bin marker
	int *mBin = ss->mBinLoc;
	int i, a, n = ss->mBinLen;
	for (i=0; i<n; i++){
		a = mBin[i];
		if (hm[a][indi]->isMulti != 2){
			hm[a][indi]->hapM = ss->mSamHap;
		}else if (hm[a][indi]->hapM != ss->mSamHap){
			hm[a][indi]->hapM = ss->mSamHap;
			hm[a][indi]->hapP = 1 - hm[a][indi]->hapP;
		}	
	}
}

/*
 *  Function: calculate_paternal_condition_prob()
 *
 *  Purpose:  calculate conditional probability based on paternal map 
 *
 *  Args:     ...
 *
 *  Return:   None
 *
 */
inline double calculate_paternal_condition_prob(HapElem ***hm, RecLOD **rm, Coordinate *ss, char hap, int pLen, gsl_rng *r){
	// get indi and loci info 
	int indi = ss->indi;
//	int pPos = ss->pMapPos;
	int loci = ss->loci;
	int pPreLoc = ss->pPreLoc;
	int pNextLoc = ss->pNextLoc;
	// sampling 
	char pPreHap, pNextHap;
	double pPreProb, pNextProb, p = 1;
	if (pPreLoc == -1 && pNextLoc != -1){ // head of p map
		pNextHap = hm[pNextLoc][indi]->hapP;
		pNextProb = rm[loci][pNextLoc].pRec;
		p = calculate_conditional_prob_ht(hap, pNextHap, pNextProb);
	}else if (pPreLoc != -1 && pNextLoc == -1){ // tail of p map 
		pPreHap = hm[pPreLoc][indi]->hapP;
		pPreProb = rm[loci][pPreLoc].pRec;
		p = calculate_conditional_prob_ht(hap, pPreHap, pPreProb);	
	}else if (pPreLoc != -1 && pNextLoc != -1){ // middle position of p map 
		pPreHap = hm[pPreLoc][indi]->hapP;
		pPreProb = rm[loci][pPreLoc].pRec;
		pNextHap = hm[pNextLoc][indi]->hapP;
		pNextProb = rm[loci][pNextLoc].pRec;
		p = calculate_conditional_prob_middle(hap, pPreHap, pNextHap, pPreProb, pNextProb);	
	}
//	else{
//		BioDie("can't find adjacent loci during calculation of paternal conditional probability");
//	}
//	if (p == -1) BioDie("Error in calculation of paternal conditional probability");
	return p;
}
/*
 *  Function: calculate_maternal_condition_prob()
 *
 *  Purpose:  calculate conditional probability based on maternal map 
 *
 *  Args:     ...
 *
 *  Return:   None
 *
 */

inline double calculate_maternal_condition_prob(HapElem ***hm, RecLOD **rm, Coordinate *ss, char hap, int mLen, gsl_rng *r){
	// get loci and indi info 
	int indi = ss->indi;
//	int mPos = ss->mMapPos;
	int loci = ss->loci;
	int mPreLoc = ss->mPreLoc;
	int mNextLoc = ss->mNextLoc;
	// sampling 
	char mPreHap, mNextHap;
	double mPreProb, mNextProb, p = 1;
	if (mPreLoc == -1 && mNextLoc != -1){ // head of m map
		mNextHap = hm[mNextLoc][indi]->hapM;
		mNextProb = rm[loci][mNextLoc].mRec;
		p = calculate_conditional_prob_ht(hap, mNextHap, mNextProb);
	}else if (mPreLoc != -1 && mNextLoc == -1){ // tail of m map 
		mPreHap = hm[mPreLoc][indi]->hapM;
		mPreProb = rm[loci][mPreLoc].mRec;
		p = calculate_conditional_prob_ht(hap, mPreHap, mPreProb);
	}else if (mPreLoc != -1 && mNextLoc != -1){ // middle position of m map 
		mPreHap = hm[mPreLoc][indi]->hapM;
		mPreProb = rm[loci][mPreLoc].mRec;
		mNextHap = hm[mNextLoc][indi]->hapM;
		mNextProb = rm[loci][mNextLoc].mRec;
		p = calculate_conditional_prob_middle(hap, mPreHap, mNextHap, mPreProb, mNextProb);
	}
//	else{
//		BioDie("can't find adjacent loci during calculation of maternal conditional probability");
//	}
//	if (p == -1) BioDie("Error in calculation of maternal conditional probability");
	return p;
}
/* reset all elem of recMatrixSum to 0 */
inline void reset_sum_recMatrix(RecLOD **rm, int *map, int maplen){
	int i, j;
	for (i=0; i<maplen; i++){
		for (j=0; j<maplen; j++){
			rm[i][j].pRec = 0;
			rm[j][i].pRec = 0;
			rm[i][j].mRec = 0;
			rm[j][i].mRec = 0;
		}
	}
}
/*
 *  Function: init_miss_observation()
 *
 *  Purpose:  imput missing data(random) to generate a complete haplo matrix
 *
 *  Args:     GibbsSAM *gSam
 *
 *  Return:   None
 *
 */

inline void init_miss_observation(GibbsSAM *gSam){
	// get sample space 
	Coordinate **ss = gSam->sampleSpace;
	int samLen = gSam->sampleLen;
	// get haplo matrix 
	MCEM *mcem = gSam->mcem;
	HapElem ***hm = mcem->hapMatrix;
	// random number generator 
	gsl_rng *r = gSam->r;

	// imputation of missing data randomly 
	int i, indi, loci;
	double p;
	HapElem *hap;
	for (i=0; i<samLen ; i++){
		if (ss[i]->isMulti == 0){ // miss data 
			p = gsl_rng_uniform(r);
			// locate miss data 
			indi = ss[i]->indi;
			loci = ss[i]->loci;
			hap = hm[loci][indi];
			if (p >= 0 && p < 0.25){
				hap->hapP = 0;
				hap->hapM = 0;
			}else if (p >= 0.25 && p < 0.5){
				hap->hapP = 0;
				hap->hapM = 1;
			}else if (p >= 0.5 && p < 0.75){
				hap->hapP = 1;
				hap->hapM = 0;
			}else{
				hap->hapP = 1;
				hap->hapM = 1;
			}	
		}else if (ss[i]->isMulti == 2){ // patial informative 
			p = gsl_rng_uniform(r);
			// locate miss data 
			indi = ss[i]->indi;
			loci = ss[i]->loci;
			hap = hm[loci][indi];
			if ( p >= 0.5){
				hap->hapP = 1 - hap->hapP;
				hap->hapM = 1 - hap->hapM;
			}	
		}else{
			BioDie("not missing observation and partial informative");
		}
	}
}
/*
 *  Function: init_map_and_sample_space()
 *
 *  Purpose:  init p and m maps based on the best map order obtained by SA, init gibbs sample space 
 *
 *  Args:     GibbsSAM *gSam  
 *
 *  Return:   None
 *
 */
inline void init_map_and_sample_space(GibbsSAM *gSam){
	// get gMatrix 
	GenoMatrix **gm = gSam->gMatrix;
	// get haoloMatrix 
	HapElem ***hm = gSam->hapMatrix;
	// get rec.freq matrix 
	RecLOD **rm = gSam->mcem->recMatrix;
	// get map and maplen
	int *map = gSam->sMap;
	int maplen = gSam->sLen;
	int *pMap = gSam->pMap;
	int *mMap = gSam->mMap;
	// pm map pos
	int **pmMapPos = gSam->pmMapPos;
	// get sample space 
	Coordinate **ss = gSam->sampleSpace;

	// iterate init sample space and map info 
	int i, j, a, mPos = 0, pPos = 0, sampleLen = 0; 
	int pPreLoc, pNextLoc, mPreLoc, mNextLoc;

	for (i=0; i<maplen ;i++){
		a = map[i];
		if (strcmp(gm[a]->type,"abxcd") == 0 
			|| strcmp(gm[a]->type,"efxeg") == 0
			|| strcmp(gm[a]->type,"hkxhk") == 0){
			// P and M map info 
			pMap[pPos] = a;   // male map
			mMap[mPos] = a; // female map
			pmMapPos[i][0] = pPos;  // map postion info: 0 for p and 1 for m
			pmMapPos[i][1] = mPos;  // map postion info: 0 for p and 1 for m
			
			// get prefix/next sex-specific loci in sMap
//			pPreLoc = get_pre_sexAver_loc(rm, gm, map, i, 0);
//			pNextLoc = get_next_sexAver_loc(gm, map, i, maplen, 0);
//			mPreLoc = get_pre_sexAver_loc(gm, map, i, 1);
//			mNextLoc = get_next_sexAver_loc(gm, map, i, maplen, 1);
		
			// init sample sapce 
			for (j=0;j<gSam->nind ;j++ ){
				if (hm[a][j]->isMulti == 0 || hm[a][j]->isMulti == 2){
					ss[sampleLen]->type = hm[a][j]->type; // loci type 
					ss[sampleLen]->isMulti = hm[a][j]->isMulti; // is multiple haplo 
					ss[sampleLen]->indi = j;  // indi id
					ss[sampleLen]->loci = a;  // loci id 
					
//					ss[sampleLen]->pPreLoc = pPreLoc;  // pre paternal loci in sexAver map 
//					ss[sampleLen]->pNextLoc = pNextLoc; // next paternal loci in sexAver map
//					ss[sampleLen]->mPreLoc = mPreLoc;  // pre maternal loci in sexAver map
//					ss[sampleLen]->mNextLoc = mNextLoc; // next maternal loci in sexAver map

					ss[sampleLen]->pPreLoc = get_pre_adjacent_loc(rm, gm, hm, map, i, j, gSam->nind, ss[sampleLen]->pBinLoc, &(ss[sampleLen]->pBinLen), 0);
					ss[sampleLen]->pNextLoc = get_next_adjacent_loc(rm, gm, hm, map, i, maplen, j, gSam->nind, ss[sampleLen]->pBinLoc, &(ss[sampleLen]->pBinLen), 0);
					ss[sampleLen]->mPreLoc = get_pre_adjacent_loc(rm, gm, hm, map, i, j, gSam->nind, ss[sampleLen]->mBinLoc, &(ss[sampleLen]->mBinLen), 1);
					ss[sampleLen]->mNextLoc = get_next_adjacent_loc(rm, gm, hm, map, i, maplen, j, gSam->nind, ss[sampleLen]->mBinLoc, &(ss[sampleLen]->mBinLen), 1);
					
					ss[sampleLen]->sMapPos = i;        // sexAver map position 
					ss[sampleLen]->pMapPos = pPos;  // sexAver 
					ss[sampleLen]->mMapPos = mPos;
					// sample space size + 1 
					sampleLen++;
				}
			}

			// male and female map pos 
			pPos++;
			mPos++;
		}else if (strcmp(gm[a]->type,"lmxll") == 0){
			pMap[pPos] = a;   // male map
			pmMapPos[i][0] = pPos;  // map postion info: 0 for p and 1 for m
			pmMapPos[i][1] = -1;  // means not exist in female map

			// get prefix/next sex-specific loci in sMap
//			pPreLoc = get_pre_sexAver_loc(gm, map, i, 0);
//			pNextLoc = get_next_sexAver_loc(gm, map, i, maplen, 0);
//			mPreLoc = -1;
//			mNextLoc = -1;

			// init sample sapce 
			for (j=0;j<gSam->nind ;j++ ){
				if (hm[a][j]->isMulti == 0 || hm[a][j]->isMulti == 2){
					ss[sampleLen]->type = hm[a][j]->type; // loci type 
					ss[sampleLen]->isMulti = hm[a][j]->isMulti; // is multiple haplo 
					ss[sampleLen]->indi = j;  // indi id
					ss[sampleLen]->loci = a;  // loci id 
					
//					ss[sampleLen]->pPreLoc = pPreLoc;  // pre paternal loci in sexAver map 
//					ss[sampleLen]->pNextLoc = pNextLoc; // next paternal loci in sexAver map
//					ss[sampleLen]->mPreLoc = mPreLoc;  // pre maternal loci in sexAver map
//					ss[sampleLen]->mNextLoc = mNextLoc; // next maternal loci in sexAver map
					
					ss[sampleLen]->pPreLoc = get_pre_adjacent_loc(rm, gm, hm, map, i, j, gSam->nind, ss[sampleLen]->pBinLoc, &(ss[sampleLen]->pBinLen), 0);
					ss[sampleLen]->pNextLoc = get_next_adjacent_loc(rm, gm, hm, map, i, maplen, j, gSam->nind, ss[sampleLen]->pBinLoc, &(ss[sampleLen]->pBinLen), 0);
					ss[sampleLen]->mPreLoc = -1;
					ss[sampleLen]->mNextLoc = -1;

					
					ss[sampleLen]->sMapPos = i;        // sexAver map position 
					ss[sampleLen]->pMapPos = pPos;  // sexAver 
					ss[sampleLen]->mMapPos = mPos;
					// sample space size + 1 
					sampleLen++;
				}
			}
			pPos++;
		}else if (strcmp(gm[a]->type,"nnxnp") == 0){
			mMap[mPos] = a;   // male map
			pmMapPos[i][0] = -1;  // means not exist in female map
			pmMapPos[i][1] = mPos;  // map postion info: 0 for p and 1 for m

			// get prefix/next sex-specific loci in sMap
//			mPreLoc = get_pre_sexAver_loc(gm, map, i, 1);
//			mNextLoc = get_next_sexAver_loc(gm, map, i, maplen, 1);
//			pPreLoc = -1;
//			pNextLoc = -1;

			// init sample sapce 
			for (j=0;j<gSam->nind ;j++ ){
				if (hm[a][j]->isMulti == 0 || hm[a][j]->isMulti == 2){
					ss[sampleLen]->type = hm[a][j]->type; // loci type 
					ss[sampleLen]->isMulti = hm[a][j]->isMulti; // is multiple haplo 
					ss[sampleLen]->indi = j;  // indi id
					ss[sampleLen]->loci = a;  // loci id 
					
//					ss[sampleLen]->pPreLoc = pPreLoc;  // pre paternal loci in sexAver map 
//					ss[sampleLen]->pNextLoc = pNextLoc; // next paternal loci in sexAver map
//					ss[sampleLen]->mPreLoc = mPreLoc;  // pre maternal loci in sexAver map
//					ss[sampleLen]->mNextLoc = mNextLoc; // next maternal loci in sexAver map
					
					ss[sampleLen]->pPreLoc = -1;
					ss[sampleLen]->pNextLoc = -1;
					ss[sampleLen]->mPreLoc = get_pre_adjacent_loc(rm, gm, hm, map, i, j, gSam->nind, ss[sampleLen]->mBinLoc, &(ss[sampleLen]->mBinLen), 1);
					ss[sampleLen]->mNextLoc = get_next_adjacent_loc(rm, gm, hm, map, i, maplen, j, gSam->nind, ss[sampleLen]->mBinLoc, &(ss[sampleLen]->mBinLen), 1);

					ss[sampleLen]->sMapPos = i;        // sexAver map position 
					ss[sampleLen]->pMapPos = pPos;  // sexAver 
					ss[sampleLen]->mMapPos = mPos;
					// sample space size + 1 
					sampleLen++;
				}
			}
			mPos++;
		}else{
			BioDie("Error: unknow type");
		}
	}

	// map len and sample space size 
	gSam->sampleLen = sampleLen;
	gSam->pLen = pPos;
	gSam->mLen = mPos;

	// debug: dump sample space 
//	int k, t, r;
//	printf("sam len = %d\n", gSam->sampleLen);
//	for (i=0; i<gSam->sampleLen; i++){
//		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", i, ss[i]->loci, ss[i]->indi, ss[i]->pBinLen, ss[i]->mBinLen, ss[i]->pPreLoc, ss[i]->pNextLoc, ss[i]->mPreLoc, ss[i]->mNextLoc);
//		t = ss[i]->loci;
//		for (j=0; j<ss[i]->pBinLen; j++){
//			r = ss[i]->pBinLoc[j];
//			printf("%d(%lf  %lf),", r, rm[t][r].pRec, rm[t][r].sRec);
//		}
//		printf("\t");
//		for(k=0; k<ss[i]->mBinLen; k++){
//			r = ss[i]->mBinLoc[k];
//			printf("%d(%lf  %lf),", r, rm[t][r].mRec, rm[t][r].sRec);		
//		}
//		printf("\n");
//	}	

}

/*
 *  Function: get_pre_adjacent_loc()
 *
 *  Purpose:  get prefix paternal or maternal loci in sexAver map 
 *
 *  Args:     GenoMatrix **gMatrix, int *map, int p, int flag
 *
 *  Return:   preloc - int 
 */
inline int get_pre_adjacent_loc(RecLOD **rm, GenoMatrix **gMatrix, HapElem ***hm, int *map, int p, int indi, int nind, int *bin, int *binLen, int flag){
	/* preloc = -1 means preloc not found, two cases involved:
	   1) head of map; 
	   2)all preloc is in the same bin with current loci, and observation in corresbonding individual is incomplete.
	*/
	int preLoc = -1; 
	*binLen = 0; // empty bin 
	if (p == 0) return preLoc;

	int i, a, loci = map[p];
	double episilon = 1.0 / nind;  // resolution of map: 1/nind, rec.feq of marker pairs that smaller than this threshold are viewed as bin  
	if (flag == 0){ // flag 0 for paternal 
		for (i=p-1; i>=0; i--){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"nnxnp") != 0){
				/* preloc selection*/
				if (hm[a][indi]->isMulti == 1){
					preLoc = a; 
					return preLoc;
				}else if (rm[a][loci].pRec > episilon){ // a and loci not in the same bin 
					preLoc = a;
					return preLoc;
				}else{ // a and loci in the same bin 
					bin[*binLen] = a;
					(*binLen)++;
				}
			}
		}
	}else if (flag == 1){ // flag 1 for maternal
		for (i=p-1; i>=0; i--){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"lmxll") != 0){
				/* preloc selection*/
				if (hm[a][indi]->isMulti == 1){
					preLoc = a; 
					return preLoc;
				}else if (rm[a][loci].mRec > episilon){ // a and loci not in the same bin 
					preLoc = a;
					return preLoc;
				}else{ // a and loci in the same bin 
					bin[*binLen] = a;
					(*binLen)++;
				}
			}
		}
	}else{
		BioDie("unknow flag during getting prefix pos");
	}

	return preLoc;
}
/*
 *  Function: get_next_adjacent_loc()
 *
 *  Purpose:  get next paternal or maternal loci in sexAver map 
 *
 *  Args:     GenoMatrix **gMatrix, int *map, int p, int flag
 *
 *  Return:   next - int 
 */
inline int get_next_adjacent_loc(RecLOD **rm, GenoMatrix **gMatrix, HapElem ***hm, int *map, int p, int maplen, int indi, int nind, int *bin, int *binLen, int flag){
	/* note: i am too tired to explain*/
	int next = -1;
	*binLen = 0;
	if (p == maplen - 1) return next;

	int i, a, loci = map[p];
	double episilon = 1.0 / nind; // resolution of map: 1/nind, rec.feq of marker pairs that smaller than this threshold are viewed as bin 
	if (flag == 0){
		for (i=p+1; i<maplen; i++){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"nnxnp") != 0){
				/* next selection*/
				if (hm[a][indi]->isMulti == 1){
					next = a; 
					return next;
				}else if (rm[a][loci].pRec > episilon){ // a and loci not in the same bin 
					next = a;
					return next;
				}else{ // a and loci in the same bin 
					bin[*binLen] = a;
					(*binLen)++;
				}
			}
		}
	}else if (flag == 1){
		for (i=p+1; i<maplen; i++){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"lmxll") != 0){
				/* next selection*/
				if (hm[a][indi]->isMulti == 1){
					next = a; 
					return next;
				}else if (rm[a][loci].mRec > episilon){ // a and loci not in the same bin 
					next = a;
					return next;
				}else{ // a and loci in the same bin 
					bin[*binLen] = a;
					(*binLen)++;
				}
			}
		}
	}else{
		BioDie("unknow flag during getting prefix pos");
	}
	
	return next;
}


/*
 *  Function: get_pre_sexAver_loc()
 *
 *  Purpose:  get prefix paternal or maternal loci in sexAver map 
 *
 *  Args:     GenoMatrix **gMatrix, int *map, int p, int flag
 *
 *  Return:   None
 */
inline int get_pre_sexAver_loc(GenoMatrix **gMatrix, int *map, int p, int flag){
	int preLoc = -1; // means the head of map
	if (p == 0) return preLoc;
	
	int i, a;
	if (flag == 0){
		for (i=p-1;i>=0 ;i-- ){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"nnxnp") != 0 && strcmp(gMatrix[a]->type,"hkxhk") != 0){
				preLoc = a;
				return preLoc;
			}
		}
	}else if (flag == 1){
		for (i=p-1;i>=0 ;i-- ){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"lmxll") != 0 && strcmp(gMatrix[a]->type,"hkxhk") != 0){
				preLoc = a;
				return preLoc;
			}
		}
	}else{
		BioDie("unknow flag during getting prefix pos");
	}

	return preLoc;
}
/*
 *  Function: get_next_sexAver_loc()
 *
 *  Purpose:  get next paternal or maternal loci in sexAver map 
 *
 *  Args:     GenoMatrix **gMatrix, int *map, int p, int flag
 *
 *  Return:   None
 */
inline int get_next_sexAver_loc(GenoMatrix **gMatrix, int *map, int p, int maplen, int flag){
	int next = -1; // means the tail of map 
	if (p == maplen - 1) return next;
	
	int i, a;
	if (flag == 0){
		for (i=p+1;i<maplen;i++ ){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"nnxnp") != 0 && strcmp(gMatrix[a]->type,"hkxhk") != 0){
				next = a;
				return next;
			}
		}
	}else if (flag == 1){
		for (i=p+1;i<maplen ;i++ ){
			a = map[i];
			if (strcmp(gMatrix[a]->type,"lmxll") != 0 && strcmp(gMatrix[a]->type,"hkxhk") != 0){
				next = a;
				return next;
			}
		}
	}else{
		BioDie("unknow flag during getting next loci");
	}
	return next;
}

//===========================
// map function
//===========================

/* Function: haldane(); kosambi();
 *
 * Purpose:  haldane map function
 *           kosambi map function
 *
 * Args:     rec - rec
 *
 * Return:   d
 *
 */
inline double haldane(const double rec) {
   if( !(rec<0.5) ) BioDie("rec >= 0.5");
   return ( -0.5 * log(1-2*rec) );
}
inline double kosambi(const double rec) {
   if( !(rec<0.5) ) BioDie("rec >= 0.5");
   double r2 = 2*rec;
   return ( 0.25 * log((1+r2)/(1-r2)) );
}

//===========================
//  inverse map function 
//===========================

/* Function:  inverseHaldane(); inverseKosambi();
 *
 * Purpose:   inverse haldane map function
 *            inverse kosambi map function
 *
 * Args:      d - d
 *
 * Return:    rec
 *
 */
inline double inverseHaldane(double d) {
   if( d <= -0.000001 ) BioDie("d < -0.000001");   // d must > negativeThreshold; 
   if( d <= 0.0 ) d = 0.000001;      // if(d<=0.0) log(inverseHaldane(d)) is -inf
   return ( 0.5 * (1-exp(-2*d)) );
}
inline double inverseKosambi(double d) {
   if( d <= -0.000001 ) BioDie("d < -0.000001");   // d must > negativeThreshold; 
   if( d <= 0.0 ) d = 0.000001;      // if(d<=0.0) log(inverseKosambi(d)) is -inf
   return ( 0.5 * tanh(2*d) );
}

//===========================
// calculate condition prob 
//===========================
/*
 *  Function: calculate_conditional_prob_middle()
 *
 *  Purpose:  calculate probability of individual haplotype in middle position of one parent map conditioned on 
 *            the map distance, adjacent markers' haplotypes 
 *
 *  Args:     haplo - individual haplotype of marker on position p
 *            preHaplo - individual haplotype of marker on position p-1 
 *            sufHaplo - individual haplotype of marker on position p+1 
 *            sufRec - recombination frquency between marker p and marker p+1  
 *            preRec - recombination frquency between marker p and marker p-1  
 *  Return:   p - conditional probability 
 *
 *  Details:  Suppose the individual observation of marker on position p is missing,
 *	          given that marker on position p-1 and p+1 are not recombinant, 
 *	          then either the markers on positions p每1 and p are not recombinant and also the markers on positions p and p+1 are not
 *	          recombinant, or the markers on positions p每1 and pare recombi-nant and also the markers on positions p and p+1 are recombinant.
 *	          Given that the markers on positions p每1 and p+1 are recombi-nant, then either the markers on positions p每1 and p are recombi-nant but the markers on positions p and p+1 are not recombinant,
 *	          or the markers on positions p每1 and p are not recombinant but the
 *	          markers on positions p and p+1 are recombinant.
 */
inline double calculate_conditional_prob_middle(char haplo, char preHaplo, char sufHaplo, double preRec, double sufRec){

	double p ; // conditional probability to be calculated 
	if (preHaplo ^ sufHaplo){ // if true, then p-1, p+1 are recombinant
		if (haplo ^ preHaplo){ // if true , then p-1, p are recombinant and p, p+1 are not recombinants
			p = preRec * (1 - sufRec) / (preRec * (1 - sufRec) + (1 - preRec) * sufRec);
		}else{
			p = (1 - preRec) * sufRec / (preRec * (1 - sufRec) + (1 - preRec) * sufRec);
		}
	}else{ // p-1, p+1 are not recombinant
		if (haplo ^ preHaplo){ // if true, then p-1,p+1 are recombinant and p, p+1 are recombinant
			p = preRec * sufRec / (preRec * sufRec + (1 - preRec) * (1 - sufRec));
		}else{
			p = (1 - preRec) * (1 - sufRec) / (preRec * sufRec + (1 - preRec) * (1 - sufRec));
		}
	}
	return p;
}
/*
 *  Function: calculate_conditional_prob_ht()
 *
 *  Purpose:  calculate conditional probability when loci is the head or the tail of map
 *
 *  Args:     char haplo, char adjHap, double adjRec
 *
 *  Return:   conditional probability 
 *
 */
inline double calculate_conditional_prob_ht(char haplo, char adjHap, double adjRec){
	double  p;
	p = (haplo ^ adjHap) ? adjRec : 1 - adjRec ;
	return p;
}
