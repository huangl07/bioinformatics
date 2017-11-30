/*
 *
 *  Spatial sampling 
 *
 */

// Copyright 2012,  Ma Chouxian <macx@biomarker.com.cn>

#include "spatsam.h"

//===================================
// static function
//===================================
static SpatSAM *init_spatial_samping(LocsINFO *locsINFO, TSPMSET *tspm);
static void _spatial_sampling(SpatSAM *ss); 
static void update_T_S(SpatSAM *ss, TSPMSET *tspm);
static void get_anchor_loci(GenoMatrix **gm, int *isAnchor, int *lociList, int len, int *anchor, int *anchorLen);
static void free_spatSAM(SpatSAM **ss);

//===================================
// interface of spatial sampling 
//===================================
int spatial_sampling(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	// init 
	SpatSAM *ss = init_spatial_samping(locsINFO, tspm);

	// spatial sampling process 
	_spatial_sampling(ss);
	
	// update T and S
	update_T_S(ss, tspm);

	// output log 
	fprintf(opt->logfp,"=========================\n");
	fprintf(opt->logfp, "spatial sampling of loci: %d\n", ss->round);
	fprintf(opt->logfp, "  rec.frq. threshold: %.3f\n", ss->radius);
	fprintf(opt->logfp, "  Nr. of loci in sample: %d/%d\n", tspm->Slen, locsINFO->nloc);
	if (ss->samLen == 0){
		fprintf(opt->logfp, "  no change, continuing with next threshold\n");	
		// free SpatSAM object 
		free_spatSAM(&ss);
		return 0;
	}else{

		fprintf(opt->logfp, "start with map order:\n");
		fprintf(opt->logfp, "%-5s %-5s %-20s %-10s\n","pos.","Nr.","locus","seg.type");
		fprintf(opt->logfp, "%-5s %-5s %-20s %-10s\n","----","---","----------------","--------");

		int i;
		for (i=0;i<tspm->Slen ;i++ ){
			fprintf(opt->logfp, "%3d %5d    %-20s %-15s\n",i, tspm->S[i],locsINFO->gMatrix[tspm->S[i]]->marker,locsINFO->gMatrix[tspm->S[i]]->type);
		}
		fprintf(opt->logfp,"\n");

		free_spatSAM(&ss);

		return 1;
	}
}
/*
 *  Function: init_spatial_samping()
 *
 *  Purpose:  initilization spatial sample struct 
 *
 *  Args:     LocsINFO *locsINFO, TSPMSET *tspm
 *
 *  Return:   SpatSAM *ss
 *
 */
static SpatSAM *init_spatial_samping(LocsINFO *locsINFO, TSPMSET *tspm){
	SpatSAM *ss = (SpatSAM *) BioMalloc(sizeof(SpatSAM) * 1);
	
	// round 
	ss->round = tspm->round;
	// sampling radius 
	ss->radius = tspm->radius;
	// sample space 
	ss->samSpace = tspm->T;
	// sample space len
	ss->ssl = tspm->Tlen;
	// genotype matrix 
	ss->gMatrix = locsINFO->gMatrix;
	// pairwise rec.freq matrix 
	ss->recLOD = locsINFO->recLOD[0];
	// indictive anchor loci 
	ss->isAnchor = locsINFO->isAnchor;

	// init random number generator 
	gsl_rng_env_setup();
	ss->T = gsl_rng_default;
	ss->r = gsl_rng_alloc(ss->T);
	gsl_rng_set(ss->r,(unsigned long int)time(NULL));

	// final spatial sample 
	ss->sample = (int *) BioMalloc(sizeof(int) * ss->ssl);
	ss->samLen = 0;
	ss->unSampled = (int *) BioMalloc(sizeof(int) * ss->ssl);
	ss->unSamLen = 0;

	return ss;
}

/*
 *  Function: update_T_S()
 *
 *  Purpose:  update T and S in TSPMSET
 *
 *  Args:     SpatSAM *ss -- spatial sampling struct 
 *            TSPMSET *tspm -- tspm 
 *  Return:   None
 *
 */
static void update_T_S(SpatSAM *ss, TSPMSET *tspm){
	if (ss->samLen == 0) return;
	// get origin sample len 
	int origin_Slen = tspm->Slen;

	// update S and T
	int i, j;
	for (i=0;i<ss->samLen ;i++ ){
		j = origin_Slen + i;
		tspm->S[j] = ss->sample[i];
		tspm->Slen++;
	}
	memcpy(tspm->T, ss->unSampled, sizeof(int) * ss->unSamLen);
	tspm->Tlen = ss->unSamLen;
}

/* Function:spatial_sampling()
 *
 * Purpose: perform spatial sampling 
 * 		
 *
 * Args:    int *samSpace   - sample space 
 *          int *nelem      - the number of element in sampling space  
 *          RecLOD **recLOD - recombination frequency matrix, use recLOD[i][j].sRec
 *          double r        - sampling radius
 *
 * Return:	pointer of spatial sample struct - spatSAM *
 *
 */

static void _spatial_sampling(SpatSAM *ss){
//	if (ss->ssl <= 0) BioDie("sample size is less than zero");
	// get gMatrix 
	GenoMatrix **gm = ss->gMatrix;
	// get recLOD 
	RecLOD **recLOD = ss->recLOD;
	// samspace 
	int *samSpace = ss->samSpace;
	int nelem = ss->ssl;
	// sampled
	int *sam = ss->sample;
	// unsam
	int *unSam = ss->unSampled;

	/* simple case: sampling space has only one element */
	if (nelem == 1){
		sam[0] = samSpace[0];
		ss->samLen = 1;
		return;
	}
	
	// define spatial sampling variables 
	int i, nbuf=0, p, samLoci;   
	int *S = (int *) BioMalloc(sizeof(int) * nelem);   // sample space in an iteration
	int *buf = (int *) BioMalloc(sizeof(int) * nelem); // buffer to record element reserved in an iteration 
	memcpy(S,samSpace,sizeof(int)*nelem);
	int Slen = nelem;
	// anchor loci 
	int *anchor = (int *) BioMalloc(sizeof(int) * nelem);
	int Alen = 0;

	// iterately sampling process
	while (Slen > 0){
		get_anchor_loci(gm, ss->isAnchor, S, Slen, anchor, &Alen);
		if (Alen == 0){
			p = gsl_rng_uniform_int(ss->r,Slen);
			samLoci = S[p];
			sam[ss->samLen++] = samLoci;
		}else{
			p = gsl_rng_uniform_int(ss->r,Alen);
			samLoci = anchor[p];
			sam[ss->samLen++] = samLoci;
		}
		// exclude elements which have shorter distance with selected element from S 
		for (i=0;i<Slen ;i++ ){
			if (S[i] == samLoci) continue; // S[i] is samLoci self
			if (recLOD[S[i]][samLoci].sRec < ss->radius ){ // S[i] is close to samLoci 
				unSam[ss->unSamLen++] = S[i];
				continue;
			}
			buf[nbuf++] = S[i];   // reserved elements are stored into buf
		}
		if (nbuf == 0) break; // all elements excluded 
		// update S 
		memcpy(S,buf,sizeof(int)*nbuf);
		Slen = nbuf;
		nbuf = 0;
	}

	// free 
	free(S);S=NULL;
	free(buf);buf=NULL;
	free(anchor);anchor=NULL;

}
/* free ss */
static void free_spatSAM(SpatSAM **ss){
		(*ss)->samSpace = NULL;
		(*ss)->isAnchor = NULL;
		(*ss)->gMatrix = NULL;
		(*ss)->recLOD = NULL;
		free((*ss)->sample);(*ss)->sample = NULL;
		free((*ss)->unSampled);(*ss)->unSampled = NULL;
		gsl_rng_free((*ss)->r);
		free(*ss); *ss =NULL;
		
}
/* get anchor loci in a list */
static void get_anchor_loci(GenoMatrix **gm, int *isAnchor, int *lociList, int len, int *anchor, int *anchorLen){
	// classification: informative or partial informative 
	*anchorLen = 0;
	int i, a, il = 0, pil = 0;
	int *informative = (int *) BioMalloc(sizeof(int) * len);
	int *partial_informative = (int *) BioMalloc(sizeof(int) * len);
	for (i=0; i<len; i++){
		a = lociList[i];
		if (isAnchor[a] == 1){
			if (strcmp(gm[a]->type, "abxcd") == 0 || strcmp(gm[a]->type, "efxeg") == 0){
				informative[il] = a;
				il++;
			}else if (strcmp(gm[a]->type, "hkxhk") == 0){
				partial_informative[pil] = a; 
				pil++;
			}else{
				BioDie("not anchor loci, %s, %s\n", gm[a]->marker, gm[a]->type);
			}
			anchor[*anchorLen] = a;
			(*anchorLen)++;
		}
	}
	// if exists informative loci, use informative loci as space , else use partial informative loci
	if (il != 0){
		memcpy(anchor, informative, sizeof(int) * il);
		*anchorLen = il;
	}else if (pil != 0){
		memcpy(anchor, partial_informative, sizeof(int) * pil);
		*anchorLen = pil;
	}else{
		*anchorLen = 0;
	}

	free(informative);informative = NULL; 
	free(partial_informative);partial_informative = NULL; 
}

