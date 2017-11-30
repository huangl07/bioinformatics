/*
 *
 * Simulating annealing library for map optimization 
 *
 */

// Copyright 2012, Ma Chouxian <macx@biomarker.com.cn>

#include "siman.h"

//===========================
// static methods 
//===========================

/* static functions */
static inline int _simulated_annealing(simAN *sa);
static inline simAN *init_siman(int *map,int maplen,LocsINFO *locsINFO, MapOpt *opt);
static inline void free_siman(simAN **sa);
static inline Map *new_map (int *map,int maplen,LocsINFO *locsINFO, int *isfix);
static inline void copy_map(Map *dest,Map *src);
static inline void free_map(Map *map);
static inline void dump_simAN(simAN *sa,FILE *fp);
static inline int check_identical_order(simAN *sa, int *map);
// take step, generate new neighbor 
static inline void take_step(Map *map, LocsINFO *locsINFO, int *isfix, gsl_rng *r);
//static inline double calculate_delta_sarf(int from, int to, LocsINFO *locsINFO);
static inline void get_permissible_position(int *map, int maplen, int *isfix, int *pos, int *poslen, int rp);
static inline void get_rLocs_space(int *mapOrder, int maplen, int *isfix, int *rLocs, int *rl);
// energy function involved 
static inline void calculate_sarf(Map *map,LocsINFO *locsINFO); 
static inline void Ef(Map *map);

//===========================
// interface 
//===========================

/* Function:map_optimization()
 *
 * Purpose: perform map optimization using simulating anealing 
 * 		
 *
 * Args:    int *map        - object map to be optimized, optimized result will be store into it
 *          int maplen      - map length   
 *          LocsINFO *locsINFO - locus information 
 *          MapOpt *opt        - calculation option of ML mapping 
 *
 * Return:  None
 *
 */

inline int map_optimization (LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	// get map sample 
	int *map = tspm->S;
	int maplen = tspm->Slen;
	
	if (maplen <= 2) return 0;
	
	// create an object of simAN and init it 
	simAN *sa = init_siman(map, maplen, locsINFO, opt);
	
	// simulating annealing process 
	int is_all_loci_fixed = _simulated_annealing(sa);
	if (is_all_loci_fixed == 0){
		fprintf(opt->logfp,"no optimization( order completely fixed )\n");
		free_siman(&sa);
		return 1;
	}
	
	// check wheather order changed 
	int isIdentical = check_identical_order(sa, map);

	if (isIdentical == 0){ // not changed
		fprintf(opt->logfp,"optimization resulted in identical order\n");
		// free sa
		free_siman(&sa);
		return isIdentical;
	}else{
		/* at the end, copy the result onto the initial map, 
	   so we pass it back to the caller */
		memcpy(map,sa->bestmap->mapOrder,sizeof(int)*maplen);

		// dump sa info 
		dump_simAN(sa,opt->logfp);

		fprintf(opt->logfp, "map order after simulated annealing:\n");
		fprintf(opt->logfp, "%-5s %-5s %-20s %-10s\n","pos.","Nr.","locus","seg.type");
		fprintf(opt->logfp, "%-5s %-5s %-20s %-10s\n","----","---","----------------","--------");
		int i;
		for (i = 0;i < maplen ; i++){
			fprintf(opt->logfp, "%3d %5d    %-20s %-15s\n",i, map[i],locsINFO->gMatrix[map[i]]->marker,locsINFO->gMatrix[map[i]]->type);
		}
		fprintf(opt->logfp,"\n");

		// free sa
		free_siman(&sa);
	}
	return 2; // means order changed 
}

/* Function:_simulated_annealing()
 *
 * Purpose: simulating anealing algorithm, the objective function is (sa->bestmap->p_sarf + sa->bestmap->m_sarf)
 * 		
 * Args:    simAN *sa  - simAN struct 
 *
 * Return:  None
 *
 */
static inline int _simulated_annealing(simAN *sa){
	if (sa == NULL) BioDie("sa is NULL");
	// get locsINFO 
	LocsINFO *locsINFO = sa->locsINFO;
	// get isfix 
	int *isfix = sa->isfixed;

	Map *curMap = sa->curMap;   // the map currently investigated in problem space 
	Map *newMap = sa->newMap;   // the newly generated map 
	Map *bestMap = sa->bestmap;  // the best map found so far
	if (curMap->rLocsLen == 0) return 0; // means all loci fixed 

	double T; // temperature 
	int unimprove_t = 0;
	/* chain 0 is used to determine initial temperature T0.
	 * seen: Constructing dense genetic linkage maps,J. Jansen ?A.G. de Jong ?J.W. van Ooijen,Theor Appl Genet (2001) 102:1113?122
	 */
	double l_,lSum=0, pi = -1, bolt_P;
	int i, l_t = 0 ;  // l_ = lSum / l_t (when l < 0), expectation of l 
	for (i=0;i<sa->nstep ;i++ ){
		copy_map(newMap, curMap);
		take_step(newMap, locsINFO, isfix, sa->r);

		if (newMap->E < bestMap->E) copy_map(bestMap,newMap); // global best 
	
		sa->n_evals++;          /* keep track of Ef() evaluations */
		/* now take the crucial step: see if the new point is accepted
		or not, as determined by the boltzmann probability */
		if (newMap->E < curMap->E){
			if (newMap->E < bestMap->E) copy_map(bestMap, newMap);
			/* yay!!! take a step */
			copy_map(curMap, newMap);
			sa->n_eless++;
		}else if (newMap->E == curMap->E){
			/* case of exchanging lmxll and nnxnp, but paternal and maternal map order still unchanged ,or 
			   paternal and maternal map order changed, but have the same sarf(minor case),
			   or finement of one parental map is at cost of same extent deterioration of the other(minor case).
			   These cases are all rejected by SA system.
			 */
			sa->n_rejects++;
			continue;
		}else{
			bolt_P = boltzmann(curMap->E, newMap->E, sa->t_initial);
			if (pi == -1) pi = bolt_P; // first acceptance probability
			if (gsl_rng_uniform(sa->r) < bolt_P){
				/* yay!!! take a step */
				copy_map(curMap, newMap);
				sa->n_accepts++;
				lSum += curMap->E - newMap->E;
			}else{
				sa->n_rejects++;
				lSum += curMap->E - newMap->E;
			}
			l_t++;
		}
	}
	l_ = lSum / l_t;
    //T = l_ / log(0.2);
	T = l_ / log(pi); /* set initial temperature */
	// reset all stat to zero
//	sa->n_accepts = 0;
//	sa->n_rejects = 0;
//	sa->n_eless = 0;
//	sa->n_evals = 0;
	// annealing procedure from Chain 1 
	while (1){
		for (i=0;i<sa->nstep ;i++ ){
			copy_map(newMap, curMap);
			take_step(newMap, locsINFO, isfix, sa->r);
			if (newMap->E < bestMap->E) copy_map(bestMap,newMap);
			sa->n_evals++;          /* keep track of Ef() evaluations */
			/* now take the crucial step: see if the new point is accepted
			 or not, as determined by the boltzmann probability */
			if (newMap->E < curMap->E){
				if (newMap->E < bestMap->E) copy_map(bestMap, newMap);	
				/* yay!!! take a step */
				copy_map(curMap, newMap);
				sa->n_eless++;
				unimprove_t = 0; // reset unimprove_t to zero if objective function improved 
			}else if (newMap->E == curMap->E){
				sa->n_rejects++;
				unimprove_t++;
			}else if (gsl_rng_uniform(sa->r) < boltzmann(curMap->E, newMap->E, T)){	
				/* yay!!! take a step */
				copy_map(curMap, newMap);
				sa->n_accepts++;
				unimprove_t++;
			}else{
				sa->n_rejects++;
				unimprove_t++;
			}
			if (unimprove_t >= sa->unimproved_chain_max) break;
		}
		sa->n_iters++;
		sa->t_final = T;
		// cooling 
		T *= sa->cooling_factor;
		if (unimprove_t >= sa->unimproved_chain_max ) break; // Stop after # chains without improvement
	}
	return 1;
}

/* Function:take_step()
 *
 * Purpose: Generate a new map by replacing random locus to random position 
 * 		
 * Args:    simAN *sa  - simAN struct 
 *
 * Return:  None
 *
 */
static inline void take_step(Map *map, LocsINFO *locsINFO, int *isfix, gsl_rng *r){
	// get map length 
	int maplen = map->maplen;
	// get rLocs space 
	int *rLocs = map->rLocs;
	int rLocsLen = map->rLocsLen;

	// new map order 
	int *newMap = (int *) BioMalloc(sizeof(int) * maplen);
	memcpy(newMap, map->mapOrder, sizeof(int)*maplen);
	int *pos = (int *) BioMalloc(sizeof(int) * maplen); // permissible position 

	// generate new map by random walking
	int rpLocus, rpPos, pl;
	while (1){
		rpLocus = gsl_rng_uniform_int(r, rLocsLen); // random locus   
		// get permissible postion of rLocus 
		get_permissible_position(newMap, maplen, isfix, pos, &pl, rLocs[rpLocus]);
		if (pl == 0){
			// rLocs 
			printf("r locs len : %d\n", rLocsLen);
			printf("map: %d  %d  %d\n", newMap[rLocs[rpLocus]-1], newMap[rLocs[rpLocus]], newMap[rLocs[rpLocus]+1]);
			printf("is fixed: %d %d %d\n",isfix[newMap[rLocs[rpLocus]-1]], isfix[newMap[rLocs[rpLocus]]], isfix[newMap[rLocs[rpLocus]+1]]);
			printf("xxxx%d\t%d\n",rLocs[rpLocus], pl);
		}
		rpPos = gsl_rng_uniform_int(r, pl);  // random position 
		// generate newMap
		repalce_random_locus(newMap, rLocs[rpLocus], pos[rpPos]);
		break;
	}
	memcpy(map->mapOrder, newMap, sizeof(int)*maplen);
	// rLocs and rLocsLen 
	get_rLocs_space(newMap, maplen, isfix, map->rLocs, &(map->rLocsLen));
	// p_sarf and m_sarf 
	calculate_sarf(map, locsINFO);
	Ef(map);
	// free 
	free(newMap); newMap = NULL;
	free(pos); pos = NULL;
}

static inline void get_permissible_position(int *map, int maplen, int *isfix, int *pos, int *poslen, int rp){
	int i, pl = 0;
	if (isfix == NULL){ // no fix order 
		for (i = 0;i < rp ; i++) pos[i] = i;
		for (i = rp; i < maplen - 1 ; i++) pos[i] = i+1;
		*poslen = maplen - 1;
		return;
	}
	if (isfix[map[rp]] == 0){ // loci is free 
		for (i = 0;i < rp ; i++) pos[i] = i;
		for (i = rp; i < maplen - 1 ; i++) pos[i] = i+1;
		*poslen = maplen - 1;
		return;
	}
	// loci is in fix order 
	// forward search
	for (i = rp+1; i < maplen ; i++){
		if (isfix[map[i]] == 0)
			pos[pl++] = i;
		else
			break;
	}
	// backward search 
	for (i = rp-1; i >= 0  ; i--){
		if (isfix[map[i]] == 0)
			pos[pl++] = i;
		else
			break;
	}
	*poslen = pl;
}

/* raplace random locus to random position  */
inline void repalce_random_locus(int *newMap, int rLocus, int rPos){
	int i, swap;
	if (rLocus > rPos){
		for (i=rLocus;i>rPos ;i--){ // swap adjacent loci
			swap = newMap[i];
			newMap[i] = newMap[i-1];
			newMap[i-1] = swap;
		}
	}else{
		for (i=rLocus;i<rPos ;i++ ){
			swap = newMap[i];
			newMap[i] = newMap[i+1];
			newMap[i+1] = swap;
		}
	}
}

/* Function:boltzmann()
 *
 * Purpose: calculate boltzmann distribution 
 * 		
 * Args:    double E    -  former energy 
 *          double new_E - new step energy
 *          double T    - current temperature 
 *
 * Return:  boltzmann probability 
 *
 */
inline double boltzmann(double E, double new_E, double T){
	double x = -(new_E - E) / T;	
	/* avoid underflow errors for large uphill steps */
	return (x < GSL_LOG_DBL_MIN) ? 0.0 : exp(x);
}

/* Function: Ef() - energy function */
static inline void Ef(Map *map){
	map->E = map->p_sarf + map->m_sarf;
}

/* Function:init_siman()
 *
 * Purpose: allocate memory and initialization for simAN struct 
 * 		
 * Args:    int *map        - object map to be optimized, optimized result will be store into it
 *          int maplen      - map length   
 *          LocsINFO *locsINFO - locus information 
 *          MapOpt *opt        - calculation option of ML mapping 
 *
 * Return:	sa - simAN *
 *
 */

static inline simAN *init_siman(int *map,int maplen,LocsINFO *locsINFO, MapOpt *opt){
	// init 
	simAN *sa = (simAN *) BioMalloc(sizeof(simAN));
	// locusINFO 
	sa->locsINFO = locsINFO;

	// init rLocs, rLocsLen and isfixed
	if (locsINFO->fixOrderSum != 0)
		sa->isfixed = locsINFO->isfix[locsINFO->fixOrderSum - 1];
	else
		sa->isfixed = NULL;

	// init bestmap to original map
	sa->curMap = new_map(map, maplen, locsINFO, sa->isfixed);
	sa->bestmap = new_map(map, maplen, locsINFO, sa->isfixed);
	// allocate memory for new map 
	sa->newMap = (Map *) BioMalloc (sizeof(Map) * 1);
	sa->newMap->mapOrder = (int *) BioMalloc (sizeof(int) * maplen);
	sa->newMap->rLocs = (int *) BioMalloc (sizeof(int) * maplen);

	// gsl random number generator 
	gsl_rng_env_setup();
	sa->T = gsl_rng_default;
	sa->r = gsl_rng_alloc(sa->T);
	gsl_rng_set(sa->r,(unsigned long int)time(NULL));
	
	// stat 
	sa->n_evals = 0;
	sa->n_iters = 0;
	sa->n_rejects = 0;
	sa->n_accepts = 0;
	sa->n_eless = 0;
	
	// simulated annealing parameters 
	sa->nstep = opt->SA_nstep;
	sa->t_initial = opt->SA_initial_acceptance;
	sa->cooling_factor = 1/(1 + opt->SA_cooling_control); // Note: cooling_factor is not cooling parameter
	sa->unimproved_chain_max = opt->SA_unimproved_chain_max;

	return sa;
}
/* destructor of simAN */ 
static inline void free_siman(simAN **sa){
	free((*sa)->curMap);(*sa)->curMap = NULL;
	free((*sa)->newMap);(*sa)->newMap = NULL;
	free((*sa)->bestmap);(*sa)->bestmap = NULL;
	(*sa)->isfixed = NULL;
	gsl_rng_free((*sa)->r);
	free(*sa);*sa =NULL;
}
/* dump simAN */
static inline void dump_simAN(simAN *sa,FILE *fp){
	// get nind 
//	int nind = sa->locsINFO->nind;

	fprintf(fp,"map optimization:\n");
	fprintf(fp,"  iterations: %d (final T: %.5f)\n", sa->n_iters, sa->t_final);
	fprintf(fp,"  stopped at chain number: %d\n", sa->n_evals);
	fprintf(fp,"  Nr. of downhills: %d\n", sa->n_eless);
	fprintf(fp,"  Nr. of rejects: %d\n", sa->n_rejects);
	fprintf(fp,"  Nr. of accepts: %d\n\n", sa->n_accepts);
	fprintf(fp,"  sum of adjacent rec.freq: %.3f (P: %.3f, M: %.3f)\n\n", sa->bestmap->p_sarf + sa->bestmap->m_sarf, sa->bestmap->p_sarf , sa->bestmap->m_sarf );
//	fprintf(fp,"  sum of adjacent rec.freq: %.3f (P: %.3f, M: %.3f)\n\n", (sa->bestmap->p_sarf + sa->bestmap->m_sarf)/(double) nind, sa->bestmap->p_sarf /(double) nind, sa->bestmap->m_sarf /(double) nind);

}
/* check whether order changed */
static inline int check_identical_order(simAN *sa, int *map){
	// get map order 
	int *bestOrder = sa->bestmap->mapOrder;
	// get map len 
	int maplen = sa->bestmap->maplen;

	// check whether order changed 
	int i, isChange = 0;
	for (i = 0;i < maplen ; i++){
		if (bestOrder[i] != map[i] && bestOrder[i] != map[maplen-i-1]){
			isChange = 1;
			break;
		}
	}

	return isChange;
}

/* Function:new_map()
 *
 * Purpose: constructor of Map struct 
 * 		
 * Args:    int *map        - map 
 *          int maplen      - map length   
 *          LocsINFO *locsINFO - locus information 
 *
 * Return:	newMap - Map *
 *
 */

static inline Map *new_map (int *map, int maplen, LocsINFO *locsINFO, int *isfix){
	Map *newMap = (Map *) BioMalloc(sizeof(Map) * 1);
	int *mapOrder = (int *) BioMalloc(sizeof(int) * maplen);
	int *rLocs = (int *) BioMalloc (sizeof(int) * maplen);
	memcpy(mapOrder,map,sizeof(int)*maplen);  // use map to init new generated Map
	newMap->maplen = maplen;
	int rLocsLen ;
	get_rLocs_space(mapOrder, maplen, isfix, rLocs, &rLocsLen);

	// init 
	newMap->mapOrder = mapOrder;
	newMap->rLocs = rLocs;
	newMap->rLocsLen = rLocsLen;

	// init sarf, newmap->p_sarf and newmap->m_sarf and 
	calculate_sarf(newMap, locsINFO);
	Ef(newMap);   // energy of new configure

	return newMap;
}

static inline void get_rLocs_space(int *mapOrder, int maplen, int *isfix, int *rLocs, int *rl){
	// get rLocs space 
	int i, fa, a, na, rLocsLen = 0; 
	if (isfix == NULL){ // no fix order 
		for (i = 0;i < maplen ; i++) rLocs[i] = i;
		*rl = maplen;
		return;
	}
	for (i = 0;i < maplen ; i++){
		a = mapOrder[i];
		if (isfix[a] == 0){ // loci not in fix order 
			rLocs[rLocsLen++] = i;
		}else{ // loci is in fix order
			if (i == 0){
				na = mapOrder[1];
				if (isfix[na] == 0) rLocs[rLocsLen++] = i; // adjacent loci is not fixed 
			}else if (i == maplen - 1){
				fa = mapOrder[maplen - 2];
				if (isfix[fa] == 0) rLocs[rLocsLen++] = i;
			}else{
				fa = mapOrder[i-1];
				na = mapOrder[i+1];
				if (isfix[fa] == 0 || isfix[na] == 0) rLocs[rLocsLen++] = i;
			}
		}
	}
	*rl = rLocsLen; 
}

/* Function: copy_map() */ 
static inline void copy_map(Map *dest,Map *src){
	dest->maplen = src->maplen;
	memcpy(dest->mapOrder,src->mapOrder,sizeof(int)*src->maplen);
	dest->rLocsLen = src->rLocsLen;
	memcpy(dest->rLocs,src->rLocs,sizeof(int)*src->rLocsLen);
	dest->p_sarf = src->p_sarf;
	dest->m_sarf = src->m_sarf;
	dest->E = src->E;
}

/* destructor of map */
static inline void free_map (Map *map){
	free(map->mapOrder);map->mapOrder = NULL;
	free(map->rLocs);map->rLocs = NULL;
	free(map);map = NULL;
}

/* Function:calculate_sarf()
 *
 * Purpose: calculate sarf of map struct, 
 * 		
 * Args:    Map *map        - Map struct 
 *          LocsINFO *locsINFO - locus information 
 *
 * Return:  None
 */
static inline void calculate_sarf(Map *map, LocsINFO *locsINFO){
	// get map order 
	int *sOrder = map->mapOrder;
	if (sOrder == NULL) BioDie("error : null map during calculation of sarf");
	int maplen = map->maplen;
	map->p_sarf = 0;
	map->m_sarf = 0;
	// get gMatrix
	GenoMatrix **gMatrix = locsINFO->gMatrix;

	// get recLOD
	RecLOD **recLOD = locsINFO->recLOD[0];

	// calculate sarf in O(n), here it is unneccessary to generate female and male maps 
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
				map->p_sarf += recLOD[loci][pLast].pRec;
				pLast = loci;
			}
		}
		if (isMstart == 1  && mLast != loci){
			if (strcmp(gMatrix[loci]->type,"lmxll") != 0){
				map->m_sarf += recLOD[loci][mLast].mRec;
				mLast = loci;
			}
		}
	}
//	map->p_sarf *= locsINFO->nind;
//	map->m_sarf *= locsINFO->nind;
}

//#define SIMAN_TEST
#ifdef SIMAN_TEST

#include "mapio.h"
#include "mapopt.h"
int main (int argc, char *argv[]){
	// init spatial sample 

	MapOpt *opt = get_map_opt(argc, argv);
	write_map_opt_into_logfile(opt);
	LocsINFO *locsINFO = read_locs_info(opt);

	fprintf(opt->logfp, "\nMarkerID\tNr.\tMarkerID\tNr.\tMarkerID\tNr.\n");
	fprintf(opt->logfp, "----------\t-----\t----------\t-----\t----------\t-----\n");
	locsINFO->locs->showrecords(locsINFO->locs, opt->logfp);
	
	TSPMSET *tspm = (TSPMSET *) BioMalloc(sizeof(TSPMSET) * 1);
	tspm->T = (int *) BioMalloc(sizeof(int) * locsINFO->nloc); 
	tspm->S = (int *) BioMalloc(sizeof(int) * locsINFO->nloc); 
	tspm->Tlen = locsINFO->nloc;
	tspm->Slen = 0;
	tspm->round = 1;
	tspm->radius = 0;
	
	int i;
	for (i=0;i<locsINFO->nloc ;i++ ){
		tspm->T[i] = i;
	}
	
	int flag;
	flag = spatial_sampling(locsINFO, tspm, opt);

	double cpubegin = clock();
	long runbegin = time(NULL);
	map_optimization(locsINFO, tspm, opt);

	long runtime = time(NULL) - runbegin;
	printf("Done.Total elapsed time :%1d seconds or %.2lf mins or %.2lf hours\n",runtime, runtime/60.0, runtime/3600.0);
	
	return 0;
}
#endif
