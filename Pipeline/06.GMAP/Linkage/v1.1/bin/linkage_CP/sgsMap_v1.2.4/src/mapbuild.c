/*
 *
 *  mapbuild.c: construct linkage maps using maximum likelihood algorithm
 *
 *
 */

// Copyright 2012, Ma Chouxian <macx@biomarker.com.cn>

#include "mapbuild.h"

static void output_haplotype_matrix(LocsINFO *locsINFO, int *map, int maplen, FILE *fp);

//===================================
// utility
//===================================

/*
 *  Function: parse_command_line()
 *
 *  Purpose:  parse command line 
 *
 *  Args:    int  argc, char *argv[]
 *
 *  Return:   MapOpt *opt 
 *
 */
MapOpt *parse_command_line(int argc, char *argv[]){
	// get Map option 
	MapOpt *opt = get_map_opt(argc, argv);

	// write command line parameters into log file 
	write_map_opt_into_logfile(opt);

	return opt;
}

/*
 *  Function: read_loc_pwd_info()
 *
 *  Purpose:  read input files
 *
 */

LocsINFO *read_loc_pwd_info(MapOpt *opt){
	// get loc and pwd info 
	LocsINFO *locsINFO = read_locs_info(opt);

	// write marker index info to log file 
	fprintf(opt->logfp, "\nMarkerID\tNr.\tMarkerID\tNr.\tMarkerID\tNr.\n");
	fprintf(opt->logfp, "----------\t-----\t----------\t-----\t----------\t-----\n");
	locsINFO->locs->showrecords(locsINFO->locs, opt->logfp);

	return locsINFO;
}

/*
 *  Function: map_building()
 *
 *  Purpose: map building process 
 *
 */
void map_building(LocsINFO *locsINFO, MapOpt *opt){
	// init tspmset 
	TSPMSET *tspm = init_tspm(locsINFO, opt);

	// init tspm->T and tspm->S with start and fix order and build submap of start sample 
	map_building_start_sample(locsINFO, tspm, opt);

	// sequentially build map by spatial sampling, simulated annealing and  gibbs sampling 
	map_building_spatial_sample(locsINFO, tspm, opt);

	// output map file, pwd file, singleton statistic 
	output_result(locsINFO, tspm, opt);

	// pp matrix 
	generate_plausible_position_matrix(locsINFO, tspm, opt);
}

/*
 *  Function: map_building_spatial_sample()
 *
 *  Purpose:  map sort of a spatial sample
 *
 *  Args:     locsINFO, tspm, opt   
 *
 *  Return:   None
 *
 */

void map_building_spatial_sample(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	
	int sample_stage = 0, isempty;
	// sequentially building map 
	for (sample_stage = 0;sample_stage < 6; sample_stage++ ){
		// if all loci mapped but not the first sample, then quit 
		if (tspm->Tlen == 0 && tspm->isStartSample == 0) break;
		// spatial sampling 
		tspm->radius = opt->SS_rec_thresholds[sample_stage]; // set sampling radius 
		tspm->round = sample_stage + 1;
		isempty = spatial_sampling(locsINFO, tspm, opt);
		if (isempty == 0 && tspm->isStartSample == 0) continue;  // empty sample at this threshold of rec.freq, but not the first sample

		// simulated annealing + gibbs sampling (multi-cycle) to get best map of this sample 
		map_sort_one_sample(locsINFO, tspm, opt);

		if (tspm->isStartSample == 1) tspm->isStartSample = 0; // important 
	}
	// cancel fix order of last sample for generating plausible position using metropolis-hasting 
	fix_next_sample_first_round(locsINFO, tspm, 1);  // Note: important
}

/*
 *  Function: map_building_start_sample()
 *
 *  Purpose:  map sort of start sample, always initialized by start and fix order
 *
 */
void map_building_start_sample(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	/* All cases of start order and fix order:
		(a). only start order 
		(b). only fix order 
		(c). start order and fix order co-exist 
			(1). start order and fix order have no overlap 
			(2). start order and fix order have overlap 
	*/
	// case (a) case (c)
	if (locsINFO->startOrder != NULL){ // exists start order 
		if (locsINFO->fixOrderSum == 0){ // case (a)
			update_T_with_order(locsINFO->startOrder, locsINFO->startOrderLen, tspm);
			// map sort of start sample 
			map_sort_one_sample(locsINFO, tspm, opt);
			tspm->isStartSample = 0;
		}else{ // case (c)
			// merge info of start and fix order 
			merge_start_fix_order(locsINFO,tspm);
			// map sort of merged order 
			map_sort_one_sample(locsINFO, tspm, opt);
			tspm->isStartSample = 0;
		}
	}else{ // case (b): incorperate fix order to first spatial sample, here it is only need to update tspm->T
		if (locsINFO->fixOrderSum != 0) 
			update_T_with_order(locsINFO->fixOrder[0], locsINFO->fixOrderLen[0], tspm);		
	}
	
}
/*
 *  Function: update_T_with_order()
 *
 *  Purpose:  update T and S using map order, called when specified start and fix order 
 *
 */

void update_T_with_order(int *order, int orderLen, TSPMSET *tspm){
	// get T 
	int *T = tspm->T;
	int Tlen = tspm->Tlen;
	// get S
	int *S = tspm->S;
	
	int *tmp = (int *) BioMalloc(sizeof(int) * Tlen);
	// update T 
	int i, j, isExists = 0, count = 0;
	for (i=0;i<Tlen ;i++ ){
		isExists = 0;
		for (j=0;j<orderLen ;j++) 
			if (T[i] == order[j]) {isExists = 1; break;}

		if (isExists == 0)  tmp[count++] = T[i];
	}
	memcpy(T, tmp, sizeof(int) * count); // update T
	tspm->Tlen = count;

	memcpy(S, order, sizeof(int) * orderLen); // update S
	tspm->Slen = orderLen;
	
	// over and free 
	free(tmp);tmp = NULL;
}
/*
 *  Function: map_sort_one_sample()
 *
 *  Purpose:  map sort of one sample, muti-cycle simulated annealing plus gibbs sampling  
 *
 */

void map_sort_one_sample (LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	// map building round 
	int round = 0, flag; 
	tspm->cycle = 0;
	while (round < opt->MB_nopt){
		fprintf(opt->logfp, "---------------------------------\n");
		fprintf(opt->logfp, "rec. freq. threshold: %.3f, loci: %d/%d, opt.round: %d/%d\n\n", tspm->radius, tspm->Slen, locsINFO->nloc, round+1, opt->MB_nopt);
//		printf(">rec. freq. threshold: %.3f, loci: %d/%d, opt.round: %d/%d\n\n", tspm->radius, tspm->Slen, locsINFO->nloc, round+1, opt->MB_nopt);
		// map optimization 
		flag = map_optimization(locsINFO, tspm, opt);
		// after round 1, cancel fix order of preceding sample
		// then following rounds have no restriction of fixing preceding samples' order 
		if (tspm->isStartSample == 0 && round == 0) fix_next_sample_first_round(locsINFO, tspm, 1);

		if (flag == 0){
			fprintf(opt->logfp,"multi-point estimation of rec.frequencies:\n");
			fprintf(opt->logfp,"  identical order, rec.freq. not re-estimated\n\n");
			if (round != 0) {round++; continue;} // first sample also need multi-estimation of rec.freq
		}else if (flag == 1 && round != 0){ // order is totally fixed and not not first round  
			round++;
			continue;
		}
		// multi_estimation of reconmbination frequency
		multi_estimation_rf(locsINFO, tspm, opt);
		round++;
		tspm->cycle++;
	}

	// after map sort of this sample, the order of this sample is fixed in the first optimization round of next sample 
	fix_next_sample_first_round(locsINFO, tspm, 0);  // Note: this order is always appended in locsINFO->fixOrder  
	// over 
}

/*
 *  Function: fix_next_sample_first_round()
 *
 *  Purpose:  fix preceding samples' order or cancel it, 0 for fix and 1 for cancel 
 *
 *  Args:     
 *
 *  Return:   
 *
 */

void fix_next_sample_first_round(LocsINFO *locsINFO, TSPMSET *tspm, int flag){
	if (tspm->Slen < 3) return;
	// Note: 0 for adding fix order, 1 for cancel last fix order -- fix order of preceding samples  
	if (flag == 0){
		// get cur fix order sum 
		int n = locsINFO->fixOrderSum;
		locsINFO->fixOrderLen[n] = tspm->Slen;
		locsINFO->fixOrder[n] = (int *) BioMalloc(sizeof(int) * locsINFO->nloc);
		memcpy(locsINFO->fixOrder[n], tspm->S, sizeof(int) * tspm->Slen);
		memset(locsINFO->isfix[n], 0, sizeof(int) * locsINFO->nloc);
		int i;
		for (i = 0;i < tspm->Slen ; i++) locsINFO->isfix[n][tspm->S[i]] = 1;
		
		locsINFO->fixOrderSum += 1;
	}else{
		if (locsINFO->fixOrderSum != 0) locsINFO->fixOrderSum -= 1;
	}
	
}
/*
 *  Function: merge_start_fix_order();
 *
 *  Purpose:  merge start and fix order to a start sample
 *
 *  Args:     
 *
 *  Return:   
 *
 */
void merge_start_fix_order(LocsINFO *locsINFO, TSPMSET *tspm){
	// get start order info 
	int *startOrder = locsINFO->startOrder;
	int sLen = locsINFO->startOrderLen;

	// get fix order info 
	int *fixOrder = locsINFO->fixOrder[0];
	int fLen = locsINFO->fixOrderLen[0];

	// incorperate start order and fix order 
	int *mergeOrder = (int *) BioMalloc(sizeof(int) * locsINFO->nloc);
	int i, j, isExists = 0, mergeLen = 0;	
	memcpy(mergeOrder, fixOrder, sizeof(int) * fLen);
	mergeLen = fLen;
	for (i=0;i<sLen ;i++ ){
		for (j=0;j<fLen ;j++ )
			if (startOrder[i] == fixOrder[j]){ isExists = 1;break;}
		// not found in fix order 
		if (isExists == 0) mergeOrder[mergeLen++] = startOrder[i]; 
	}

	// update 
	update_T_with_order(mergeOrder, mergeLen, tspm);
	
	// free 
	free(mergeOrder); mergeOrder = NULL;
}
/*
 *  Function: int_tspm()
 *
 *  Purpose:  initialize struct TSPMSET 
 *
 *  Args:     
 *
 *  Return:   
 *
 */
TSPMSET *init_tspm(LocsINFO *locsINFO, MapOpt *opt){
	// get the number of loci
	int nloc = locsINFO->nloc;
	
	// new tspm
	TSPMSET *tspm = (TSPMSET *) BioMalloc(sizeof(TSPMSET) * 1);
	tspm->T = (int *) BioMalloc(sizeof(int) * nloc);
	tspm->S = (int *) BioMalloc(sizeof(int) * nloc);
	tspm->SD = (double *) BioMalloc(sizeof(double) * nloc); 
	tspm->P = (int *) BioMalloc(sizeof(int) * nloc);
	tspm->PD = (double *) BioMalloc(sizeof(double) * nloc); 
	tspm->M = (int *) BioMalloc(sizeof(int) * nloc);
	tspm->MD = (double *) BioMalloc(sizeof(double) * nloc); 
	// T, S, P, M 
	int i;
	for (i=0;i<nloc ;i++ ){
		tspm->T[i] = i; 
		tspm->S[i] = i; tspm->SD[i] = -1.0;
		tspm->P[i] = i; tspm->PD[i] = -1.0;
		tspm->M[i] = i; tspm->MD[i] = -1.0;
	}
	// Tlen, Slen, Plen, Mlen 
	tspm->Tlen = nloc;
	tspm->Slen = 0;
	tspm->Plen = 0;
	tspm->Mlen = 0;
	// sarf and singleton
	tspm->S_sarf = 0.0; tspm->S_singleton = 0.0;
	tspm->P_sarf = 0.0; tspm->P_singleton = 0.0;
	tspm->M_sarf = 0.0; tspm->M_singleton = 0.0;
	// round
	tspm->round = 0;
	tspm->radius = 0;
	tspm->cycle = 0;
	tspm->isStartSample = 1;

	return tspm;
}

//===================================
// output 
//===================================

void output_result(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt){
	// output map 
	FILE *sfp, *mfp, *pfp, *fp;
	if( (sfp = fopen(opt->sexAverMapFile, "w")) == NULL )
		BioDie("open or create sexAver map file %s is failed!", sfp);
	if( (mfp = fopen(opt->femaleMapFile, "w")) == NULL )
		BioDie("open or create male map file %s is failed!", mfp);
	if( (pfp = fopen(opt->maleMapFile, "w")) == NULL )
		BioDie("open or create female map file %s is failed!", pfp);

	output_map_file(locsINFO, tspm->S, tspm->SD, tspm->Slen, sfp);
	output_map_file(locsINFO, tspm->P, tspm->PD, tspm->Plen, pfp);
	output_map_file(locsINFO, tspm->M, tspm->MD, tspm->Mlen, mfp);

	// output pwd 
//	output_pwd_files();

	// output singleton statistic
	if ((fp = fopen(opt->singletonFile, "w")) == NULL){
		BioDie("open or create singleton statistic file %s is failed!", mfp);
	}
	fprintf(fp, "%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\n", opt->fKey, locsINFO->nind, tspm->Slen, tspm->Plen, tspm->Mlen, tspm->S_singleton, tspm->P_singleton, tspm->M_singleton);
	fclose(fp);
	
	// output haplomatrix 
	if ((fp = fopen(opt->haplotypeFile, "w")) == NULL){
		BioDie("open or create haplotype file %s is failed!", mfp);
	}
	output_haplotype_matrix(locsINFO, tspm->S, tspm->Slen, fp);
	fclose(fp);
}

void output_map_file(LocsINFO *locsINFO, int *map,double *mapdist, int maplen, FILE *fp){
	// get gMatrix 
	GenoMatrix **gm = locsINFO->gMatrix;
	
	fprintf(fp,"group 1\n\n");
	int i;
	for (i=0; i<maplen; i++){
		fprintf(fp, "%s\t%.3f\n", gm[map[i]]->marker, 100*mapdist[i]);
	}
}
static void output_haplotype_matrix(LocsINFO *locsINFO, int *map, int maplen, FILE *fp){
	// get nind 
	int nind = locsINFO->nind;
	// get genotype matrix 
	GenoMatrix **gm = locsINFO->gMatrix;
	// haplomatrix 
	HapElem ***hm = locsINFO->hapMatrix;
	
	int i, j, a;
	for (i=0; i<maplen; i++){
		a = map[i];
		fprintf(fp, "%s\t<%s>\t{%s}", gm[a]->marker, gm[a]->type, gm[a]->phase);
		if (strcmp(gm[a]->type, "abxcd") == 0 || strcmp(gm[a]->type, "efxeg") == 0 || strcmp(gm[a]->type, "hkxhk") == 0){
			for (j=0; j<nind; j++){
				fprintf(fp, "\t%d%d", hm[a][j]->hapP, hm[a][j]->hapM);
			}
		}else if (strcmp(gm[a]->type, "lmxll") == 0){
			for (j=0; j<nind; j++){
				fprintf(fp, "\t%d-", hm[a][j]->hapP);
			}
		}else if (strcmp(gm[a]->type, "nnxnp") == 0){
			for (j=0; j<nind; j++){
				fprintf(fp, "\t-%d", hm[a][j]->hapM);
			}
		}else{
			BioDie("unknown segregation type");
		}
		
		fprintf(fp, "\n");
	}
}


//#define MAPBUILD_TEST
#ifdef MAPBUILD_TEST

int main(int argc, char *argv[]) {
	// stat the time of map sort
	double cpubegin = clock();
	long runbegin = time(NULL);
	
	// parse command line option 
	MapOpt *opt = parse_command_line(argc, argv);
	
	// read loc and pwd info 
	LocsINFO *locsINFO = read_loc_pwd_info(opt);

	// get linear arrangement of loci in a linkage group 
	map_building(locsINFO, opt);

	// over and calculate elapsed time 
	double cputime = (clock() - cpubegin) / CLOCKS_PER_SEC;
	long runtime = time(NULL) - runbegin;
	fprintf(opt->logfp, "CPU time using: %.0lf seconds or %.3lf mins or %.3lf hours\n", cputime, cputime/60, cputime/3600);
	fprintf(opt->logfp, "Run time using: %ld seconds or %.3lf mins or %.3lf hours\n", runtime, runtime/60.0, runtime/3600.0);
	fclose(opt->logfp); // close log file 

	printf("Done.Total elapsed time :%1d seconds or %.2lf mins or %.2lf hours\n",runtime, runtime/60.0, runtime/3600.0);
	return 0;
}

#endif
