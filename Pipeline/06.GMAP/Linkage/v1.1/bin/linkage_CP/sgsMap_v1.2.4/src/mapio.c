/* 
 *	mapio.c: IO library of mlMap 
 * 
 */

// CopyRight 2012, Ma Chouxian <macx@biomarker.com.cn>

#include "mapio.h"

//===========================
//  static methods 
//===========================

static void read_linkage_group_locs(const char *filename, LocsINFO *locsINFO);
static void allocate_memory_ghMatrix(LocsINFO *locsINFO, int nloc, int nind);
static void init_locus(char **str, LocsINFO *locsINFO, int nind, int index);
static RecLOD ***allocate_memory_for_RecLOD(size_t n);
static RecLOD ***read_rec_LOD(const char *filename, LocsINFO *locsINFO);
static void init_start_order(LocsINFO *locsINFO, MapOpt *opt);
static void read_start_order_from_file(LocsINFO *locsINFO, MapOpt *opt);
static void read_fix_order_from_file(LocsINFO *locsINFO, MapOpt *opt);
static void init_fix_order(LocsINFO *locsINFO, MapOpt *opt);
static void check_start_fix_conflict(LocsINFO *locsINFO);
static inline void calculate_one_haplo(const char *type,const char *phase, const char *g, HapElem *h); // will be used in gibbs sampling 

//====================================
// interface of input 
//====================================

/*
 *  Function: read_locs_info()
 *
 *  Purpose:  interface function of input 
 *
 *  Args:     MapOpt *opt -- map option 
 *
 *  Return:   LocsINFO *locsINFO
 */
LocsINFO *read_locs_info(MapOpt *opt) {
	// init LocsINFO
	LocsINFO *locsINFO;
	locsINFO = (LocsINFO *) BioMalloc(sizeof(LocsINFO) * 1);

	// locs information
	// init BioStrIntHASH
	locsINFO->locs = (BioStrIntHASH *) BioMalloc(sizeof(BioStrIntHASH));
	newBioStrIntHASH(locsINFO->locs);
	read_linkage_group_locs(opt->locFile, locsINFO);
	
	// rec and LOD information
	// init recLOD
	locsINFO->recLOD = read_rec_LOD(opt->pwdFile, locsINFO);

	// init fix order
	init_fix_order(locsINFO, opt);

	// init start order
	init_start_order(locsINFO, opt);

	// check conflict of start order and fix order 
	check_start_fix_conflict(locsINFO);

	// over
	return locsINFO;
}

//===============================
//  Read loc file
//===============================
/*
 *  Function: read_linkage_group_locs()
 *
 *  Purpose:  read loc file 
 *
 *  Args:     filename, locsINFO 
 *
 *  Return:   None
 *
 */

static void read_linkage_group_locs(const char *filename, LocsINFO *locsINFO){
	// get locs
	BioStrIntHASH *locs = locsINFO->locs;

	// get filenames and open file handle
	FILE *fp;
	if( (fp = fopen(filename, "r")) == NULL ) BioDie("open loc file %s is failed!", filename);

	// get locs in the linkage group
	int nloc = 0, tmp, locsCount=0, nind = 0, lmxllCount = 0, nnxnpCount=0;
	char **str, *line, popt[3];
	regmatch_t *reg;

	// get loc head
	while ( (line=bio_readline(fp)) != NULL ){
		bio_chomp(line);
		// check blank line
		if( bio_is_blank_line(line) || bio_regex_m(line, "^;") || bio_regex_m(line, "^name") ) {
			free(line); line=NULL;
			continue;
		}	
		// skip population type
		if(bio_regex_m(line, "^popt")){
			sscanf(line, "popt = %s", popt);
			if( strcmp(popt,"CP") != 0 ) BioDie("population type error, %s", popt);
			free(line); line=NULL;
			continue;
		}
		// get the sum number of locs
		if(bio_regex_m(line, "^nloc")){
			sscanf(line, "nloc = %d", &nloc);
			if( nloc<3 ) BioDie("nloc(%d) < 3", nloc);
			free(line); line=NULL;
			if (nind != 0) break;
			else continue;
		}
		// get the sum number of individual
		if(bio_regex_m(line, "^nind")){
			sscanf(line, "nind = %d", &nind);
			free(line); line=NULL;
			if (nloc != 0) break;
			else continue;
		}
		// check 
		if (nloc == 0 || nind == 0) 
			BioDie("read loc error, nloc = %d, nind = %d",nloc, nind);
		else
			break;  // already get head info, quit temporarily
	}
	// allocate memory for locsINFO->gMatrix
	locsINFO->nloc = nloc;
	locsINFO->nind = nind;
	allocate_memory_ghMatrix(locsINFO, nloc, nind);
	locsINFO->isAnchor = (int *) BioMalloc(sizeof(int) * nloc);
	memset(locsINFO->isAnchor, 0, sizeof(int) * nloc);
	
	// continue reading loc 
	while ( (line=bio_readline(fp)) != NULL ){
		if (locsCount > nloc) BioDie("read loc error: locscount != nloc");
		bio_chomp(line);
		// check blank line and comment line 
		if( bio_is_blank_line(line) || bio_regex_m(line, "^;")) {
			free(line); line=NULL;
			continue;
		}
		// get locus info : MarkerID,type,phase,progeny genotypes 
		tmp = bio_regex_gms(line, "[^[:space:]]+", &reg, 1, &str);
		if( tmp < nind + 3) BioDie("read locs id and type error!loci: %s",str[0]);
		// add them into hash
		if(locs->findrecord(locs,str[0])!=NULL) BioDie("locs id repeat: %s", str[0]);
		locs->addrecord(locs, str[0], locsCount);
		if(bio_regex_m(line, "lmxll")){
			lmxllCount++;
		}else if(bio_regex_m(line, "nnxnp")){
			nnxnpCount++;
		}else{ // is a anchor loci 
			locsINFO->isAnchor[locsCount] = 1;
		}
		init_locus(str,locsINFO,nind, locsCount);
		locsCount++;

		// free
		bio_free_2D_array((void **)str, tmp); str = NULL;
		
		if(reg==NULL){free(reg); reg = NULL;}
	}

	// check the number of locs
	if( locs->recordsnum(locs) != nloc ) BioDie("read locs error: locs->recordsnum(locs) != nloc");

	// over
	fclose(fp);
	locsINFO->lmxllSum = lmxllCount;
	locsINFO->nnxnpSum = nnxnpCount;
}

/*
 *  Function: allocate_memory_ghMatrix()
 *
 *  Purpose:  allocate memory for GenoMatrix **gMatrix and hapMatrix 
 *
 *  Args:     locsINFO->nloc, locsINFO->nind
 *
 *  Return:   GenoMatrix **
 *
 */

static void allocate_memory_ghMatrix(LocsINFO *locsINFO, int nloc, int nind){

	GenoMatrix **gm = (GenoMatrix **) BioMalloc(sizeof(GenoMatrix *) * nloc);
	HapElem ***hapMatrix = (HapElem ***) BioMalloc(sizeof(HapElem **) * nloc);
	
	int i,j;
	for (i=0;i<nloc ;i++ ){
		gm[i] = (GenoMatrix *) BioMalloc(sizeof(GenoMatrix));
		// NOTE : Marker ID will be allacated memory dynamically  
		// type
		gm[i]->type = (char *) BioMalloc(sizeof(char) * 6);
		// phase 
		gm[i]->phase = (char *) BioMalloc(sizeof(char) * 3);
		// genotypes and haplotypes
		gm[i]->genotype = (char **) BioMalloc(sizeof(char *) * nind);
		hapMatrix[i] = (HapElem **) BioMalloc(sizeof(HapElem *) * nind);
		for (j=0;j<nind ;j++ ){
			gm[i]->genotype[j] = (char *) BioMalloc(sizeof(char) * 3); // only for dioploid species 
			hapMatrix[i][j] = (HapElem *) BioMalloc(sizeof(HapElem) * nind);
		}
	}
	
	// init 
	locsINFO->gMatrix = gm;
	locsINFO->hapMatrix = hapMatrix;	
}

/*
 *  Function: init_locus();
 *
 *  Purpose:  initialize gMatrix using a string 
 *
 *  Args:     str -- a record of loc, gm, locsINFO->nind 
 *
 *  Return:   None
 *
 */

static void init_locus(char **str, LocsINFO *locsINFO, int nind, int index){
	// get gMatrix and hapMatrix 
	GenoMatrix *gm = locsINFO->gMatrix[index];
	HapElem **hapMatrix = locsINFO->hapMatrix[index];

	// NOTE: marker ID,allocate memory here 
	gm->marker = bio_strcpy(str[0]);
	int i;
	// type 
	if ((strlen(str[1])!=5 && strlen(str[1])!=7)||!bio_regex_m(str[1], "x")) BioDie("segregation type error, %s, %s\n",gm->marker, str[1]);
	if (bio_regex_m(str[1],"^<")){
		for (i=0;i<5 ;i++ ) gm->type[i] = str[1][i+1];
		gm->type[5] = '\0';
	}else{
		for (i=0;i<5 ;i++ ) gm->type[i] = str[1][i];
		gm->type[5] = '\0';
	}
	// phase 
	if (strlen(str[2]) != 4 || str[2][0] != '{' || str[2][3] != '}') BioDie("phase error: %s, %s\n",gm->marker, str[2]);
	for (i=0;i<2 ;i++ ) gm->phase[i] = str[2][i+1];
	gm->phase[2] = '\0';
	// indi genotypes and haplotypes 
	for (i=0;i<nind ;i++ ){
		strcpy(gm->genotype[i],str[i+3]); // g
		calculate_one_haplo(gm->type,gm->phase,gm->genotype[i],hapMatrix[i]); // calculate haplotype 
	}
}
/*
 *  Function: calculate_one_haplo()
 *
 *  Purpose:  calculate haploSource using related info 
 *
 *  Args:     type, phase, genotype(per loci per individual)
 *
 *  Return:   None;
 *
 */
static inline void calculate_one_haplo(const char *type, const char *phase, const char *g, HapElem *h){
	// init 
	h->isMulti = 0;
	h->hapP = '-';
	h->hapM = '-';
	// type
	if (strcmp(type,"abxcd") == 0 || strcmp(type,"efxeg") == 0){
		h->type = 0;	
	}else if (strcmp(type,"hkxhk") == 0){
		h->type = 1;
	}else if (strcmp(type,"lmxll") == 0){
		h->type = 2;
	}else if (strcmp(type,"nnxnp") == 0){
		h->type = 3;
	}else{
		BioDie("error: unknow segregation type, %s", type);
	}

	if (strcmp(g,"--") == 0) return;  // miss observation
	
	// calculte if not missing data 
	char hapIndex[4][2] = {{0,0},{0,1},{1,0},{1,1}};
	char p[2],m[2],gamete[4][3];
	int i, j, n_gamete = 0;
	for (i=0; i<2; i++ ){
		p[i] = (phase[0] == '0') ? type[1-i] : type[i];
		m[i] = (phase[1] == '0') ? type[4-i] : type[3+i];
	}
	// gamete 
	for (i=0;i<2 ;i++ ){
		for (j=0;j<2 ;j++ ){
			gamete[n_gamete][0] = p[i];
			gamete[n_gamete][1] = m[j];
			gamete[n_gamete][2] = '\0';
			n_gamete++;
		}
	}

	// determine haplotype: here assume that genotype code is in the alphabetic order
	// for the sake of safety, reorder genotype code if it is not in the alphabetic order
	char swap,g_cpy[3];
	strcpy(g_cpy,g);
	if (g_cpy[0] > g_cpy[1]){
		swap = g_cpy[0];
		g_cpy[0] = g_cpy[1];
		g_cpy[1] = swap;
	}
	for (i=0;i<4 ;i++ ){
		if (gamete[i][0] > gamete[i][1]){
			swap = gamete[i][0];
			gamete[i][0] = gamete[i][1];
			gamete[i][1] = swap;
		}

		if (strcmp(gamete[i],g_cpy) == 0){
			if (p[0] == p[1]){ // case: nnxnp
				if (h->hapM == '-'){
					h->hapM = hapIndex[i][1];
					h->isMulti++;
				}
			}else if (m[0] == m[1]){ // case: lmxll
				if (h->hapP == '-'){
					h->hapP = hapIndex[i][0];
					h->isMulti++;
				}
			}else if (p[0] == m[0] && p[1] == m[1] && g_cpy[0] != g_cpy[1]){ // case: hkxhk->hk
				h->hapP = hapIndex[i][0];  // Note: the antithesis of haplotype can be calculated for a diploid species 
				h->hapM = hapIndex[i][1];
				h->isMulti++;
			}else{ // case: abxcd & efxeg & hkxhk->hh & hkxhk->kk
				h->hapP = hapIndex[i][0];
				h->hapM = hapIndex[i][1];
				h->isMulti++;
			}
		}
	}
}
//===================================
// Read pwd file
//===================================

/*
 *  Function: read_rec_LOD()
 *
 *  Purpose:  read pwd file, get pairwise data
 *
 *  Args:     pwd filename, LocsINFO *locsINFO
 *
 *  Return:   RecLOD **
 *
 */

static RecLOD ***read_rec_LOD(const char *filename, LocsINFO *locsINFO){
	BioStrIntHASH *locs = locsINFO->locs;
	// get the number of locs
	size_t n = locs->recordsnum(locs);
	// get gMatrix 
	GenoMatrix **gMatrix = locsINFO->gMatrix;

	// allocate memory for RecLOD
	RecLOD ***recLOD;
	recLOD = allocate_memory_for_RecLOD(n);

	// get filenames and open file handle
	FILE *fp;
	if( (fp = fopen(filename, "r")) == NULL ) BioDie("open rec LOD file %s is failed!", filename);

	// get locs in the linkage group
	size_t recordCount = 0, tmp ;
	char **str, *line;
	regmatch_t *reg;
	double rec=0, LOD=0 ;
	BioStrIntRecord *ra, *rb;
	int a, b;
	while ( (line=bio_readline(fp)) != NULL ){
		// check empty line
		bio_chomp(line);
		if( bio_is_blank_line(line) || bio_regex_m(line, "^;") || bio_regex_m(line, "^name")) {
			free(line); line=NULL;
			continue;
		}
		if (bio_regex_m(line, "^locus")){free(line);line=NULL;break;}
		// read locs pair
		// get ida idb rec and LOD
		tmp = bio_regex_gms(line, "[^[:space:]]+", &reg, 1, &str);
		if( tmp < 4) BioDie("read rec and LOD error: tmp(%d) < 4; %s\n", tmp, line);
		sscanf(str[2], "%lf", &rec);
		sscanf(str[3], "%lf", &LOD);
		// storing
		if( (ra = locs->findrecord(locs, str[0])) == NULL ) continue;
		if( (rb = locs->findrecord(locs, str[1])) == NULL ) continue;
		a = ra->value;
		b = rb->value;
		if( ((int)recLOD[0][a][b].sRec) != 1 ) BioDie("pwd file error: repeat record (%s, %s)", str[0], str[1]);
		if ( (strcmp(gMatrix[a]->type,"nnxnp") == 0 && strcmp(gMatrix[b]->type,"lmxll") == 0) || (strcmp(gMatrix[a]->type,"lmxll") == 0 && strcmp(gMatrix[b]->type,"nnxnp") == 0)) 
			BioDie("pwd file error: lmxll and nnxnp pair (%s, %s)",str[0], str[1]);
		if( !(rec<0.5) ) rec = 0.4999;	// 0.5 => 0.4999
		recLOD[0][a][b].sRec = rec;
		recLOD[0][b][a].sRec = rec;
		if (strcmp(gMatrix[a]->type,"nnxnp") != 0 && strcmp(gMatrix[b]->type,"nnxnp") != 0){ // not pair (nnxnp, nnxnp) series
			recLOD[0][a][b].pRec = rec;  // here p and m recombination frequency are equal to sexAver initially 
			recLOD[0][b][a].pRec = rec;  // still 1 signifies record not exists
		}
		if (strcmp(gMatrix[a]->type,"lmxll") != 0 && strcmp(gMatrix[b]->type,"lmxll") != 0){ // not pair (lmxll, lmxll) series
			recLOD[0][a][b].mRec = rec;
			recLOD[0][b][a].mRec = rec;
		}
		recLOD[0][a][b].iLOD = LOD;
		recLOD[0][b][a].iLOD = LOD;
		recLOD[0][a][b].pnRec = 0;
		recLOD[0][b][a].pnRec = 0;
		recLOD[0][b][a].mnRec = 0;
		recLOD[0][b][a].mnRec = 0;

		// free
		bio_free_2D_array((void **)str, tmp); str = NULL;
		if(reg==NULL){free(reg); reg = NULL;}
		if(line!=NULL){free(line); line=NULL;}
		recordCount++;

	}

	// check loc file and pwd file
	long r = n;
	if( recordCount != r*(r-1)/2 - locsINFO->lmxllSum * locsINFO->nnxnpSum ) 
		BioDie("loc and pwd file error: recordCount != r*(r-1)/2-lmxll*nnxnp");

	// over
	fclose(fp);
	return recLOD;
}
/*
 *  Function: allocate_memory_for_RecLOD()
 *
 *  Purpose:  allocate memory for pairwise matrix
 *
 *  Args:     locsINFO->nloc
 *
 *  Return:   RecLOD **
 *
 */
static RecLOD ***allocate_memory_for_RecLOD(size_t n){
	// check n
	if(n<3) BioDie("allocate_memory_for_RecLOD: n < 3");

	// allocate memory
	// recLOD[0] is used to store original pairwise data and updated pairwise data in gibbs sampling 
	// recLOD[1] is used to store sum of sampled pairwise matrix, so all elements are initialized to zero
	int i, j;
	RecLOD ***recLOD = (RecLOD ***) BioMalloc(sizeof(RecLOD **) * 2);
	recLOD[0] = (RecLOD **) BioMalloc(sizeof(RecLOD *) * n);
	recLOD[1] = (RecLOD **) BioMalloc(sizeof(RecLOD *) * n);  // RecLOD *
	for (i=0; i<n; i++){
		recLOD[0][i] = (RecLOD *) BioMalloc(sizeof(RecLOD) * n); // RecLOD
		recLOD[1][i] = (RecLOD *) BioMalloc(sizeof(RecLOD) * n); // RecLOD
	}

	// init RecLOD **;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			recLOD[0][i][j].sRec = 1;	// NOTE: 0<=rec<=0.5; so 1>0.5 mean the loc pair do not exist
			recLOD[0][i][j].pRec = 1;
			recLOD[0][i][j].mRec = 1;
			recLOD[0][i][j].iLOD = -1;	// nothing

			recLOD[1][i][j].sRec = 0;	// NOTE: 0<=rec<=0.5; so 1>0.5 mean the loc pair do not exist
			recLOD[1][i][j].pRec = 0;
			recLOD[1][i][j].mRec = 0;
			recLOD[1][i][j].pnRec = 0;
			recLOD[1][i][j].mnRec = 0;
			recLOD[1][i][j].iLOD = -1;	// nothing
		}
	}

	// over
	return recLOD;
}

//===================================
// Read fix order file 
//===================================
/*
 *  Function: init_fix_order()
 *
 *  Purpose:  init fix order info by reading fix order file 
 *
 *  Args:     LocsINFO->locsINFO, MapOpt->opt
 *
 *  Return:   None
 *
 */
static void init_fix_order(LocsINFO *locsINFO, MapOpt *opt) {
	// read fix order from file
	if( opt->fixOrderFile != NULL ) {
		read_fix_order_from_file(locsINFO, opt);
	// init the fix order to NULL
	}else{
		locsINFO->fixOrder = (int **) BioMalloc(sizeof(int *) * 1);;		// int **; a matrix
		locsINFO->fixOrder[0] = (int *) BioMalloc(sizeof(int) * locsINFO->nloc); 
		locsINFO->fixOrderLen = (int *) BioMalloc(sizeof(int) * 1);;	// int *; each order length or the number of loc in each order
		locsINFO->isfix = (int **) BioMalloc(sizeof(int *) * 1);
		locsINFO->isfix[0] = (int *) BioMalloc (sizeof(int) * locsINFO->nloc);
		memset(locsINFO->isfix[0], 0, sizeof(int)*locsINFO->nloc);
		locsINFO->fixOrderSum = 0;		// int ; the number of order
	}
}

/*
 *  Function: read_fix_order_from_file()
 *
 *  Purpose:  read fix order from file 
 *
 *  Args:     LocsINFO *locsINFO, MapOpt *opt
 *
 *  Return:   None
 */
static void read_fix_order_from_file(LocsINFO *locsINFO, MapOpt *opt) {
	// get the number of locs
	size_t n = locsINFO->nloc;
	// get locs
	BioStrIntHASH *locs = locsINFO->locs;

	// allocate memory for startOrder
	int order_count = 0, i, locs_count = 0;
	locsINFO->fixOrder = (int **) BioMalloc(sizeof(int *) * 2);	// default 
	locsINFO->fixOrderLen = (int *) BioMalloc(sizeof(int) * 2);
	locsINFO->isfix = (int **) BioMalloc(sizeof(int *) * 2);
	locsINFO->isfix[0] = (int *) BioMalloc (sizeof(int) * n);
	locsINFO->isfix[1] = (int *) BioMalloc (sizeof(int) * n);
	memset(locsINFO->isfix[0], 0, sizeof(int)*n);
	memset(locsINFO->isfix[1], 0, sizeof(int)*n); // zero 
	int *oneFixOrder = (int *) BioMalloc(sizeof(int) * n);			// one fix order
	int *order_check = (int *) BioMalloc(sizeof(int) * n);			// check if there are repeat locs

	// get filenames and open file handle
	FILE *fp;
	if( (fp = fopen(opt->fixOrderFile, "r")) == NULL ) BioDie("open fix order file %s is failed!", opt->fixOrderFile);

	// get locs in the linkage group
	size_t tmp, isStart=0;
	char **str, *line;
	regmatch_t *reg;
	BioStrIntRecord *r;
	while ( (line=bio_readline(fp)) != NULL ){
		// check empty line
		bio_chomp(line);
		if( bio_is_blank_line(line) ) { free(line); line=NULL; continue; }

		// check @
		if(bio_regex_m(line, "^@")){
			free(line); line=NULL;
			// store one fix order
			if(isStart == 1 && locs_count >= 3) {
//				if(locs_count < 3) BioDie("fix order error: locs_count < 3");		// the number of loc in one order must > 3
				if(order_count > MAXFIXORDER) BioDie("MAXFIXORDER is too large!");	// the number of fix order is too large
				locsINFO->fixOrder[order_count] = oneFixOrder;
				locsINFO->fixOrderLen[order_count] = locs_count;
				for (i = 0;i < locs_count ; i++) locsINFO->isfix[order_count][oneFixOrder[i]] = 1;
				order_count ++;
				oneFixOrder = (int *) BioMalloc(sizeof(int) * n);			// new a fix order
			}
			// update and init
			isStart = 1;
			locs_count = 0;
			for(i=0; i<n; i++) order_check[i] = 0;		// no 0; have 1; repeat 2;
			continue;
		}
		if(isStart==0) continue;

		// read locs
		tmp = bio_regex_gms(line, "[^[:space:]]+", &reg, 1, &str);
		if( tmp < 1) BioDie("unknown error: tmp(%d) < 1; %s\n", tmp, line);
		// store each locs
		for(i=0; i<tmp; i++){
//			if( (r = locs->findrecord(locs, str[i])) == NULL ) BioDie("fix order file error: id %s not exists!", str[i]);
			if( (r = locs->findrecord(locs, str[i])) == NULL ) continue;
			if(order_check[r->value] > 0) BioDie("repeat loc %s", str[i]);
			order_check[r->value] ++;
			oneFixOrder[locs_count] = r->value;
			locs_count ++;
		}

		// read one line over
		bio_free_2D_array((void **)str, tmp); str = NULL;
		if(reg==NULL){free(reg); reg = NULL;}
		if(line!=NULL){free(line); line=NULL;}
	}

	// process last fix order
	if(isStart == 1 && locs_count >= 3){
//		if(locs_count < 3) BioDie("fix order error: locs_count < 3");		// the number of loc in one order must > 3
		if(order_count > MAXFIXORDER) BioDie("MAXFIXORDER is too large!");	// the number of fix order is too large
		locsINFO->fixOrder[order_count] = oneFixOrder;
		locsINFO->fixOrderLen[order_count] = locs_count;
		for (i = 0;i < locs_count ; i++) locsINFO->isfix[order_count][oneFixOrder[i]] = 1;
		order_count ++;
	}

	// check and init start loc number
//	if(order_count < 1) BioDie("fix order file error: no fix order");  // 
	locsINFO->fixOrderSum = order_count;		// init
	if (order_count >=2) BioDie("Warning: Because multiple fix orders would seriously deteriorate the speed of algorithm,so they are not allowed\n"); 

	// over
	fclose(fp);
	free(order_check); order_check=NULL;
}
//===================================
// Read start order file 
//===================================

/*
 *  Function: init_start_order()
 *
 *  Purpose:  initialize start order info 
 *
 *  Args:     LocsINFO *locsINFO, MapOpt *opt
 *
 *  Return:   None
 */

static void init_start_order(LocsINFO *locsINFO, MapOpt *opt) {
	// read start order from file
	if( opt->startOrderFile != NULL ) {
		read_start_order_from_file(locsINFO, opt);
	}else{
		locsINFO->startOrder = NULL; 
	}
}

/* Function: read_start_order_from_file()
 *
 * Purpose:  reading start order file
 *
 * Args:     locsINFO - LocsINFO *
 *           opt - MapOpt *
 *
 * Return:   None
 *
 */
static void read_start_order_from_file(LocsINFO *locsINFO, MapOpt *opt) {
	// get the number of locs
	size_t n = locsINFO->nloc;
	// get locs
	BioStrIntHASH *locs = locsINFO->locs;

	// allocate memory for startOrder
	int order_count = 0, i;
	locsINFO->startOrder = (int *) BioMalloc(sizeof(int) * n);
	int *order_check = (int *) BioMalloc(sizeof(int) * n);		// check if there are repeat locs
	for(i=0; i<n; i++) order_check[i] = 0;				// no 0; have 1; repeat 2;

	// get filenames and open file handle
	FILE *fp;
	if( (fp = fopen(opt->startOrderFile, "r")) == NULL ) BioDie("open start order file %s is failed!", opt->startOrderFile);

	// get locs in the linkage group
	size_t tmp, isStart=0;
	char **str, *line;
	regmatch_t *reg;
	BioStrIntRecord *r;
	while ( (line=bio_readline(fp)) != NULL ){
		// check empty line
		bio_chomp(line);
		if( bio_is_blank_line(line) ) { free(line); line=NULL; continue; }

		// check @
		if(bio_regex_m(line, "^@")){
			free(line); line=NULL;
			isStart ++;
			if(isStart >= 2) BioDie("mutilple @ in the start order file: %s", opt->startOrderFile);
			continue;
		}
		if(isStart==0) continue;

		// read locs
		tmp = bio_regex_gms(line, "[^[:space:]]+", &reg, 1, &str);
		if( tmp < 1) BioDie("unknown error: tmp(%d) < 1; %s\n", tmp, line);
		// store each locs
		for(i=0; i<tmp; i++){
//			if( (r = locs->findrecord(locs, str[i])) == NULL ) BioDie("start order file error: id %s not exists!", str[i]);
			if( (r = locs->findrecord(locs, str[i])) == NULL ) continue;
			if(order_check[r->value] > 0) BioDie("repeat loc %s", str[i]);
			order_check[r->value] ++;
			locsINFO->startOrder[order_count] = r->value;
			order_count++;
		}

		// read one line over
		bio_free_2D_array((void **)str, tmp); str = NULL;
		if(reg==NULL){free(reg); reg = NULL;}
		if(line!=NULL){free(line); line=NULL;}
	}

	// check and init start loc number
	if(order_count < 2) BioDie("start order file error: order_count < 2");
	locsINFO->startOrderLen = order_count;

	// over
	fclose(fp);
	free(order_check); order_check=NULL;
}
/*
 *  Function: check_start_fix_conflict()
 *
 *  Purpose:  check confliction of fix order and start order
 *
 *  Args:     LocsINFO *locsINFO
 *
 *  Return:   None
 */
static void check_start_fix_conflict(LocsINFO *locsINFO){
	if (locsINFO->fixOrderSum == 0 || locsINFO->startOrder == NULL) return;
	int isConflict;  // 0: conflict, 1: ok
	// calling function check_map_order() to check,
	// this function is defined in siman.c
	isConflict = check_map_order(locsINFO,locsINFO->startOrder,locsINFO->startOrderLen); 
	
	if (!isConflict) BioDie("start order is conflict with fix order");
}

/* Function:check_map_order()
 *
 * Purpose: check map order with fix order
 *
 * Args:    locsINFO - LocsINFO *
 *          map - int *
 *          maplen - the number of locs in map
 *
 * Return:  0: confict with fix order; 1: ok
 *
 */
inline int check_map_order(LocsINFO *locsINFO, int *map, int maplen) {
	if(locsINFO->fixOrderSum == 0) return 1;

	// locsSum
	int locsSum = locsINFO->nloc;
	// fix order
	int orderSum = locsINFO->fixOrderSum;
	int *orderLen = locsINFO->fixOrderLen;
	int **fixOrder = locsINFO->fixOrder;

	// init mapIndex
	int i,j,a;
	int *mapIndex = (int *) BioMalloc(sizeof(int) * locsSum);
	for(i=0; i<locsSum; i++) mapIndex[i] = -1;		// -1: this loc does not exist
	for(i=0; i<maplen; i++) mapIndex[map[i]] = i;		// store loc index

	// init submap
	int *submap = (int *) BioMalloc(sizeof(int) * maplen);
	int submaplen = 0;

	// iteration check
	int result = 1;
	for(i=0; i<orderSum; i++){
		submaplen = 0;
		// get sub map according the one fix order
		for(j=0; j<orderLen[i]; j++){
			a = fixOrder[i][j];
			if(mapIndex[a] == -1) continue;
			submap[submaplen] = mapIndex[a];
			submaplen ++;
		}
		// check submaplen
		if(submaplen < 3) continue;	// submaplen must > 3
		// check submaplen
		if(submap[0] > submap[1]){	// check >>>>>>>
			for(j=1; j<submaplen-1; j++)
				if(submap[j] < submap[j+1]) {result=0; break; }	// check conflict
		}else{				// check <<<<<<
			for(j=1; j<submaplen-1; j++)
				if(submap[j] > submap[j+1]) {result=0; break; }	// check conflict
		}
		// check result
		if(result == 0) break;	// check conflict
	}

	// over and free
	free(mapIndex); mapIndex=NULL;
	free(submap); submap=NULL;
	return result;
}

//===================================
//  interface of output (loc, pwd, map)
//===================================


//#define MAPIO_TEST
#ifdef MAPIO_TEST
int main(void) {
	printf("Hello world\n");
	// init LocsINFO

	LocsINFO *locsINFO;
	char loc[] = "./td/demo.loc";
	char pwd[] = "./td/demo.pwd";

	MapOpt *opt = (MapOpt *) BioMalloc(sizeof(MapOpt) * 1);
	opt->locFile = bio_strcpy(loc);
	opt->pwdFile = bio_strcpy(pwd);
	opt->fixOrderFile = NULL;
	opt->startOrderFile = NULL;

	locsINFO = read_locs_info(opt);
	// get locs 
	BioStrIntHASH *locs = locsINFO->locs;

	printf("nloc = %d\n",locsINFO->nloc);
	printf("nind = %d\n",locsINFO->nind);
	

	printf("test over\n");
	return 0;
}

#endif
