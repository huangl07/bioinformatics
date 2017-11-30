/*
 * Get maximum likelihood mapping options 
 *
 */

// copyright 2012, Ma Chouxian <macx@biomarker.com.cn> <machouxian@mail.nankai.edu.cn>

#include "mapopt.h"

// static methos 
static void map_usage();
static void opt_default(MapOpt *opt);

MapOpt *get_map_opt(int argc, char *argv[]){
	// init opt
	// allocate memory 
	MapOpt *opt = (MapOpt *) BioMalloc(sizeof(MapOpt) * 1);
	opt->SS_rec_thresholds = (double *) BioMalloc(sizeof(double)*6);
	// default option 
	opt_default(opt);
	// define long option struct
	static struct option long_options[] =
	{
		{"loc", 1, 0, 'L'},
		{"pwd", 1, 0, 'P'},
		{"strts", 1, 0, 's'},
		{"fix", 1, 0, 'f'},
		{"mf", 1, 0, 't'},
		{"s1", 1, 0, '1'},
		{"s2", 1, 0, '2'},
		{"s3", 1, 0, '3'},
		{"s4", 1, 0, '4'},
		{"s5", 1, 0, '5'},
		{"round", 1, 0, 'r'},
		{"SA_chain_length", 1, 0, 'c'},
		{"SA_initial_acceptance", 1, 0, 'p'},
		{"SA_cooling_control", 1, 0, 'C'},
		{"SA_unimproved_chain_max", 1, 0, 'm'},
		{"Gibbs_burn_in", 1, 0, 'b'},
		{"Gibbs_ncycle_MCEM", 1, 0, 'M'},
		{"Gibbs_chain_length", 1, 0, 'l'},
		{"Gibbs_sample_period", 1, 0, 'I'},
		{"PP_acceptance_control", 1, 0, 'A'},
		{"PP_burn_in", 1, 0, 'B'},
		{"PP_nsample", 1, 0, 'N'},
		{"PP_sample_period", 1, 0, 'S'},
		{0,0,0,0}
	};
	//===================================
	// Get option
	//===================================
	int opti;
	while ((opti = getopt_long_only(argc, argv,"k:h",long_options,NULL)) != -1){
		switch (opti){
		// key of output files 
		case 'k':opt->fKey = bio_strcpy(optarg); break;
		// loc file 
		case 'L':opt->locFile = bio_strcpy(optarg); break;
		// pwd file 
		case 'P':opt->pwdFile = bio_strcpy(optarg); break;
		// start order file 
		case 's':opt->startOrderFile = bio_strcpy(optarg); break;
		// fix order file 
		case 'f':opt->fixOrderFile = bio_strcpy(optarg); break;
		// map function: 0 for Haldane and 1 for kosambi
		case 't':opt->mapFunc = atoi(optarg); break;
		// spatial sampling rec.freq. threshold 1
		case '1':sscanf(optarg, "%lf", &(opt->SS_rec_thresholds[0])); break;
		// spatial sampling rec.freq. threshold 1
		case '2':sscanf(optarg, "%lf", &(opt->SS_rec_thresholds[1])); break;
		// spatial sampling rec.freq. threshold 1
		case '3':sscanf(optarg, "%lf", &(opt->SS_rec_thresholds[2])); break;
		// spatial sampling rec.freq. threshold 1
		case '4':sscanf(optarg, "%lf", &(opt->SS_rec_thresholds[3])); break;
		// spatial sampling rec.freq. threshold 1
		case '5':sscanf(optarg, "%lf", &(opt->SS_rec_thresholds[4])); break;
		// nr. of map optimiztion rounds per sample
		case 'r':opt->MB_nopt = atoi(optarg); break;
		// chain length (with constant acc.prob)
		case 'c':opt->SA_nstep = atoi(optarg); break;
		// SA_initial_acceptance
		case 'p':sscanf(optarg, "%lf", &(opt->SA_initial_acceptance)); break;
		// SA_cooling_control
		case 'C':sscanf(optarg, "%lf", &(opt->SA_cooling_control)); break;
		// SA_unimproved_chain_max
		case 'm':opt->SA_unimproved_chain_max = atoi(optarg); break;
		// Gibbs_burn_in
		case 'b':opt->Gibbs_burn_in = atoi(optarg); break;
		// Gibbs_ncycle_MCEM
		case 'M':opt->Gibbs_ncycle_MCEM = atoi(optarg); break;
		// Gibbs_chain_length
		case 'l':opt->Gibbs_chain_length = atoi(optarg); break;
		// Gibbs_sample_period
		case 'I':opt->Gibbs_sample_period = atoi(optarg); break;
		// PP_acceptance_control
		case 'A':sscanf(optarg, "%lf", &(opt->PP_acceptance_control)); break;
		// PP_burn_in
		case 'B':opt->PP_burn_in = atoi(optarg); break;
		// PP_nsample
		case 'N':opt->PP_nsample = atoi(optarg); break;
		// PP_sample_period
		case 'S':opt->PP_sample_period = atoi(optarg); break;
		// help
		case 'h': map_usage(); break;
		default : map_usage();
		}
	} 
	// check opt
	
	// loc file,pwd file and map file are essentially required
	if (opt->locFile == NULL || opt->pwdFile == NULL || opt->fKey == NULL)
		map_usage();
	//===================================
	// check map building parameters
	//===================================
	// check spatial sampling threshold
	int i, j;
	for (i=0;i<4 ;i++ ){
		for (j=i+1;j<5 ;j++ ){
			if (opt->SS_rec_thresholds[i] < opt->SS_rec_thresholds[j])
				BioDie("Option Error: spatial sampling threshod %d <= spatial sampling threshold %d", i, j);
		}
	}
	for (i=0;i<5 ;i++ )
		if (opt->SS_rec_thresholds[i] < 0) BioDie("Option Error: spatial sampling threshod %d: %lf < 0",i, opt->SS_rec_thresholds[i]);
	
	// MB_nopt
	if (opt->MB_nopt <= 0 || opt->MB_nopt > 20) 
		BioDie("Option Error: -MB_nopt must be in range: (0,20]");
	//===================================
	// check simulating annealing 
	//===================================
	// SA_nstep
	if (opt->SA_nstep < 100)
		BioDie("Option Error: -SA_nstep must be in range: [100,+inf]");
	// SA_initial_acceptance
	if (opt->SA_initial_acceptance <= 0 || opt->SA_initial_acceptance > 1)
		BioDie("Option Error: -SA_initial_acceptance must be in range: (0,1]");
	// SA_cooling_control
	if (opt->SA_cooling_control <= 0 || opt->SA_cooling_control > opt->SA_initial_acceptance)
		BioDie("Option Error: -SA_cooling_control must be in range: (0,%lf]",opt->SA_initial_acceptance);
	// SA_unimproved_chain_max
	if (opt->SA_unimproved_chain_max <= 0)
		BioDie("Option Error: -SA_unimproved_chain_max must be in range: (0,+inf]");

	//===================================
	// check Gibbs option
	//===================================
	//Gibbs_burn_in
	if (opt->Gibbs_burn_in < 100)
		BioDie("Option Error: -Gibbs_burn_in must be in range: [100,+inf]");
	//Gibbs_ncycle_MCEM
	if (opt->Gibbs_ncycle_MCEM <= 0)
		BioDie("Option Error: -Gibbs_ncycle_MCEM must be in range: (0,+inf]");
	//Gibbs_chain_length
	if (opt->Gibbs_chain_length <= 0)
		BioDie("Option Error: -Gibbs_chain_length must be in range: (0,+inf]");
	//Gibbs_sample_period
	if (opt->Gibbs_sample_period <= 0 || opt->Gibbs_sample_period > opt->Gibbs_chain_length)
		BioDie("Option Error: -Gibbs_sample_period must be in range: (0,%d]",opt->Gibbs_chain_length);
	//Gibbs_global_sampling
	if (opt->Gibbs_global_sampling != 0 && opt->Gibbs_global_sampling != 1)
		BioDie("Option Error: -Gibbs_global_sampling must be in {0,1}");
	//===================================
	// check plausible position option 
	//===================================
	// PP_acceptance_control;
	if (opt->PP_acceptance_control <= 0 || opt->PP_acceptance_control > 1)
		BioDie("Option Error: -PP_acceptance_control must be in range: (0,1]");
	// PP_burn_in;
	if (opt->PP_burn_in < 100)
		BioDie("Option Error: -PP_burn_in must be in range: [100,+inf]");
	// PP_nsample;
	if (opt->PP_nsample <= 0)
		BioDie("Option Error: -PP_nsample must be in range: (0,+inf]");
	// PP_sample_period;
	if (opt->PP_sample_period <= 0 || opt->PP_sample_period > opt->PP_nsample)
		BioDie("Option Error: -Gibbs_sample_period must be in range: (0,%d]",opt->PP_nsample);

	// check other options
	if((opt->mapFunc  != 0) && (opt->mapFunc != 1) )
		BioDie("Option Error: map function (-m) must be 1 or 2");

	//===================================
	// init non-command line options
	//===================================
	// map
	if (opt->sexAverMapFile == NULL)
		opt->sexAverMapFile = bio_strcat(opt->fKey, ".sexAver.map");
	if (opt->maleMapFile == NULL)
		opt->maleMapFile = bio_strcat(opt->fKey, ".male.map");
	if (opt->femaleMapFile == NULL)
		opt->femaleMapFile = bio_strcat(opt->fKey, ".female.map");
	// pwd 
	if (opt->sexAverPwdFile == NULL)
		opt->sexAverPwdFile = bio_strcat(opt->fKey, ".sexAver.pwd");
	if (opt->malePwdFile == NULL)
		opt->malePwdFile = bio_strcat(opt->fKey, ".male.pwd");
	if (opt->femalePwdFile == NULL)
		opt->femalePwdFile = bio_strcat(opt->fKey, ".female.pwd");
	// plausible position matrix 
	if (opt->malePlausiblePosFile == NULL){
		opt->malePlausiblePosFile = bio_strcat(opt->fKey, ".male.ppm");
	}
	if (opt->femalePlausiblePosFile == NULL){
		opt->femalePlausiblePosFile = bio_strcat(opt->fKey, ".female.ppm");
	}
	// singleton stat file
	if (opt->singletonFile == NULL){
		opt->singletonFile = bio_strcat(opt->fKey, ".singleton.txt");
	}
	// haplotype File
	if (opt->haplotypeFile == NULL){
		opt->haplotypeFile = bio_strcat(opt->fKey, ".haplotype.txt");
	}
	// log file 
	if( opt->logFile == NULL )
		opt->logFile = bio_strcat(opt->fKey, ".log");

	// open log file
	if( (opt->logfp = fopen(opt->logFile, "w")) == NULL )
		BioDie("open or create log file %s is failed!", opt->logFile);

	// over
	return opt;
}

void write_map_opt_into_logfile(MapOpt *opt) {

	struct tm *startTime = (struct tm *) BioMalloc(sizeof(struct tm) * 1);
	time_t timep;
	time (&timep);
	startTime = localtime(&timep);
	char *t = (char *) BioMalloc (sizeof(char)* 26);
	t=asctime(startTime);

	// software name \\ version \\ start time \\ license 
	fprintf(opt->logfp,"sgsMap v1.0.0\t%s\n",t);
	fprintf(opt->logfp,"license to :  Ma Chouxian <macx@biomarker.com.cn>\n\n");

	// input files 
	fprintf(opt->logfp,"%-50s  %s\n", "locus genotype file:", opt->locFile);
	fprintf(opt->logfp,"%-50s  %s\n", "pwd file:", opt->pwdFile);
	// output files
	fprintf(opt->logfp,"%-50s  %s\n\n", "key of output files:", opt->fKey);
	// fix and start order 
	if(opt->startOrderFile != NULL) fprintf(opt->logfp,"%-50s %s\n\n", "start order file:", opt->startOrderFile);
	if(opt->fixOrderFile != NULL) fprintf(opt->logfp,"%-50s  %s\n\n", "fix order file:", opt->fixOrderFile);
	
	// calculate options
	if(opt->mapFunc == 1)
		fprintf(opt->logfp,"%-50s  Kosambi\n\n","map function:");
	else
		fprintf(opt->logfp,"%-50s  Handale\n\n","map function:");
	// map building 
	fprintf(opt->logfp,"map building parameters:\n");
	fprintf(opt->logfp,"  %-48s  %.3f\n", "spatial sampling rec.feq. threshold 1:", opt->SS_rec_thresholds[0]);
	fprintf(opt->logfp,"  %-48s  %.3f\n", "spatial sampling rec.feq. threshold 2:", opt->SS_rec_thresholds[1]);
	fprintf(opt->logfp,"  %-48s  %.3f\n", "spatial sampling rec.feq. threshold 3:", opt->SS_rec_thresholds[2]);
	fprintf(opt->logfp,"  %-48s  %.3f\n", "spatial sampling rec.feq. threshold 4:", opt->SS_rec_thresholds[3]);
	fprintf(opt->logfp,"  %-48s  %.3f\n", "spatial sampling rec.feq. threshold 5:", opt->SS_rec_thresholds[4]);
	fprintf(opt->logfp,"  %-48s  %d\n", "nr. of map optimization rounds per sample:", opt->MB_nopt);
	// simulating annealing 
	fprintf(opt->logfp,"simulated annealing parameters (map order optimization):\n");
	fprintf(opt->logfp,"  %-48s  %d\n","chain length (with constant acc.prob.):",opt->SA_nstep);
	fprintf(opt->logfp,"  %-48s  %.3f\n","initial acceptance probability:",opt->SA_initial_acceptance);
	fprintf(opt->logfp,"  %-48s  %.3f\n","cooling control parameter:",opt->SA_cooling_control);
	fprintf(opt->logfp,"  %-48s  %d\n","stop after # chains without improvement:",opt->SA_unimproved_chain_max);
	// Gibbs sampling 
	fprintf(opt->logfp,"Gibbs sampling parameters (MC-EM multipoint ML estimation of rec. frequencies):\n");
	fprintf(opt->logfp,"  %-48s  %d\n","length of burn-in chain:",opt->Gibbs_burn_in);
	fprintf(opt->logfp,"  %-48s  %d\n","nr. of EM cycles:",opt->Gibbs_ncycle_MCEM);
	fprintf(opt->logfp,"  %-48s  %d\n","chain length per EM cycle:",opt->Gibbs_chain_length);
	fprintf(opt->logfp,"  %-48s  %d\n","period between rec.freq. matrix samples:",opt->Gibbs_sample_period);

	fprintf(opt->logfp,"======================================================\n");

}

static void opt_default(MapOpt *opt){
	// map bulding option: default 
	opt->SS_rec_thresholds[0] = 0.1;
	opt->SS_rec_thresholds[1] = 0.05;
	opt->SS_rec_thresholds[2] = 0.03;
	opt->SS_rec_thresholds[3] = 0.02;
	opt->SS_rec_thresholds[4] = 0.01;
	opt->SS_rec_thresholds[5] = 0.0;
	opt->MB_nopt = 3;
	
	// SA option: default 
	opt->SA_nstep = 1000;
	opt->SA_initial_acceptance = 0.25;
	opt->SA_cooling_control = 0.001;
	opt->SA_unimproved_chain_max = 10000;

	// Gibbs option: default 
	opt->Gibbs_burn_in = 10000;
	opt->Gibbs_ncycle_MCEM = 4;
	opt->Gibbs_chain_length = 1000;
	opt->Gibbs_sample_period = 5;
	opt->Gibbs_global_sampling = 0;

	// Metropolis Hasting option: defaut 
	opt->PP_acceptance_control =  1.0;
	opt->PP_burn_in =  10000;
	opt->PP_nsample =  1000;
	opt->PP_sample_period =  1000;
	
	opt->mapFunc = 0;	// default map function is Haldane
	
	// IO option
	opt->locFile = NULL;
	opt->pwdFile = NULL;
	opt->startOrderFile = NULL;
	opt->fixOrderFile = NULL;
	opt->fKey = NULL;
	opt->sexAverMapFile = NULL;
	opt->maleMapFile = NULL;
	opt->femaleMapFile = NULL;
	opt->sexAverPwdFile = NULL;
	opt->malePwdFile = NULL;
	opt->femalePwdFile = NULL;
	opt->malePlausiblePosFile = NULL;
	opt->femalePlausiblePosFile = NULL;
	opt->singletonFile = NULL;
	opt->haplotypeFile = NULL;
	opt->logFile = NULL;
}

static void map_usage(){

	printf("NAME \n\tsgsMap - Get linear arangement of loci in a linkage group and construct genetic maps of CP population\n\n");
	
	printf("SYNOPSIS \n\tsgsMap -L <LOC FILE> -P <PWD FILE> -k <key> [OPTIONS]\n\n");
	
	printf("DESCRIPTION \n\tThis program uses a combination of several techniques to order loci and compute their mutual distance:\n");
	printf("\tsimulated annealing, Gibbs sampling and spatial sampling.This method is introduced from joinmap.\n");
	printf("\tGibbs sampling is used to estimate multipoint recombination frequencies that can be used to calculate the likelihoods.\n");
	printf("\tSimulated annealing searches for order that has the maximum likelihood.\n");
	printf("\tSpatial sampling is a technique that is needed to prevent getting trapped at local optima\n");
	printf("\trather than arriving at the global optimum solution due to missing genotype information and genotyping errors.\n\n");
	
	printf("OPTIONS \n");
	printf("\t%-30s\t%-10s\t%s\n\n","-k","<str>","Key of output file, required");
	printf("\t%-30s\t%-10s\t%s\n\n","-L, -loc/--loc","<file>","Locus genotype file(loc), required");
	printf("\t%-30s\t%-10s\t%s\n\n","-P, -pwd/--pwd","<file>","PairWise data file,(pwd), required");
	printf("\t%-30s\t%-10s\t%s\n\n","-s, -strts/--strts","<file>","Start order file ,optional");
	printf("\t%-30s\t%-10s\t%s\n\n","-f, -fix/--fix","<file>","Fix order file ,optional");
	
	printf("\tMap building parameters:\n");
	printf("\t  %-28s\t%-10s\t%s\n\n","-1, -s1/--s1","<double>","spatial sampling rec.feq. threshold 1 ,optional, default[0.1]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-2, -s2/--s2","<double>","spatial sampling rec.feq. threshold 2 ,optional, default[0.05]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-3, -s3/--s3","<double>","spatial sampling rec.feq. threshold 3 ,optional, default[0.03]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-4, -s4/--s4","<double>","spatial sampling rec.feq. threshold 4 ,optional, default[0.02]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-5, -s5/--s5","<double>","spatial sampling rec.feq. threshold 5 ,optional, default[0.01]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-r, -round/--round","<int>","Nr. of map optimization rounds per sample, optional, default [3]");
	
	printf("\tMap order optimization:\n");
	printf("\t  %-28s\t%-10s\t%s\n\n","-c, -SA_chain_length","<int>","Chain length(with constant acc.prob), optional, default [1000]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-p, -SA_initial_acceptance","<double>","Initial acceptance probability, optional, default [0.25]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-C, -SA_cooling_control","<double>","Cooling control parameter,optional, default [0.001]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-m, -SA_unimproved_chain_max","<int>","Stop after # chains without improvement, optional, default [10000]");
	
	printf("\tMultipoint estimation of recombination frequencies:\n");
	printf("\t  %-28s\t%-10s\t%s\n\n","-b, -Gibbs_burn_in","<int>","Length of burn-in chain, optional, default [10000]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-M, -Gibbs_ncycle_MCEM","<int>","Nr. of Monte Carlo EM cycles, optional, default [4]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-l, -Gibbs_chain_length","<int>","Chain length per Monte Carlo EM cycle, optional, default [1000]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-I, -Gibbs_sample_period","<int>","sampling period for rec.freq matrix samples, optional, default [5]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-t, -mf/--mf","<int>","Map function, 0 for Haldane and 1 for Kosambi, optional, default [0]");
	
	printf("\tPlausible map positions:\n");
	printf("\t  %-28s\t%-10s\t%s\n\n","-A, -PP_acceptance_control","<double>","Acceptance control parameter, optional, default [1]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-B, -PP_burn_in","<int>","Burn-in chain length, optional, default [10000]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-N, -PP_nsample","<int>","Number of map samples to draw, optinal, default [1000]");
	printf("\t  %-28s\t%-10s\t%s\n\n","-S, -PP_sample_period","<int>","Sampling period for map samples, optional, default [1000]");
	
	printf("\t-h\t\t\t\t\t\tHelp\n\n");
	
	printf("AUTHOR \n\tWritten by Ma Chouxian\n\n");
	printf("REPORTING BUGS \n\tReport bugs to <macx@biomarker.com.cn>\n\n");
	printf("COPYRIGHT \n\tThis is non-open source software.All copyright reserved.\n\n");
	
	exit(0);
}

//#define MAPOPT_TEST
#ifdef MAPOPT_TEST
int main(int argc,char *argv[]){
	
	MapOpt *opt = get_map_opt(argc, argv);
	
	write_map_opt_into_logfile(opt);
	printf("over\n");

	return 0;
}

#endif

