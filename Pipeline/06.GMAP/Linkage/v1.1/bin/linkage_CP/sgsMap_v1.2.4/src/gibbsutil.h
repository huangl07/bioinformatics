#ifndef GIBBSUTIL_H
#define GIBBSUTIL_H

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "bioerror.h"
#include "biomemory.h"
#include "biostring.h"
#include "map.h"
#include "gibbssam.h"


typedef double (*MapFun)(const double);
extern MapFun mapfun[2];

// extern functions
extern inline void update_rec_matrix(MCEM *mcem, int *map, int maplen, int n);
extern inline void calculate_pairwise_rec_info(GibbsSAM *gSam);
extern inline void calculate_singleton_and_rn(GibbsSAM *gSam, double *ps, double *ms, double *pn, double *mn);
extern inline void generate_next_sample(GibbsSAM *gSam);
extern inline void sampling_one_observation(GibbsSAM *gSam, Coordinate *ss);
extern inline void sampling_partial_informative(GibbsSAM *gSam, Coordinate *ss);
extern inline void sampling_miss_observation(GibbsSAM *gSam, Coordinate *ss);
extern inline void sample_anchor_haplo(HapElem ***hm, RecLOD **rm, Coordinate *ss, int pLen, int mLen, gsl_rng *r);
extern inline void sample_paternal_haplo(HapElem ***hm, RecLOD **rm, Coordinate *ss, int pLen, gsl_rng *r);
extern inline void sample_maternal_haplo(HapElem ***hm, RecLOD **rm, Coordinate *ss, int mLen, gsl_rng *r);
extern inline double calculate_paternal_condition_prob(HapElem ***hm, RecLOD **rm, Coordinate *ss, char hap, int pLen, gsl_rng *r);
extern inline double calculate_maternal_condition_prob(HapElem ***hm, RecLOD **rm, Coordinate *ss, char hap, int mLen, gsl_rng *r);
extern inline void reset_sum_recMatrix(RecLOD **rm, int *map, int maplen);
extern inline void init_miss_observation(GibbsSAM *gSam);
extern inline void init_map_and_sample_space(GibbsSAM *gSam);
extern inline int get_pre_sexAver_loc(GenoMatrix **gMatrix, int *map, int p, int flag);
extern inline int get_next_sexAver_loc(GenoMatrix **gMatrix, int *map, int p, int maplen, int flag);

extern inline int get_pre_adjacent_loc(RecLOD **rm, GenoMatrix **gMatrix, HapElem ***hm, int *map, int p, int indi, int nind, int *bin, int *binLen, int flag);
extern inline int get_next_adjacent_loc(RecLOD **rm, GenoMatrix **gMatrix, HapElem ***hm, int *map, int p, int maplen, int indi, int nind, int *bin, int *binLen, int flag);

extern inline double haldane(const double rec);
extern inline double kosambi(const double rec);
extern inline double inverseHaldane(double d);
extern inline double inverseKosambi(double d);

extern inline double calculate_conditional_prob_middle(char haplo, char preHaplo, char sufHaplo, double preRec, double sufRec);
extern inline double calculate_conditional_prob_ht(char haplo, char adjHap, double adjRec);

extern inline void print_gibbs_info_into_log(GibbsSAM *gSam, FILE *fp);
extern inline void print_gibbs_stat(GibbsSAM *gSam, FILE *fp);
extern inline void print_best_map(GibbsSAM *gSam, FILE *fp);

#endif
