#ifndef MAPGUILD_H
#define MAPBUILD_H

#include "map.h"
#include "mapio.h"
#include "mapopt.h"
#include "siman.h"
#include "spatsam.h"
#include "gibbssam.h"
#include "gibbsutil.h"
#include "metrohast.h"

// extern functions 
extern MapOpt *parse_command_line(int argc, char *argv[]);
extern LocsINFO *read_loc_pwd_info(MapOpt *opt);
extern void map_building(LocsINFO *locsINFO, MapOpt *opt);

extern void map_building_spatial_sample(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);
extern void map_building_start_sample(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);
extern void update_T_with_order(int *order, int orderLen, TSPMSET *tspm);
extern void map_sort_one_sample (LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);
extern void fix_next_sample_first_round(LocsINFO *locsINFO, TSPMSET *tspm, int flag);
extern void merge_start_fix_order(LocsINFO *locsINFO, TSPMSET *tspm);
extern TSPMSET *init_tspm(LocsINFO *locsINFO, MapOpt *opt);
extern void output_result(LocsINFO *locsINFO, TSPMSET *tspm, MapOpt *opt);
extern void output_map_file(LocsINFO *locsINFO, int *map,double *mapdist, int maplen, FILE *fp);


#endif
