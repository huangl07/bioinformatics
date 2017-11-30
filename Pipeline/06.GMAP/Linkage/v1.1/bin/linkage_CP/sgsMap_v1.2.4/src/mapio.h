#ifndef MAPIO_H
#define MAPIO_H


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioio.h"
#include "biohash.h"
#include "bioerror.h"
#include "biomemory.h"
#include "biostring.h"
#include "bioregex.h"
#include "map.h"
#include "mapopt.h"


// extern function
extern LocsINFO *read_locs_info(MapOpt *opt);
extern inline int check_map_order(LocsINFO *locsINFO, int *map, int maplen);

#endif
