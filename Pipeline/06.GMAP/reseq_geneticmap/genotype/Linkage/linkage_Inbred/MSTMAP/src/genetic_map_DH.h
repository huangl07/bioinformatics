/*
 *  single_mapping_population_raw_data.h
 *  ApproxMap
 *
 *  Created by yonghui on 4/7/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __GENETIC_MAP_DH_H__
#define __GENETIC_MAP_DH_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include "constants.h"
#include <queue>
#include "linkage_group_DH.h"
#include "genetic_map.h"


using namespace std;

class genetic_map_DH: public genetic_map {
    private:
    	vector<pair<string, string> > suspicious_data;
        /*calculate the pair-wise distance*/
        void calculate_pair_wise_distance();
 
        /*generate a linkage group*/
        linkage_group_DH* construct_linkage_group(int group_id);
        linkage_group_DH* construct_linkage_group_whole_map();
        void print_suspicious_data();
    public:
        genetic_map_DH():genetic_map(){};
        ~genetic_map_DH();
        virtual void generate_map();
        void print_double_cross_overs();
};

#endif

