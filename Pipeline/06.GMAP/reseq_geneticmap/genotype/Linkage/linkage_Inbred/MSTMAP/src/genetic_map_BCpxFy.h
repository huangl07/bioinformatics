/*
 *  genetic_map_RIL.h
 *  ApproxMap
 *
 *  Created by yonghui on 12/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __GENETIC_MAP_BCPXFY_H__
#define __GENETIC_MAP_BCPXFY_H__

#include "genetic_map.h"
#include "linkage_group_BCpxFy.h"

class genetic_map_BCpxFy: public genetic_map {
    public:
        genetic_map_BCpxFy():genetic_map(){};
      
        virtual void generate_map();
    private:
        void gen_raw_prob_data();
        void calculate_pair_wise_distance();

        /*generate a linkage group*/
        linkage_group_BCpxFy* construct_linkage_group(int group_id);
        linkage_group_BCpxFy* construct_linkage_group_whole_map();
        
        // ----------------------------
        // private data members section
        // ----------------------------
        int BC_generation_index_;
        int F_generation_index_;
        vector<vector<allel_state> > raw_prob_data_;
};
#endif

