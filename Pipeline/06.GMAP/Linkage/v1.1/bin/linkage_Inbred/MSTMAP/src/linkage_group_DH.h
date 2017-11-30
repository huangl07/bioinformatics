/*
 *  linkage_group.h
 *  ApproxMap
 *
 *  Created by yonghui on 4/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __LINKAGE_GROUP_DH_H__
#define __LINKAGE_GROUP_DH_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <cassert>
#include <queue>
#include <utility>
#include <set>
#include <algorithm>
#include "constants.h"
#include "MSTOpt.h"
#include "linkage_group.h"


using namespace std;

class linkage_group_DH: public linkage_group {
    private:
        vector<vector<float> > raw_data;
        
        // 0: if it is the same as the original input
        // 1: if it is missing
        // 2 and up: if it is deleted 
        vector<vector<int> > data_status;
        int iteration_number;
        
        vector<double> suspicious_data_backup;
        
        /*Calculate the pair_wise distance*/
        void calculate_pair_wise_distance();
        void calculate_pair_wise_distance_initialize();
        int detect_bad_markers();
        void revert_suspicious_data();
        
        void estimate_missing_data();
        // a supportive function to be called by estimate_missing_data
        double estimate_one_missing(const vector<int>& rev_perm, int marker_id, int indi_id);
        
    public:
        linkage_group_DH(int _number_of_bins, 
                         int _number_of_individuals,
                         bool _detect_bad_data,
                         ObjFunc _objective_function,
                         DF* _df,
                         const vector<vector<float> >& _raw_data, 
                         const vector<int>& _current_order, 
                         const vector<pair<int,int> >& _missing_data,
                         const vector<int>& _bin_sizes);
        ~linkage_group_DH();
        
        void dump() const;
    
        void order_markers();
};


#endif


