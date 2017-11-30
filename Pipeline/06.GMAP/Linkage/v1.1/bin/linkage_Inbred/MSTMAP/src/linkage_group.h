/*
 *  linkage_group.h
 *  ApproxMap
 *
 *  Created by yonghui on 4/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __LINKAGE_GROUP_H__
#define __LINKAGE_GROUP_H__

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

using namespace std;

bool cmp(pair<double, pair<int,int> > element1, pair<double, pair<int,int> > element2);

class DF{
    public:
        virtual double CM(double rp) const {
            cout << "ERROR, specific DF should be used instead" << endl;
            assert(false);
            return 0;
        };
        virtual double RP(double cm) const {
            cout << "ERROR, specific DF should be used instead" << endl;
            assert(false);
            return 0;
        };
        virtual void print_df_name()  const {
            cout << "generic df" << endl;
        };
};

class DF_Haldane:public DF{
    public:
        virtual double CM(double rp) const {
            if (rp >= 0.5) rp = 0.5 - 0.000001;
            return -50.0 * log(1 - 2 * rp);
        };
        virtual double RP(double cm) const {
            return -0.5 * (exp(-cm / 50.0) - 1.0);
        };
        virtual void print_df_name() const {
            cout << "Haldane" << endl;
        };
};

class DF_Kosambi:public DF{
    public:
        virtual double CM(double rp) const {
            if (rp >= 0.5) rp = 0.5 - 0.000001;
            return 25.0 * log((1 + 2.0 * rp) / (1 - 2.0 * rp));
        };
        virtual double RP(double cm) const {
            return 0.5 * ((exp(cm / 25) - 1) / (exp(cm / 25) + 1));
        };
        virtual void print_df_name() const {
            cout << "Kosambi" << endl;
        };
};

class linkage_group {
    public:
        const vector<vector<double> >& get_pair_wise_distance() const;
        void return_order(vector<int>& out_order, 
                          double & _lowerbound, 
                          double & _upper_bound, 
                          double & _cost_after_initialization, 
                          vector<double> & _distances) const;
        void dump_common() const;
        void bad_genotypes(vector<pair<int,int> >& bad_genotypes) const;
        void dump_distance_matrix();
    protected:
    
        // produce pairwise_distance in cM
        void generate_distance_in_cM(vector<vector<double> >& distance_in_cM);
        void generate_distance_in_ML(vector<vector<double> >& distance_in_ML);
        
        bool detect_bad_data;
        ObjFunc objective_function;
        int number_of_bins;
        int number_of_individuals;
        
        // the distance is normalized to be the expected number of cross-overs per meiosis 
        // given the number of individuals 
        // namely, d_{i,j} = number_of_individuals * r
        // where r is the recombination probability
        vector<vector<double> > pair_wise_distances;
        vector<pair<int,int> > missing_data;
        vector<int> bin_sizes;
        
        // added by yonghui on March 8th, 2008. This data structure is updated 
        // after each iteration
        vector<pair<int,int> > suspicious_data;

        vector<int> current_order; //concatenation of the markers in the order of the path
        vector<int> MST; // the current MST
        double MST_lower_bound;
        double current_upper_bound;
        double cost_after_initialization;
        
        // the wrapper of distance functions (i.e. haldane or kosambi)
        // this is just a reference to an object owned by another object (an genetic_map instance)
        // no need to delete this object
        DF* df;
};


#endif


