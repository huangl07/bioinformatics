/*
 *  linkage_group_RIL.h
 *  ApproxMap
 *
 *  Created by yonghui on 12/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __LINKAGE_GROUP_BCPXFY_H__
#define __LINKAGE_GROUP_BCPXFY_H__

#include <cassert>
#include <vector>
#include "linkage_group.h"
#include "MSTOpt.h"
#include "def.h"

using namespace std;

/*
 * BCpxFy: advanced backcross inbred line family: starting from the BC1 repeatedly backcrossing to the same parent
 * (as used for the BC1) of each individual resulting in a single offspring per individual, followed by selfing with
 * single seed descent; the backcross parent p and the generation x and y must be specified: p = A or B, x is the number
 * of backcrosses including the one for creating the BC1: 1 <= x <= 99, y is the number of selfings: 0 <= y <= 99, BCa1F0
 * is equivalent to BC1.
 */

class BCpxFy_dist_cal{
    public:
        BCpxFy_dist_cal(int _BC_generation_index, int _F_generation_index, const vector<allel_state>& _marker1, const vector<allel_state>& _marker2)
            :BC_generation_index_(_BC_generation_index),
            F_generation_index_(_F_generation_index),
            marker1_(_marker1),
            marker2_(_marker2){
            num_of_eff_individuals_ = 0;
            assert(marker1_.size() == marker2_.size());
            num_of_individuals_ = marker1_.size();
            count_class();
            upper_bound();
            lower_bound();
        };
        double Dist() const {
            double opt_delta = find_opt_delta();
            return opt_delta * num_of_individuals_;
        };
    private:
        void count_class();
        void upper_bound();
        void lower_bound();
        
        void expected_CDEFG(double delta,
                            double& C,
                            double& D,
                            double& E,
                            double& F,
                            double& G) const; 
                            
        double squared_error(double delta) const;
        
        double find_opt_delta() const;
        
        // ----------------------------
        // private data members section
        // ----------------------------
        
        // p = A|B, 1 <= x <= 99, 0 <= y <= 99
        int BC_generation_index_;
        int F_generation_index_;
        // char p_;
        
        int num_of_individuals_;
        int num_of_eff_individuals_;
        const vector<allel_state>& marker1_;
        const vector<allel_state>& marker2_;
        
        /*    Status       Diplotype               Proportion
         *      1          AABB, aabb              CC_
         *      2          AAbb, aaBB              DD_
         *      3          AAbb,AaBB,Aabb,aaBb     EE_
         *      4          [AB][ab]                FF_
         *      5          [Ab][aB]                GG_     
         *    NOTE: O(4,5) = FG_ = FF_ + GG_ 
         */
        double CC_;
        double DD_;
        double EE_;
        double FG_;
        double delta_upper_bound_;
        double delta_lower_bound_;
};

/*
 * BCpxFy: advanced backcross inbred line family: starting from the BC1 repeatedly backcrossing to the same parent
 * (as used for the BC1) of each individual resulting in a single offspring per individual, followed by selfing with
 * single seed descent; the backcross parent p and the generation x and y must be specified: p = A or B, x is the number
 * of backcrosses including the one for creating the BC1: 1 <= x <= 99, y is the number of selfings: 0 <= y <= 99, BCa1F0
 * is equivalent to BC1.
 */
class linkage_group_BCpxFy: public linkage_group {
    public:
        linkage_group_BCpxFy(int _number_of_bins, 
                          int _number_of_individuals,
                          int _BC_generation_index,
                          int _F_generation_index,
                          DF* _df,
                          const vector<vector<allel_state> >& _raw_data, 
                          const vector<int>& _current_order, 
                          const vector<pair<int,int> >& _missing_data);
        ~linkage_group_BCpxFy();
        
        void dump() const;
    
        void order_markers();
        
    private:

        /*Calculate the pair_wise distance*/
        void calculate_pair_wise_distance();
        
        void estimate_missing_data();
        
        // ----------------------------
        // private data members section
        // ----------------------------
        vector<vector<allel_state> > raw_data;
        // p = A|B, 1 <= x <= 99, 0 <= y <= 99
        int BC_generation_index_;
        int F_generation_index_; 
        char p_;                 

};

#endif

