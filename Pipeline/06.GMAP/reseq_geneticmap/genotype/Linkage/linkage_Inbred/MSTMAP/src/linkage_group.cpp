/*
 *  linkage_group_DH.cpp
 *  ApproxMap
 *
 *  Created by yonghui on 4/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "linkage_group.h"
#include <limits>

bool cmp(pair<double, pair<int,int> > element1, pair<double, pair<int,int> > element2)
{
    return element1.first < element2.first;
};


const vector<vector<double> >& linkage_group::get_pair_wise_distance() const
{
    return pair_wise_distances;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void linkage_group::generate_distance_in_cM(vector<vector<double> >& distance_in_cM){
    distance_in_cM.resize(number_of_bins);
    for (int ii = 0; ii < number_of_bins; ii++) {
        distance_in_cM[ii].resize(number_of_bins);
    }
    for (int ii = 0; ii < number_of_bins; ii++) {
        for (int jj = 0; jj < number_of_bins; jj++) {
            double r = pair_wise_distances[ii][jj] / number_of_individuals;
            if (r >= 0.5) {
                r = r - ZERO_PLUS;
            }
            distance_in_cM[ii][jj] = df->CM(r);
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group::generate_distance_in_ML(vector<vector<double> >& distance_in_ML){
    distance_in_ML.resize(number_of_bins);
    for (int ii = 0; ii < number_of_bins; ii++) {
        distance_in_ML[ii].resize(number_of_bins);
    }
    for (int ii = 0; ii < number_of_bins; ii++) {
        for (int jj = 0; jj < number_of_bins; jj++) {
            double r = pair_wise_distances[ii][jj] / number_of_individuals;
            if (r >= 0.5) {
                r = r - ZERO_PLUS;
            }
            if (r == 0.0) {
                distance_in_ML[ii][jj] = 0.0;
            } else {
                distance_in_ML[ii][jj] = -(r * log(r) + (1 - r) * log(1 - r));
            }
        }
    }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group::dump_common() const {
    cout << "number of bins:" << number_of_bins << endl;
    cout << "number of individuals:" << number_of_individuals << endl;
    cout << "current_order:" << endl;
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        cout << current_order[ii] << ',' ;
    }    
    cout << endl;

    cout << "lowerbound: " << MST_lower_bound << " the upperbound:" << current_upper_bound << endl;
    
    /*Print out the information regarding the MST*/
    cout << "The MST:" << endl;
    for (int ii = 0 ; ii < number_of_bins ; ii++)
    {
        cout << MST[ii] << ',' ;
    }
    
    vector<int> tmp_count(number_of_bins, 0);
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        tmp_count[MST[ii]] = tmp_count[MST[ii]] + 1;
    }
    cout << endl;
    
    cout << "The indegree for each of the vertices: " << endl;
    for (int ii = 0 ; ii < number_of_bins; ii++)
    {
        cout << tmp_count[ii] << ',';
    }
    
    cout << endl;    

    cout << "df function:";
    df->print_df_name();
    cout << endl;
    cout << "the distance between consecutive pairs:" << endl;
    for (int ii = 1 ; ii < number_of_bins; ii++)
    {
        cout << pair_wise_distances[current_order[ii]][current_order[ii-1]] << ',';
    }
    cout << endl;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group::bad_genotypes(vector<pair<int,int> >& bad_genotypes) const{
    bad_genotypes.clear();
    for(int ii = 0; ii < suspicious_data.size(); ii++) {
        bad_genotypes.push_back(suspicious_data[ii]);
    } 
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group::dump_distance_matrix() {
    char buffer[10];
    cout << "distance matrix within linkage_group" << endl;
    cout << "matrix dimension:" << pair_wise_distances.size() << endl;
    for (int ii = 0; ii < pair_wise_distances.size(); ii++) {
        for (int jj = 0; jj < pair_wise_distances[ii].size(); jj++) {
            sprintf(buffer, "%.2f ", pair_wise_distances[ii][jj]);
            cout << buffer;
        }
        cout << endl;
    } 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void linkage_group::return_order(vector<int>& out_order, 
                                 double & _lowerbound, 
                                 double & _upper_bound, 
                                 double & _cost_after_initialization, 
                                 vector<double> & _distances) const {
    out_order = current_order;
    _lowerbound = MST_lower_bound;
    _upper_bound = current_upper_bound;
    _cost_after_initialization = cost_after_initialization;
    _distances.clear();
    _distances.resize(number_of_bins-1);
    for (int ii = 1 ; ii < number_of_bins; ii++)
    {
        _distances[ii-1] = df->CM(pair_wise_distances[current_order[ii]][current_order[ii-1]] /number_of_individuals);
    }

};

//end


