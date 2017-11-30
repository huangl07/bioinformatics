/*
 *  single_mapping_population_raw_data.cpp
 *  ApproxMap
 *
 *  Created by yonghui on 4/7/07.linkage_group_DH
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "genetic_map_DH.h"


genetic_map_DH::~genetic_map_DH(){
};

void genetic_map_DH::calculate_pair_wise_distance()
{
    pair_wise_distances.resize(number_of_loci);
    for (int ii = 0 ; ii < number_of_loci; ii++)
    {
        pair_wise_distances[ii].resize(number_of_loci, 0.0);
    }

    cout << "start calculating pair-wise distance" << time(NULL) << endl;
    for (int ii = 0; ii < number_of_loci; ii++)
    {
        for (int jj = ii ; jj < number_of_loci; jj++)
        {
            double distance_ii_jj = 0;
            double none_missing = 0;
            for (int kk = 0 ; kk < number_of_individual; kk++)
            {
                
                if ((raw_mapping_data[ii][kk] != '-') and 
                    (raw_mapping_data[jj][kk] != '-')) {
                    none_missing = none_missing + 1.0;
                    if (raw_mapping_data[ii][kk] != raw_mapping_data[jj][kk]) {
                        distance_ii_jj = distance_ii_jj + 1.0;
                    }
                }
            }
            if (none_missing < 0.5 * number_of_individual) { 
                cout << "caution: too many missing for pair:(" 
                     << marker_names[ii] << "," << marker_names[jj] << ")" << endl;
            }
            
            if (none_missing < 0.25 * number_of_individual) { // almost everything is missing, adjust the estimate 
                distance_ii_jj = 0.5 * number_of_individual;
                none_missing = number_of_individual;
            }
            pair_wise_distances[ii][jj] = (distance_ii_jj / none_missing) * number_of_individual;
            pair_wise_distances[jj][ii] = pair_wise_distances[ii][jj];
        }
    }    
    cout << "finished calculating pair-wise distance:" << time(NULL) << endl;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_DH* genetic_map_DH::construct_linkage_group(int group_id)
{
    int _number_of_bins = (linkage_group_bins[group_id]).size();
    int _number_of_individuals = number_of_individual;
    
    /*Store the probability for each allele to be A*/
    vector<vector< float > > _raw_data ;
    
    vector<pair<int, int> >  _missing_data;
    
    vector<int> _current_order;
    
    _raw_data.resize(_number_of_bins);
    for (int ii = 0 ; ii < _number_of_bins; ii++)
    {
        _raw_data[ii].resize(_number_of_individuals);
        for (int jj = 0; jj < _number_of_individuals; jj++)
        {
            if (raw_mapping_data[linkage_group_bins[group_id][ii][0]][jj] == 'A') 
            {
                /*If an allele is A, then its probability being A is A*/
                _raw_data[ii][jj] = 1.0;
            }
            else if (raw_mapping_data[linkage_group_bins[group_id][ii][0]][jj] == 'B') 
            {
                /*If an allele is B, then its probability of being A is 0*/
                _raw_data[ii][jj] = 0.0;
            }
            else 
            {
                /*If an allele is missing, then, assign probability 0.5 for it to be A*/
                _raw_data[ii][jj] = 0.5;
                _missing_data.push_back(make_pair(ii,jj)); /*ii is the id for the marker, and jj is the id for the individual*/
            }
        }
    }
    for (int ii = 0 ; ii < _number_of_bins; ii ++)
    {
        _current_order.push_back(ii);
    }
    vector<int> bin_sizes;
    for (int ii = 0; ii < _number_of_bins; ii++) {
        bin_sizes.push_back(linkage_group_bins[group_id][ii].size());
    }
    linkage_group_DH * to_be_returned = new linkage_group_DH(_number_of_bins, 
                                                             _number_of_individuals,
                                                             detect_bad_data, 
                                                             objective_function,
                                                             df_,
                                                             _raw_data, 
                                                             _current_order, 
                                                             _missing_data,
                                                             bin_sizes);
    return to_be_returned;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

linkage_group_DH* genetic_map_DH::construct_linkage_group_whole_map()
{
    int _number_of_bins = number_of_loci;
    int _number_of_individuals = number_of_individual;
    
    /*Store the probability for each allele to be A*/
    vector<vector< float > > _raw_data ;
    
    vector<pair<int, int> >  _missing_data;
    
    vector<int> _current_order;
    
    _raw_data.resize(_number_of_bins);
    for (int ii = 0 ; ii < _number_of_bins; ii++)
    {
        _raw_data[ii].resize(_number_of_individuals);
        for (int jj = 0; jj < _number_of_individuals; jj++)
        {
            if (raw_mapping_data[ii][jj] == 'A') 
            {
                /*If an allele is A, then its probability being A is A*/
                _raw_data[ii][jj] = 1.0;
            }
            else if (raw_mapping_data[ii][jj] == 'B') 
            {
                /*If an allele is B, then its probability of being A is 0*/
                _raw_data[ii][jj] = 0.0;
            }
            else 
            {
                /*If an allele is missing, then, assign probability 0.5 for it to be A*/
                _raw_data[ii][jj] = 0.5;
                _missing_data.push_back(make_pair(ii,jj)); /*ii is the id for the marker, and jj is the id for the individual*/
            }
        }
    }
    for (int ii = 0 ; ii < _number_of_bins; ii ++)
    {
        _current_order.push_back(ii);
    }
    vector<int> bin_sizes;
    for (int ii = 0; ii < _number_of_bins; ii++) {
        bin_sizes.push_back(1);
    }
    linkage_group_DH* to_be_returned = new linkage_group_DH(_number_of_bins, 
                                                             _number_of_individuals,
                                                             false, // this is fixed to be false for whole map 
                                                             OBJF_COUNT,
                                                             df_,
                                                             _raw_data, 
                                                             _current_order, 
                                                             _missing_data,
                                                             bin_sizes);
    return to_be_returned;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genetic_map_DH::print_suspicious_data(){
    cout << endl;
    for (int ii = 0; ii < suspicious_data.size(); ii++) {
        cout << suspicious_data[ii].first;
        cout << '\t';
        cout << suspicious_data[ii].second;
        cout << endl;
    }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_DH::generate_map()
{
    pair_wise_distances.resize(number_of_loci);
    for (int ii = 0 ; ii < number_of_loci; ii++)
    {
        pair_wise_distances[ii].resize(number_of_loci, 0.0);
    }    
    /*
      if the total number of missing observations exceeds certain threshold, 
      we need to estimate the missing data before clustering
    */
    if ((total_number_of_missing_obs >= ESTIMATION_BEFORE_CLUSTERING * number_of_loci * number_of_individual) and 
        (estimation_before_clustering)) {
        linkage_group_DH * linkage_group_whole_map = construct_linkage_group_whole_map();
        linkage_group_whole_map->order_markers();
        const vector<vector<double> > & new_dist = linkage_group_whole_map-> get_pair_wise_distance();
        for (int ii = 0 ; ii < number_of_loci; ii++)
        {
            for (int jj = 0 ; jj < number_of_loci; jj++)
            {
                pair_wise_distances[ii][jj] = new_dist[ii][jj];
            }
        }
        // linkage_group_whole_map->dump_distance_matrix();
        delete linkage_group_whole_map;
    } else {
        cout << "calculating the pair-wise hamming distance" << endl;
        calculate_pair_wise_distance();
        cout << "finished calculating the pair-wise hamming distance" << endl;
    }
    // dump_distance_matrix();
    cluster();
    // dump_distance_matrix();
    cout << "found " << number_of_connected_components << " connected components" << endl;
    
    condense_markers_into_bins();
    // dump_distance_matrix();
    
    orders.resize(number_of_connected_components);
    upperbounds.resize(number_of_connected_components);
    lowerbounds.resize(number_of_connected_components);
    approx_bounds.resize(number_of_connected_components);
    distance_between_adjacent_pairs.resize(number_of_connected_components);
    
    for (int ii = 0 ; ii < number_of_connected_components; ii++)
    {
        linkage_group_DH * current_linkage_group = construct_linkage_group(ii);

        current_linkage_group->order_markers();
        current_linkage_group->return_order(orders[ii], 
                                            lowerbounds[ii], 
                                            upperbounds[ii], 
                                            approx_bounds[ii], 
                                            distance_between_adjacent_pairs[ii]);
        vector<pair<int, int> > bad_data_ii;                                    
        current_linkage_group->bad_genotypes(bad_data_ii);
        for (int jj = 0; jj < bad_data_ii.size(); jj++) {
            int bin_id = bad_data_ii[jj].first;
            int indi_id = bad_data_ii[jj].second;
            for (int kk = 0; kk < linkage_group_bins[ii][bin_id].size(); kk++) {
                string marker_name = marker_names[linkage_group_bins[ii][bin_id][kk]];
                string indi_name = individual_names[indi_id];
                suspicious_data.push_back(make_pair(marker_name, indi_name));
            }
        }
        current_linkage_group->dump();
        delete current_linkage_group;
        cout <<"finished the " << ii+1 << " linkage group" << endl; 
    }
    
    // Added by Yonghui on Oct 20, 2007 
    // The last step is to condense adjacent bins if they are too close to each other
    condense_bin();
    cout << "suspicious data detected by our algorithm" << endl;
    print_suspicious_data();
    cout << "double cross overs based on the current order" << endl;
    print_double_cross_overs();
    // dump the distance matrix
    // dump_connected_components_edges();
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void genetic_map_DH::print_double_cross_overs(){
    for (int ii = 0; ii < lg_bins_condensed.size(); ii++) {
        if (lg_bins_condensed[ii].size() < 3) {
            continue;
        }
        for (int jj = 0; jj < lg_bins_condensed[ii].size(); jj++){
            if (lg_bins_condensed[ii][jj].size() > 1) {
                continue;
            }
            int marker_id = lg_bins_condensed[ii][jj][0];           
            int pre_marker_id = -1;
            if (jj == 0){ // it is the first bin
                pre_marker_id = lg_bins_condensed[ii][1][0];
            } else {
                pre_marker_id = lg_bins_condensed[ii][jj - 1][0];
            }
            
            int next_marker_id = -1;
            if (jj == lg_bins_condensed[ii].size() - 1) { // it is the last bin
                next_marker_id = lg_bins_condensed[ii][lg_bins_condensed[ii].size() - 2][0];
            } else {
                next_marker_id = lg_bins_condensed[ii][jj + 1][0];
            }
            
            for (int kk = 0; kk < number_of_individual; kk++) {
                if (raw_mapping_data[marker_id][kk] == '-') { // ignore missing 
                    continue;
                }
                if ((raw_mapping_data[marker_id][kk] != raw_mapping_data[pre_marker_id][kk]) and 
                    (raw_mapping_data[marker_id][kk] != raw_mapping_data[next_marker_id][kk])) {
                    // this is a double cross-over
                    cout << marker_names[marker_id] << "," << individual_names[kk] << endl;
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

