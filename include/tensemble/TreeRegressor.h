/* * * * *
 *  TreeRegressor.h 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_TreeRegressor_h
#define libTM_TreeRegressor_h
#include "Tree.h"

class TreeRegressor:public Tree {
    //DecisionTreeRegressor
    
public:
    //constructor,regression set n_classes=1; 
    TreeRegressor(uint n_features,\
                  uint max_features,\
                  uint min_sample_leaf,\
                  uint max_depth,\
                  uint find_split_algorithm=FIND_BEST,\
                  uint random_seed=0,\
                  uint n_threads=1);
    void compute_importance(REAL* importance);
};

TreeRegressor::TreeRegressor(uint n_features,\
                             uint max_features,\
                             uint min_sample_leaf,\
                             uint max_depth,\
                             uint find_split_algorithm,\
                             uint random_seed,\
                             uint n_threads)\
:Tree(CRITERION_MSE,1,n_features,max_features,min_sample_leaf,max_depth,find_split_algorithm,random_seed,n_threads){
}
void TreeRegressor::compute_importance(REAL* importance){
    //node from 1....numNodes
    uint f;
    REAL normalize=0;
    for (uint i=0; i<n_features; i++) {
        importance[i]=0.0;
    }
    for (uint i=ROOT; i<=numNodes; i++) {
        // don't calculate in redundant nodes or leaf node
        if (nodes[i]->leaf || nodes[i]->nSamples<=0) {
            continue;
        }
        f=nodes[i]->feature_split;
        importance[f]+=\
        (nodes[i]->ini_error-nodes[i]->best_error)*(nodes[i]->ini_error-nodes[i]->best_error);
    }
    for (uint i=0; i<n_features; i++) {
        normalize+=importance[i];
    }
    if (normalize!=0) {
        for (uint i=0; i<n_features; i++) {
            importance[i]/=normalize;
        }
    }
}
#endif
