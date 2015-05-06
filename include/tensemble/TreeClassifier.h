/* * * * *
 *  TreeClassifier.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */


#ifndef libTM_TreeClassifier_h
#define libTM_TreeClassifier_h

#include "Tree.h"

class TreeClassifier:public Tree {
    /* DecisionTreeClassifier
     * label must be [0,1,...,n_classes-1] before calling TreeClassifier
     */
    
public:
    TreeClassifier(int split_criterion,\
                   uint n_classes,\
                   uint n_features,\
                   uint max_features,\
                   uint min_sample_leaf,\
                   uint max_depth,\
                   uint find_split_algorithm=FIND_BEST,\
                   uint random_seed=0,\
                   uint n_threads=1);
    
    void compute_importance(REAL* importance);
};

TreeClassifier::TreeClassifier(int split_criterion,\
                               uint n_classes,\
                               uint n_features,\
                               uint max_features,\
                               uint min_sample_leaf,\
                               uint max_depth,\
                               uint find_split_algorithm,\
                               uint random_seed,\
                               uint n_threads)\
:Tree(split_criterion,n_classes,n_features,max_features,min_sample_leaf,max_depth,find_split_algorithm,random_seed,n_threads){
}
void TreeClassifier::compute_importance(REAL* importance){
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
        nodes[i]->nSamples*(nodes[i]->ini_error-nodes[i]->best_error);
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
