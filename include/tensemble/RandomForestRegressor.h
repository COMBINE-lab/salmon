/* * * * *
 *  RandomForestRegressor.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_RandomForestRegressor_h
#define libTM_RandomForestRegressor_h
#include "BaseForest.h"

extern boost::mutex RF_mutex;
//global function for parallel build random forest regressor
void RFG_build_trees_range(REAL** X,REAL* y,uint n_samples,BaseForest* forest,uint n_trees_beg,uint n_trees_end,uint* oob_count,REAL* oob_prediction,REAL* importances);

class RandomForestRegressor:public BaseForest {
    //A random forest regressor
    
public:
    RandomForestRegressor():BaseForest(){};
    
    RandomForestRegressor(uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,bool bootstrap,bool oob,bool compute_importance,uint random_seed,uint n_jobs,bool verbose);
    
    ~RandomForestRegressor();
    
    int build(REAL** X,REAL* y,uint n_samples);
    
    void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures);
    
    int save_model(const char* filename);
    
};
RandomForestRegressor::RandomForestRegressor(uint n_trees,\
                                             uint n_features,\
                                             uint max_depth,\
                                             uint min_sample_leaf,\
                                             REAL max_features_ratio,\
                                             bool bootstrap,\
                                             bool oob,\
                                             bool compute_importance,\
                                             uint random_seed,\
                                             uint n_jobs,\
                                             bool verbose)\
:BaseForest(CRITERION_MSE,n_trees,n_features,max_depth,min_sample_leaf, max_features_ratio,FIND_BEST, bootstrap,oob,compute_importance,random_seed,n_jobs,verbose){    
}
RandomForestRegressor::~RandomForestRegressor(){
    
}

int RandomForestRegressor::build(REAL** X,REAL* y,uint n_samples){
    
    //set n_classes=1 for regression
    this->n_classes=1;
    
    time_t beg,end;
    beg=time(NULL);
    uint i,j,**oob_count,n_trees_beg,n_trees_end;
    REAL** oob_prediction,**importance_tmp;
    oob_count=new uint* [n_jobs];
    oob_prediction=new REAL* [n_jobs];
    importance_tmp=new REAL*[n_jobs];
    std::pair<uint,uint> n_trees_range;
    if (oob) {
        for (j=0; j<n_jobs; j++) {
            oob_count[j]=new uint[n_samples];
            oob_prediction[j]=new REAL[n_samples];
        }
    }
    if (compute_importance) {
        this->importances=new REAL[n_features];
        for (j=0; j<n_jobs; j++) {
            importance_tmp[j]=new REAL[n_features];
        }
    }
    boost::thread** thread=new boost::thread*[n_jobs];
    for (j=0; j<n_jobs; j++) {
        n_trees_beg=(n_trees/n_jobs)*j;
        n_trees_end=(n_trees/n_jobs)*(j+1);
        if (j==n_jobs-1) {
            //set the rest trees to the last job.
            n_trees_end=n_trees;
        }
        n_trees_range.first=n_trees_beg;
        n_trees_range.second=n_trees_end;
        thread[j]=new boost::thread(bind(RF_build_trees_range, X, y, n_samples, this, n_trees_range, oob_count[j], oob_prediction[j], importance_tmp[j]));
    }
    for (j=0; j<n_jobs; j++) {
        thread[j]->join();
        delete thread[j];
    }
    if (verbose) {
        fprintf(stderr, "\n");
    }
    if (oob) {
        bool warn_flag=false;
        for (i=0; i<n_samples; i++) {
            for (j=1; j<n_jobs; j++) {
                oob_count[0][i]+=oob_count[j][i];
                oob_prediction[0][i]+=oob_prediction[j][i];
            }
            if (oob_count[0][i]==0) {
                warn_flag==true?warn_flag=true:\
                (fprintf(stderr, "WARN: Some inputs do not have OOB scores.This probably means too few trees were used to compute any reliable oob estimates.\n"),warn_flag=true);
                oob_count[0][i]=1;
            }
            oob_prediction[0][i]/=oob_count[0][i];
        }
        REAL oob_R2=R2(oob_prediction[0], y, n_samples);
        oob_scores=rmse(oob_prediction[0], y, n_samples);
        if (1) {
            fprintf(stderr, "Out-of-bag score(RMSE,Correlation Coefficient)=(%lf,%lf)\n",oob_scores,oob_R2);
        }
        for (i=0; i<n_jobs; i++) {
            delete []oob_count[i];
            delete []oob_prediction[i];
        }
        delete []oob_count;
        delete []oob_prediction;
    }
    if (compute_importance) {
        for (i=0; i<n_features; i++) {
            importances[i]=0;
            for (j=0; j<n_jobs; j++) {
                importances[i]+=importance_tmp[j][i];
            }
            importances[i]/=n_trees;
        }
        for (i=0; i<n_jobs; i++) {
            delete []importance_tmp[i];
        }
        delete []importance_tmp;
    }
    end=time(NULL);
    fprintf(stderr, "|Random Forest training done. | Using time: %.0lf secs|\n",difftime(end, beg));
    return ENSEMBLE_SUCCESS;
}
void RandomForestRegressor::predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures){
    uint i,j;
    REAL* single_tree_pred=new REAL[nSamples];
#ifdef DEBUG
    assert(nFeatures==this->n_features);
#endif
    for (i=0; i<nSamples; i++) {
        pred[i]=0;
    }
    for (i=0; i<n_trees; i++) {
#ifdef DEBUG
        assert(tree[i]!=NULL);
#endif
        tree[i]->predict(X, single_tree_pred, nSamples, nFeatures);
        for (j=0; j<nSamples; j++) {
            pred[j]+=single_tree_pred[j];
        }
    }
    
    //average prediction
    for (i=0; i<nSamples; i++) {
        pred[i]/=n_trees;
    }
    delete [] single_tree_pred;
}


void RFG_build_trees_range(REAL** X,REAL* y,uint n_samples,BaseForest* forest,uint n_trees_beg,uint n_trees_end,uint* oob_count,REAL* oob_prediction,REAL* importances){
    Tree** tree=forest->tree;
    uint n_features=forest->n_features;
    uint max_features=forest->max_features;
    uint max_depth=forest->max_depth;
    uint min_sample_leaf=forest->min_sample_leaf;
    uint find_split_algorithm=forest->find_split_algorithm;
    uint random_seed=forest->random_seed;
    bool oob=forest->oob;
    bool compute_importance=forest->compute_importance;
    bool bootstrap=forest->bootstrap;
    
    uint i,j;
    REAL *oob_prediction_tmp,*importance_tmp;
    REAL** sub_X;
    REAL* sub_y;
    bool* mask;//for oob prediction
    if (oob) {
        oob_prediction_tmp=new REAL[n_samples];
        for (i=0; i<n_samples; i++) {
            oob_prediction[i]=0;
            oob_count[i]=0;
        }
        mask=new bool[n_samples];
    }
    if (compute_importance) {
        importance_tmp=new REAL[n_features];
        for (j=0; j<n_features; j++) {
            importances[j]=0;
        }
    }
    if (bootstrap) {
        sub_X=new REAL* [n_samples];
        sub_y=new REAL[n_samples];
    }
    for (i=n_trees_beg; i<n_trees_end; i++) {
        if (oob) {
            for (j=0; j<n_samples; j++) {
                mask[j]=false;
            }
        }
        if (bootstrap) {
            for (j=0; j<n_samples; j++) {
                //                uint idx=rand()%n_samples;
                uint idx=randomMT()%n_samples;
                sub_X[j]=X[idx];
                sub_y[j]=y[idx];
                if (oob) {
                    mask[idx]=true;
                }
            }
        }else {
            sub_X=X;
            sub_y=y;
        }
        if (1) {
            fprintf(stderr, "Random forest: building tree %d\n",i);
        }
        //build tree
        tree[i]=new TreeRegressor(n_features,\
                                  max_features,\
                                  min_sample_leaf,\
                                  max_depth,\
                                  find_split_algorithm,\
                                  random_seed,\
                                  1);
        tree[i]->build(sub_X, sub_y, n_samples);
//        tree[i]->resize();
        if (oob) {
            tree[i]->predict(X, oob_prediction_tmp, mask, n_samples, n_features);
            for (j=0; j<n_samples; j++) {
                if (mask[j]) {
                    continue;
                }
                oob_count[j]++;
                oob_prediction[j]+=oob_prediction_tmp[j];
            }
        }
        if (compute_importance) {
            tree[i]->compute_importance(importance_tmp);
            for (j=0; j<n_features; j++) {
                importances[j]+=importance_tmp[j];
            }
            
        }
    }
    if (oob) {
        delete []mask;
        delete []oob_prediction_tmp;
    }
    if (compute_importance) {
        delete []importance_tmp;
    }
    if (bootstrap) {
        delete []sub_X;
        delete []sub_y;
    }
}

int RandomForestRegressor::save_model(const char *filename){
    FILE* fp=fopen(filename,"w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for save RandomForestRegressor model.Please check your file path is correct.\n",filename);
        return false;
    }
    fprintf(fp, "RandomForestRegressor\n");
    BaseForest::save_model(fp);
    fclose(fp);
    return true;
}
#endif
