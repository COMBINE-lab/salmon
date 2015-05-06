/* * * * *
 *  GBMRegressor.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_GBMRegressor_h
#define libTM_GBMRegressor_h
#include "BaseGBM.h"


class GBMRegressor:public BaseGBM {
    //gradient boosting regressor
    
public:
    GBMRegressor(){};
    
    GBMRegressor(int loss_function,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,REAL subsample,REAL learning_rate,bool oob,bool compute_importance,uint random_seed,uint n_jobs,int verbose);
    
    ~GBMRegressor();
    
    int build(REAL **X, REAL *y, uint n_samples,\
              REAL** val_X=NULL,REAL* val_y=NULL,uint n_val_samples=0);
    
    void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures,uint k_trees=0);
    
    //may useful for validation
    void predict_verbose(REAL** X,REAL* y,REAL* pred,uint nSamples,uint nFeatures,\
                         uint k);
    int save_model(const char* filename);
};

GBMRegressor::GBMRegressor(int loss_function,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,REAL subsample,REAL learning_rate,bool oob,bool compute_importance,uint random_seed,uint n_jobs,int verbose)\
:BaseGBM(loss_function,n_trees,n_features,max_depth,min_sample_leaf,max_features_ratio,subsample,learning_rate,oob,compute_importance,random_seed,n_jobs,verbose){
}

GBMRegressor::~GBMRegressor(){
    
}

int GBMRegressor::build(REAL **X, REAL *y, uint n_samples,REAL** val_X,REAL* val_y,uint n_val_samples){
    /* BUILD GRADIENT BOOSTING MODEL FOR REGRESSION
     * For classification, labels must correspond to classes 0, 1, ..., n_classes_-1
     * For regression and two-classes classification,just need fit one tree in each GBM round
     * For multi-classes classification,need fit n_classes trees in each GBM round
     */
    uint i,j,k,idx,n_subsamples,*sample_index,count;
    REAL** sub_X,*sub_y;
    REAL* y_pred,*val_pred,*val_pred_stage,*residual,train_score,oob_score,valid_score;
    bool *mask;
    val_pred=NULL;
#ifdef DEBUG
    assert(n_trees>=1);
#endif
    //initialize loss function and set n_classes=1.
    n_classes=1;
    if (loss_function==SQUARE_LOSS) {
        loss=new SquareLoss();
    }else {
        fprintf(stderr, "Error: Unknow loss function\n");
        return ENSEMBLE_FAIL;
    }
    //allocate memory for prior and trees
    prior_pred=new REAL[n_classes];
    this->tree=new Tree*[n_trees];
    
    n_subsamples=n_samples;
    mask=new bool[n_samples];
    y_pred=new REAL[n_samples*n_classes];
    if (subsample<1.0) {
        sample_index=new uint[n_samples];
        n_subsamples=(uint)floor((REAL)n_samples*subsample);
        sub_X=new REAL*[n_subsamples];
        sub_y=new REAL[n_subsamples];
        for (i=0; i<n_samples; i++) {
            sample_index[i]=i;
        }
#ifdef DEBUG
        assert(n_subsamples<n_samples);
#endif
    }
    residual=new REAL[n_samples];
    //initialize prediction and save prior
    loss->get_init_estimate(y, y_pred, prior_pred, n_samples, this->n_classes);
    if (val_X && val_y) {
        val_pred=new REAL[n_val_samples];
        val_pred_stage=new REAL[n_val_samples];
        for (i=0; i<n_val_samples; i++) {
            val_pred[i]=*prior_pred;
        }
    }
    time_t beg,end;
    beg=time(NULL);
    //main iteration
    for (i=0; i<n_trees; i++) {
        //compute working residual
        loss->compute_residual(residual, y, y_pred, n_samples, n_classes, 0);
        //sub-sampling
        memset(mask, false, n_samples*sizeof(bool));
        if (subsample<1.0) {
            count=n_samples;
            for (j=0; j<n_subsamples; j++) {
//                k=rand()%count;
                k=randomMT()%count;
                idx=sample_index[k];
                sub_X[j]=X[idx];
                sub_y[j]=residual[idx];
                mask[idx]=true;
                swap(sample_index[k], sample_index[--count]);
            }
        }else {
            sub_X=X;
            sub_y=residual;
            for (j=0; j<n_samples; j++) {
                mask[j]=true;
            }
        }
        uint tree_random_seed=randomMT();
        TreeRegressor *t=new TreeRegressor(n_features,\
                                           max_features,\
                                           min_sample_leaf,\
                                           max_depth,\
                                           FIND_BEST,\
                                           tree_random_seed,\
                                           n_jobs);
        //build tree
        t->build(sub_X, sub_y, n_subsamples);
        //update
        loss->update(t, X, mask, y, y_pred, residual, n_samples, n_features, learning_rate, n_classes,0);
        //add tree to ensemble
        tree[i]=t;
        if (val_X && val_y) {
            t->predict(val_X, val_pred_stage, n_val_samples, n_features);
            for (j=0; j<n_val_samples; j++) {
                val_pred[j]+=learning_rate*val_pred_stage[j];
            }
        }
        if ( (verbose && (i+1)%verbose==0) || i+1==n_trees) {      
            //calculate train score and oob score;
            train_score=loss->loss(y, y_pred, mask, true, n_samples, n_classes);
            train_score=sqrt(train_score);
            if (oob && false) {
                oob_score=loss->loss(y, y_pred, mask, false, n_samples, n_classes);
                oob_score=sqrt(oob_score);
            }
            if (verbose<=0) {
                fprintf(stderr, "GBM train done.");
            }else {
                fprintf(stderr, "build tree %u of %u,",i+1,n_trees);
            }
            fprintf(stderr, "train RMSE=%f",train_score);
            if (oob && false) {
                fprintf(stderr, ",oob score(deviance)=%f",oob_score);
            }
            if (val_y) {
                valid_score=rmse(val_pred, val_y, n_val_samples);
                fprintf(stderr, ",validation RMSE=%f",valid_score);
            }
            fprintf(stderr, ".\n");

        }
    }
    delete []mask;
    delete []y_pred;
    delete []residual;
    if (subsample<1.0) {
        delete []sample_index;
        delete []sub_X;
        delete []sub_y;
    }
    if (val_y) {
        delete []val_pred_stage;
        delete []val_pred;
    }
    end=time(NULL);
    fprintf(stderr, "|Gradient Boosting Regressor training done. | Using time: %.0lf secs|\n",difftime(end, beg));
    return ENSEMBLE_SUCCESS;
}

void GBMRegressor::predict(REAL **X, REAL *pred, uint nSamples, uint nFeatures,uint k_trees){
#ifdef DEBUG
    assert(k_trees<=n_trees);
#endif
    if (k_trees==0) {
        k_trees=n_trees;
    }
    uint i,j;
    REAL* pred_stage=new REAL[nSamples];
    for (i=0; i<nSamples; i++) {
        pred[i]=*prior_pred;
    }
    for (i=0; i<k_trees; i++) {
        tree[i]->predict(X, pred_stage, nSamples, nFeatures);
        for (j=0; j<nSamples; j++) {
            pred[j]+=learning_rate*pred_stage[j];
        }
    }
    delete []pred_stage;
}

void GBMRegressor::predict_verbose(REAL **X, REAL *y, REAL *pred, uint nSamples, uint nFeatures,uint k_trees){
#ifdef DEBUG
    assert(k_trees<=n_trees);
#endif
    if (k_trees==0) {
        k_trees=n_trees;
    }
    uint i,j;
    REAL* pred_stage=new REAL[nSamples];
    for (i=0; i<nSamples; i++) {
        pred[i]=*prior_pred;
    }
    for (i=0; i<k_trees; i++) {
        tree[i]->predict(X, pred_stage, nSamples, nFeatures);
        for (j=0; j<nSamples; j++) {
            pred[j]+=learning_rate*pred_stage[j];
        }
        REAL MSE=mse(pred, y, nSamples);
        REAL Rscore=R2(pred,y,nSamples);
        fprintf(stderr,"GBMRegressor prediction: tree %u,MSE=%f,R^2=%f\n",i+1,MSE,Rscore);
    }
    delete []pred_stage;
}

int GBMRegressor::save_model(const char *filename){
    FILE* fp=fopen(filename,"w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for save GBMRegressor model.Please check your file path is correct.\n",filename);
        return false;
    }
    fprintf(fp, "GBMRegressor\n");
    BaseGBM::save_model(fp);
    fclose(fp);
    return true;
}
#endif
