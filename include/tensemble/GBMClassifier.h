/* * * * *
 *  GBMClassifier.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_GBMClassifier_h
#define libTM_GBMClassifier_h
#include "BaseGBM.h"
#include <map>

class GBMClassifier:public BaseGBM {
    //gradient boosting regressor
    
public:
    GBMClassifier(){};
    
    GBMClassifier(int loss_function,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,REAL subsample,REAL learning_rate,bool oob,bool compute_importance,uint random_seed,uint n_jobs,int verbose);
    ~GBMClassifier();
    
    int build(REAL **X, REAL *original_y, uint n_samples,\
              REAL** val_X=NULL,REAL* val_y=NULL,uint n_val_samples=0);
    
    void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures,uint k_trees=0);
    
    void predict_prob(REAL** X,REAL* pred_prob,uint nSamples,uint nFeatures,uint k_trees=0);
    
    void predict_stage(REAL** X,REAL* pred,uint nSamples,uint nFeatures,uint k);
    
    void score2label(REAL* score,REAL* label,uint nSamples);
    
    void score2prob(REAL* score,REAL* prob,uint nSamples);
    
    int save_model(const char* filename);
    
    void predict_verbose(REAL** X,REAL* y,REAL* pred,uint nSamples,uint nFeatures,\
                                 uint k_trees){};
    
private:
    void convert_to_unique(REAL* original_y,REAL* y,uint nSamples);
    
};

GBMClassifier::GBMClassifier(int loss_function,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,REAL subsample,REAL learning_rate,bool oob,bool compute_importance,uint random_seed,uint n_jobs,int verbose)\
:BaseGBM(loss_function,n_trees,n_features,max_depth,min_sample_leaf,max_features_ratio,subsample,learning_rate,oob,compute_importance,random_seed,n_jobs,verbose){
#ifdef DEBUG
    assert(loss_function==BINOMIAL_DEVIANCE || loss_function==MULTINOMIAL_DEVIANCE);
#endif
}

GBMClassifier::~GBMClassifier(){
}

void GBMClassifier::convert_to_unique(REAL *original_y, REAL *y,uint nSamples){
    std::map<int,uint> key_values;
    std::map<int,uint>::const_iterator iter;
    vector<int> classes_map;
    n_classes=0;
    for (uint i=0; i<nSamples; i++) {
        int origin_class=(int)original_y[i];
        if ((iter=key_values.find(origin_class))==key_values.end()) {
            key_values.insert(make_pair(origin_class, n_classes));
            y[i]=n_classes;
            classes_map.push_back(origin_class);
            n_classes++;
        }else {
            y[i]=iter->second;
        }
    }
    original_classes=new int[n_classes];
    for (uint i=0; i<n_classes; i++) {
        original_classes[i]=classes_map[i];
    }
}

int GBMClassifier::build(REAL **X, REAL *original_y, uint n_samples,REAL** val_X,REAL* val_y,uint n_val_samples){
    /* BUILD GRADIENT BOOSTING MODEL FOR CLASSIFICATION
     * For classification, labels must correspond to classes 0, 1, ..., n_classes_-1
     * For regression and two-classes classification,just need fit one tree in each GBM round
     * For multi-classes classification,need fit n_classes trees in each GBM round
     */
#ifdef DEBUG
    assert(loss_function==BINOMIAL_DEVIANCE || loss_function==MULTINOMIAL_DEVIANCE);
#endif
    uint i,j,k,idx,n_subsamples,*sample_index,count,offset;
    REAL** sub_X,*sub_y,*y,*zero_one_y;
    REAL* y_pred,*y_pred_label,*val_pred,*val_pred_stage,*val_pred_label,*residual,train_score,oob_score;
    bool *mask;
    //convert y to unique classes
    y=new REAL[n_samples];
    y_pred_label=new REAL[n_samples];
    convert_to_unique(original_y, y, n_samples);
    
    if (n_classes==1) {
        fprintf(stderr, "Error: The training file only have one label,at least two-labels are required for classification.Please check your file.\n");
        delete []y;
        return ENSEMBLE_FAIL;
    }
    
    //initialize loss function
    if (n_classes==2) {
        loss_function=BINOMIAL_DEVIANCE;
        loss=new BinomialDeviance;
    }else {
        loss_function=MULTINOMIAL_DEVIANCE;
        loss=new MultinomialDeviance;
    }
    //set n_classes=1 for two-classes problems(just need build one tree each round)
    if (loss_function==BINOMIAL_DEVIANCE) {
        n_classes=1;
    }
    //allocate memory for prior and trees
    prior_pred=new REAL[n_classes];
    this->tree=new Tree*[n_trees*n_classes];
    
    n_subsamples=n_samples;
    mask=new bool[n_samples];
    y_pred=new REAL[n_samples*n_classes];
    residual=new REAL[n_samples*n_classes];
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

    //initialize prediction and save prior
    loss->get_init_estimate(y, y_pred, prior_pred, n_samples, this->n_classes);
    
    //if have validation set,initialize prediction for validation set
    if (val_X && val_y) {
        val_pred=new REAL[n_classes*n_val_samples];
        val_pred_stage=new REAL[n_classes*n_val_samples];
        val_pred_label=new REAL[n_val_samples];
        for (j=0; j<n_classes; j++) {
            offset=j*n_val_samples;
            for (i=0; i<n_val_samples; i++) {
                val_pred[offset+i]=prior_pred[j];
            }
        }
    }
    if (loss_function==MULTINOMIAL_DEVIANCE) {
        zero_one_y=new REAL[n_samples];
    }
    time_t beg,end;
    beg=time(NULL);
    //main iteration
    for (i=0; i<n_trees; i++) {
        //sub-sampling
        memset(mask, false, n_samples*sizeof(bool));
        if (subsample<1.0) {
            count=n_samples;
            for (j=0; j<n_subsamples; j++) {
//                k=rand()%count;
                k=randomMT()%count;
                idx=sample_index[k];
                sub_X[j]=X[idx];
                mask[idx]=true;
                swap(sample_index[k], sample_index[--count]);
            }
        }else {
            sub_X=X;
            for (j=0; j<n_samples; j++) {
                mask[j]=true;
            }
        }
        //build k trees each round
        for(k=0;k<n_classes;k++){
            REAL* kth_residual=residual+k*n_samples;
            //convert to zero-one value for Multinomial Deviance
            if (loss_function==MULTINOMIAL_DEVIANCE) {
                for (j=0; j<n_samples; j++) {
                    if (y[j]==k) {
                        zero_one_y[j]=1;
                    }else {
                        zero_one_y[j]=0;
                    }
                }
            }else {
                zero_one_y=y;
            }
            //compute working residual for k-th tree
            loss->compute_residual(kth_residual, zero_one_y, y_pred, n_samples, n_classes, k);
            if (subsample<1.0) {
                count=0;
                for (j=n_samples-1; j>=n_samples-n_subsamples; j--) {
                    idx=sample_index[j];
                    sub_y[count++]=kth_residual[idx];
                }
#ifdef DEBUG
                assert(count==n_subsamples);
#endif
            }else {
                sub_y=kth_residual;
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
            loss->update(t, X, mask, zero_one_y, y_pred+k*n_samples, kth_residual, n_samples, n_features, learning_rate, n_classes,k);
            //add tree to ensemble
            tree[i*n_classes+k]=t;
        }         
        if (val_X && val_y) {
            predict_stage(val_X, val_pred_stage, n_val_samples, n_features, i);
            for (k=0; k<n_classes; k++) {
                offset=k*n_val_samples;
                for (j=0; j<n_val_samples; j++) {
                    val_pred[offset+j]+=this->learning_rate*val_pred_stage[offset+j];
                }
            }
        }
        if ((verbose && (i+1)%verbose==0) || i+1==n_trees) {
            //calculate train score and oob score;
//            train_score=loss->loss(y, y_pred, mask, true, n_samples, n_classes);
            score2label(y_pred, y_pred_label, n_samples);
            if (oob && false) {
                oob_score=loss->loss(y, y_pred, mask, false, n_samples, n_classes);
                oob_score=sqrt(oob_score);
            }
            if (verbose<=0) {
                fprintf(stderr, "GBM train done.");
            }else {
                fprintf(stderr, "build tree %u of %u,",i+1,n_trees);
            }
            fprintf(stderr, "train Acc=%lf",Accuracy(original_y, y_pred_label, n_samples));
            if (oob && false) {
                fprintf(stderr, ",oob score(deviance)=%lf",oob_score);
            }
            if (val_y) {
                score2label(val_pred, val_pred_label, n_val_samples);
                fprintf(stderr, ",validation Acc=%f",Accuracy(val_y, val_pred_label, n_val_samples));
            }
            fprintf(stderr, ".\n");
        }
    }
    delete []mask;
    delete []y;
    delete []y_pred;
    delete []y_pred_label;
    delete []residual;
    if (subsample<1.0) {
        delete []sample_index;
        delete []sub_X;
        delete []sub_y;
    }
    if (val_y) {
        delete []val_pred_stage;
        delete []val_pred;
        delete []val_pred_label;
    }
    if (loss_function==MULTINOMIAL_DEVIANCE) {
        delete []zero_one_y;
    }
    end=time(NULL);
    fprintf(stderr, "|Gradient Boosting Classifier training done. | Using time: %.0lf secs|\n",difftime(end, beg));
    return ENSEMBLE_SUCCESS;
}

void GBMClassifier::predict_stage(REAL **X, REAL *pred, uint nSamples, uint nFeatures, uint k){
    /* predict k-th round */
    for (uint i=0; i<n_classes; i++) {
        tree[k*n_classes+i]->predict(X, pred+i*nSamples, nSamples, nFeatures);
    }
}

void GBMClassifier::score2label(REAL *score, REAL *label,uint nSamples){
    uint label_idx;
    REAL max_score;
    if (loss_function==BINOMIAL_DEVIANCE) {
        for (uint i=0; i<nSamples; i++) {
            if (1.0/(1.0+exp(-score[i]))>=0.5) {
                label[i]=original_classes[1];
            }else {
                label[i]=original_classes[0];
            }
        }
    }else{
        for (uint i=0; i<nSamples; i++) {
            max_score=-HUGE_VAL;
            for (uint j=0; j<n_classes; j++) {
                if (score[j*nSamples+i]>max_score) {
                    max_score=score[j*nSamples+i];
                    label_idx=j;
                }
            }
            label[i]=original_classes[label_idx];
        }
    }
}

void GBMClassifier::predict(REAL **X, REAL *pred, uint nSamples, uint nFeatures,uint k_trees){
    if (k_trees==0 || k_trees>this->n_trees) {
        k_trees=n_trees;
    }
    uint i,j,k,offset;
    REAL* score=new REAL[nSamples*n_classes];
    REAL* score_stage=new REAL[nSamples];
    for (i=0; i<n_classes; i++) {
        offset=i*nSamples;
        for (j=0; j<nSamples; j++) {
            score[offset+j]=prior_pred[i];
        }
    }
    for (i=0; i<k_trees; i++) {
        for (k=0; k<n_classes; k++) {
            offset=k*nSamples;
            tree[i*n_classes+k]->predict(X, score_stage, nSamples, nFeatures);
            for (j=0; j<nSamples; j++) {
                score[offset+j]+=learning_rate*score_stage[j];
            }
        }
    }
    score2label(score, pred, nSamples);
    delete []score;
    delete []score_stage;
}
void GBMClassifier::score2prob(REAL *score, REAL *prob, uint nSamples){
    uint i,j;
    REAL sum_exp;
    if (loss_function==BINOMIAL_DEVIANCE) {
        for (i=0; i<nSamples; i++) {
            prob[i]=1.0/(1.0+exp(-score[i]));
        }
    }else{
        for (i=0; i<nSamples; i++) {
            sum_exp=0.0;
            for (j=0; j<n_classes; j++) {
                sum_exp+=exp(score[j*nSamples+i]);
            }
            for (j=0; j<n_classes; j++) {
                if (sum_exp<=EPS) {
                    prob[j*nSamples+i]=0.0;
                }else{
                    prob[j*nSamples+i]=exp(score[j*nSamples+i])/sum_exp;
                }
            }
        }
    }
}
void GBMClassifier::predict_prob(REAL **X, REAL *pred_prob, uint nSamples, uint nFeatures, uint k_trees){
    if (k_trees==0 || k_trees>this->n_trees) {
        k_trees=n_trees;
    }
    uint i,j,k,offset;
    REAL* score=new REAL[nSamples*n_classes];
    REAL* score_stage=new REAL[nSamples];
    for (i=0; i<n_classes; i++) {
        offset=i*nSamples;
        for (j=0; j<nSamples; j++) {
            score[offset+j]=prior_pred[i];
        }
    }
    for (i=0; i<k_trees; i++) {
        for (k=0; k<n_classes; k++) {
            offset=k*nSamples;
            tree[i*n_classes+k]->predict(X, score_stage, nSamples, nFeatures);
            for (j=0; j<nSamples; j++) {
                score[offset+j]+=learning_rate*score_stage[j];
            }
        }
    }
    score2prob(score, pred_prob, nSamples);
    delete []score;
    delete []score_stage;

}

int GBMClassifier::save_model(const char *filename){
    FILE* fp=fopen(filename,"w");
    if (!fp) {
        fprintf(stderr, "Cannot open file %s for save GBMClassifier model.Please check your file path is correct.\n",filename);
        return false;
    }
    fprintf(fp, "GBMClassifier\n");
    BaseGBM::save_model(fp);
    fclose(fp);
    return true;

}
#endif
