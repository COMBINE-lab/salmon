/* * * * *
 *  RandomForestClassifier.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_RandomForestClassifier_h
#define libTM_RandomForestClassifier_h
#include "BaseForest.h"

//global function for parallel build random forest classifier
void RFC_build_trees_range(REAL** X,REAL* y,uint n_samples,BaseForest* forest,std::pair<uint, uint> n_trees_range,uint* oob_count,REAL* oob_prediction,REAL* importances);

class RandomForestClassifier:public BaseForest {
    //A random forest regressor
    
public:
    RandomForestClassifier():BaseForest(){};
    
    RandomForestClassifier(int split_criterion,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,bool bootstrap,bool oob,bool compute_importance,uint random_seed,uint n_jobs,bool verbose);
    
    ~RandomForestClassifier();
    
    int build(REAL** X,REAL* original_y,uint n_samples);
    
    void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures);
    
    void predict_prob(REAL** X,REAL* pred_prob,uint nSamples,uint nFeatures);
    
    void score2label(REAL* score,REAL* pred_label,uint nSamples);
    
    int save_model(const char* filename);
    
private:
    void convert_to_unique(REAL* original_y,REAL* y,uint nSamples);
};
RandomForestClassifier::RandomForestClassifier(int split_criterion,\
                                             uint n_trees,\
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
:BaseForest(split_criterion, n_trees,n_features,max_depth,min_sample_leaf, max_features_ratio,FIND_BEST, bootstrap,oob,compute_importance,random_seed,n_jobs,verbose){    
    
}
RandomForestClassifier::~RandomForestClassifier(){
    
}


int RandomForestClassifier::build(REAL** X,REAL* original_y,uint n_samples){
    time_t beg,end;
    beg=time(NULL);
    uint i,j,k,**oob_count,n_trees_beg,n_trees_end;
    REAL** oob_prediction,**importance_jobs,*oob_label_prediction;
    REAL* y;
    
    //convert y to unique classes
    y=new REAL[n_samples];
    convert_to_unique(original_y, y, n_samples);
    if (n_classes==1) {
        fprintf(stderr, "The training file only have one label,at least two-labels are required for classification.Please check your file.\n");
        delete []y;
        return ENSEMBLE_FAIL;
    }
    
    std::pair<uint, uint> n_trees_range;
    oob_count=new uint* [n_jobs];
    oob_prediction=new REAL* [n_jobs];
    importance_jobs=new REAL*[n_jobs];
    if (oob) {
        oob_label_prediction=new REAL[n_samples];
        for (j=0; j<n_jobs; j++) {
            oob_count[j]=new uint[n_samples];
            oob_prediction[j]=new REAL[n_samples*n_classes];
        }
    }
    if (compute_importance) {
        this->importances=new REAL[n_features];
        for (j=0; j<n_jobs; j++) {
            importance_jobs[j]=new REAL[n_features];
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
        thread[j]=new boost::thread(bind(RF_build_trees_range, X, y, n_samples, this, n_trees_range, oob_count[j], oob_prediction[j], importance_jobs[j]));
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
                for (k=0; k<n_classes; k++) {
                    oob_prediction[0][i*n_classes+k]+=oob_prediction[j][i*n_classes+k];
                }
            }
            if (oob_count[0][i]==0) {
                warn_flag==true?warn_flag=true:\
                (fprintf(stderr, "WARN: Some inputs do not have OOB scores.This probably means too few trees were used to compute any reliable oob estimates.\n"),warn_flag=true);
                oob_count[0][i]=1;
            }
            for (k=0; k<n_classes; k++) {
                oob_prediction[0][i*n_classes+k]/=oob_count[0][i];
            }
        }
        score2label(oob_prediction[0], oob_label_prediction, n_samples);
        oob_scores=Accuracy(original_y, oob_label_prediction, n_samples);
        if (1) {
            fprintf(stderr, "Out-of-bag score(Accuracy)=%lf\n",oob_scores);
        }
        for (i=0; i<n_jobs; i++) {
            delete []oob_count[i];
            delete []oob_prediction[i];
        }
        delete []oob_count;
        delete []oob_prediction;
        delete []oob_label_prediction;
    }
    if (compute_importance) {
        for (i=0; i<n_features; i++) {
            importances[i]=0;
            for (j=0; j<n_jobs; j++) {
                importances[i]+=importance_jobs[j][i];
            }
            importances[i]/=n_trees;
        }
        for (i=0; i<n_jobs; i++) {
            delete []importance_jobs[i];
        }
        delete []importance_jobs;
    }
    delete []y;
    end=time(NULL);
    fprintf(stderr, "|Random Forest training done. | Using time: %.0lf secs|\n",difftime(end, beg));
    return ENSEMBLE_SUCCESS;
}

void RandomForestClassifier::convert_to_unique(REAL *original_y, REAL *y,uint nSamples){
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

void RandomForestClassifier::score2label(REAL *score, REAL *pred_label, uint nSamples){
    uint idx,i,j;
    for (i=0; i<nSamples; i++) {
        REAL prob=-1;
        for (j=0; j<n_classes; j++) {
            if (score[i*n_classes+j]>prob) {
                prob=score[i*n_classes+j];
                idx=j;
            }
        }
        pred_label[i]=original_classes[idx];
    }
}

void RandomForestClassifier::predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures){
    uint i,j;
    REAL* single_tree_pred=new REAL[nSamples*n_classes];
    REAL* score=new REAL[nSamples*n_classes];
#ifdef DEBUG
    assert(nFeatures==this->n_features);
#endif
    for (i=0; i<nSamples*n_classes; i++) {
        score[i]=0;
    }
    for (i=0; i<n_trees; i++) {
#ifdef DEBUG
        assert(tree[i]!=NULL);
#endif
        tree[i]->predict(X, single_tree_pred, nSamples, nFeatures);
        for (j=0; j<nSamples*n_classes; j++) {
            score[j]+=single_tree_pred[j];
        }
    }
    
    //average prediction
    for (i=0; i<nSamples*n_classes; i++) {
        score[i]/=n_trees;
    }
//    predict_prob(X, score, nSamples, nFeatures);
    score2label(score, pred, nSamples);
    delete [] score;
    delete [] single_tree_pred;
}

void RandomForestClassifier::predict_prob(REAL **X, REAL *pred_prob, uint nSamples, uint nFeatures){
    uint i,j;
    REAL* single_tree_pred=new REAL[nSamples*n_classes];
#ifdef DEBUG
    assert(nFeatures==this->n_features);
#endif
    for (i=0; i<nSamples*n_classes; i++) {
        pred_prob[i]=0;
    }
    for (i=0; i<n_trees; i++) {
#ifdef DEBUG
        assert(tree[i]!=NULL);
#endif
        tree[i]->predict(X, single_tree_pred, nSamples, nFeatures);
        for (j=0; j<nSamples*n_classes; j++) {
            pred_prob[j]+=single_tree_pred[j];
        }
    }
    //average prediction
    for (i=0; i<nSamples*n_classes; i++) {
        pred_prob[i]/=n_trees;
    }
    delete []single_tree_pred;
}


int RandomForestClassifier::save_model(const char *filename){
    FILE* fp=fopen(filename,"w");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for save RandomForestClassifier model.Please check your file path is correct.\n",filename);
        return false;
    }
    fprintf(fp, "RandomForestClassifier\n");
    BaseForest::save_model(fp);
    fclose(fp);
    return true;
}

#endif
