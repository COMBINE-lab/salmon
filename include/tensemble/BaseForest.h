/* * * * *
 *  BaseForest.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_Forest_h
#define libTM_Forest_h
#include "TreeRegressor.h"
#include "TreeClassifier.h"
#include "algorithm"
#include "EvaluateMetric.h"
#include "stdlib.h"
#define FOREST_TREE_MAX_DEPTH 20

boost::mutex RF_mutex;
boost::mutex random_number_mutex;
uint RF_n_trees_done;

class BaseForest;


class BaseForest {
    //Base class for random forests
    
public:
    BaseForest();
    
    BaseForest(int split_criterion,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,uint find_split_algorithm,bool bootstrap,bool oob,bool compute_importance,uint random_seed,uint n_jobs,bool verbose);
    
    virtual ~BaseForest();
    
    virtual int build(REAL** X,REAL* y,uint n_samples)=0;
    
    virtual void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures)=0;
    
    REAL* GetImportances();
    
    virtual int save_model(const char* filename)=0;
    
    bool load_model(const char* filename);
protected:
    void save_model(FILE* fp=NULL);
    bool check_parameters();
#ifdef DEBUG
    void print_parameters();
#endif
    
public:
    Tree** tree;
    uint n_jobs;
    uint n_trees;
    uint n_features;
    uint n_classes;
    uint max_depth;
    uint min_sample_leaf;
    uint max_features;
    uint find_split_algorithm;
    REAL oob_scores;
    REAL* importances;
    int bootstrap;
    int oob;
    int compute_importance;
    uint random_seed;
    int verbose;

    int split_criterion;
    int *original_classes;// converted classes=>original_classes
};
BaseForest::BaseForest(){
    this->tree=NULL;
    this->importances=NULL;
    this->original_classes=NULL;
}
BaseForest::BaseForest(int split_criterion,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,uint find_split_algorithm,bool bootstrap,bool oob,bool compute_importance,uint random_seed,uint n_jobs,bool verbose){
    
    RF_n_trees_done=0;
    int max_thread=boost::thread::hardware_concurrency();
    if (n_jobs>n_trees) {
        n_jobs=n_trees;
    }
    if (n_jobs>max_thread) {
        fprintf(stderr, "WARN: Number of thread is exceed the maximum number of hardware concurrency.Now set number of thread = %d\n",max_thread);
        n_jobs=max_thread;
    }
    this->tree=NULL;
    this->importances=NULL;
    this->original_classes=NULL;
    this->split_criterion=split_criterion;
    this->tree=new Tree*[n_trees];
    this->n_jobs=n_jobs;
    this->n_trees=n_trees;
    this->n_features=n_features;
    this->find_split_algorithm=find_split_algorithm;
    this->max_depth=MIN(FOREST_TREE_MAX_DEPTH, max_depth);
    this->min_sample_leaf=min_sample_leaf;
    this->bootstrap=bootstrap;
    //if compute oob scores,must bootstrap=true
    if (oob && !bootstrap) {
        fprintf(stderr, "Out of bag estimation only available if bootstrap=True,auto turn out of bag estimation off.\n");
        oob=false;
    }
    this->oob=oob;
    this->oob_scores=-1;
    this->compute_importance=compute_importance;
    this->random_seed=random_seed;
    this->verbose=verbose;
#ifdef DEBUG
    assert(max_features_ratio<=1.0);
#endif
    if (max_features_ratio==MTRY_DEFAULT) {
        this->max_features=sqrt(n_features);
    }else {
        this->max_features=MAX(1,MIN(max_features_ratio,1.0)*n_features);
    }
    if (compute_importance) {
        importances=new REAL[n_features];
    }
    //initialize all tree to null
    for (uint i=0; i<n_trees; i++) {
        this->tree[i]=NULL;
    }
}
BaseForest::~BaseForest(){
    if (tree) {
        for (uint i=0; i<n_trees; i++) {
            if (tree[i]) {
                delete tree[i];
                tree[i]=NULL;
            }
        }
        delete []tree;
        tree=NULL;
    }
    if (importances) {
        delete [] importances;
        importances=NULL;
    }
    if (original_classes) {
        delete []original_classes;
        original_classes=NULL;
    }
}
REAL* BaseForest::GetImportances(){
    return this->importances;
}

void BaseForest::save_model(FILE *fp){
    if (!fp) {
        fprintf(stderr, "Invalid FILE handler for save GBM model.\n");
        return;
    }
    fprintf(fp, "split_criterion %d\n",split_criterion);
    fprintf(fp, "n_classes %d\n",n_classes);
    if (split_criterion==CRITERION_ENTROPY || split_criterion==CRITERION_GINI) {
#ifdef DEBUG
        assert(n_classes>=2);
#endif
        fprintf(fp, "original_classes");
        for (uint i=0; i<n_classes; i++) {
            fprintf(fp, " %d ",original_classes[i]);
        }
        fprintf(fp, "\n");
    }else {
        assert(n_classes==1);
    }
    fprintf(fp, "n_jobs %d\n",n_jobs);
    fprintf(fp, "n_trees %d\n",n_trees);
    fprintf(fp, "n_features %d\n",n_features);
    fprintf(fp, "max_depth %d\n",max_depth);
    fprintf(fp, "min_sample_leaf %d\n",min_sample_leaf);
    fprintf(fp, "max_features %d\n",max_features);
    fprintf(fp, "find_split_algorithm %d\n",find_split_algorithm);
    fprintf(fp, "bootstrap %d\n",bootstrap);
    fprintf(fp, "oob %d\n",oob);
    if (oob) {
        fprintf(fp, "oob_score %lf\n",oob_scores);
    }
    fprintf(fp, "compute_importance %d\n",compute_importance);
    fprintf(fp, "random_seed %d\n",random_seed);
    fprintf(fp, "verbose %d\n",verbose);
    
    fprintf(fp, "trees\n");
    for (uint i=0; i<n_trees; i++) {
        Tree *t=tree[i];
        fprintf(fp, "nodes %d\n",t->numNodes);
        for (uint k=ROOT; k<=t->numNodes; k++) {
            fprintf(fp, "%u %.16lg %u %u %u %lg %lg %d ",\
                    t->nodes[k]->feature_split,t->nodes[k]->value_split,\
                    t->nodes[k]->nSamples,t->nodes[k]->left_child,\
                    t->nodes[k]->right_child,t->nodes[k]->ini_error,\
                    t->nodes[k]->best_error,t->nodes[k]->leaf);
            if (t->nodes[k]->leaf) {
                for (uint j=0; j<n_classes; j++) {
                    fprintf(fp, "%.16lg ",t->nodes[k]->pred[j]);
                }
            }
        }
        fprintf(fp, "\n");
    }
}

bool BaseForest::load_model(const char *filename){
    int is_leaf;
    original_classes=NULL;
    FILE* fp=fopen(filename,"r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open file %s for load RandomForest model.Please check your model file.\n",filename);
        return false;
    }
    int parameter_count=0;
    char cmd[100];
    fscanf(fp, "%*s");
    while (1) {
        fscanf(fp, "%100s",cmd);
        if (strcmp(cmd, "split_criterion")==0) {
            parameter_count++;
            fscanf(fp, "%d",&split_criterion);
        }
        else if(strcmp(cmd, "n_classes")==0){
            parameter_count++;
            fscanf(fp, "%d",&n_classes);
        }
        else if(strcmp(cmd, "n_jobs")==0){
            parameter_count++;
            fscanf(fp, "%d",&n_jobs);
        }
        else if(strcmp(cmd, "n_trees")==0){
            parameter_count++;
            fscanf(fp, "%d",&n_trees);
        }
        else if(strcmp(cmd, "n_features")==0){
            parameter_count++;
            fscanf(fp, "%d",&n_features);
        }
        else if(strcmp(cmd, "max_depth")==0){
            parameter_count++;
            fscanf(fp, "%d",&max_depth);
        }
        else if(strcmp(cmd, "min_sample_leaf")==0){
            parameter_count++;
            fscanf(fp, "%d",&min_sample_leaf);
        }
        else if(strcmp(cmd, "max_features")==0){
            parameter_count++;
            fscanf(fp, "%d",&max_features);
        }
        else if(strcmp(cmd, "find_split_algorithm")==0){
            parameter_count++;
            fscanf(fp, "%d",&find_split_algorithm);
        }
        else if(strcmp(cmd, "bootstrap")==0){
            parameter_count++;
            fscanf(fp, "%d",&bootstrap);
        }
        else if(strcmp(cmd, "oob")==0){
            parameter_count++;
            fscanf(fp, "%d",&oob);
        }
        else if(strcmp(cmd, "random_seed")==0){
            parameter_count++;
            fscanf(fp, "%d",&random_seed);
        }
        else if(strcmp(cmd, "verbose")==0){
            parameter_count++;
            fscanf(fp, "%d",&verbose);
        }
        else if(strcmp(cmd, "oob_score")==0){
//            parameter_count++;
            fscanf(fp, "%lf",&oob_scores);
        }
        else if(strcmp(cmd, "compute_importance")==0){
            parameter_count++;
            fscanf(fp, "%d",&compute_importance);
        }
        else if(strcmp(cmd, "original_classes")==0){
            parameter_count++;
            original_classes=new int[n_classes];
            for (uint i=0; i<n_classes; i++) {
                fscanf(fp, "%d",&original_classes[i]);
            }
        }
        else if(strcmp(cmd,"trees")==0){
            break;
        }
        else {
            fprintf(stderr, "RandomForest read model failed,unknowed parameter %s.model file %s is corrupted.\n",cmd,filename);
            return false;
        }
    }
    if (!check_parameters()) {
        fprintf(stderr, "RandomForest read model failed,model file %s is corrupted.\n", filename);
        return false;
    }
    //allocate trees
    tree=new Tree*[n_trees];
    for (uint i=0; i<n_trees; i++) {
        tree[i]=NULL;
    }
    for (uint i=0; i<n_trees; i++) {
        fscanf(fp, "%100s",cmd);
        if (strcmp(cmd, "nodes")!=0) {
            fprintf(stderr, "RandomForest read model failed,model file %s is corrupted.\n",filename);
            return false;
        }
        Tree *t=new Tree(split_criterion, n_classes, n_features, max_features, min_sample_leaf, max_depth,FIND_BEST,random_seed,1);
        fscanf(fp, "%d\n",&t->numNodes);
        t->max_nNodes=t->numNodes+1;
        t->nodes=new TreeNode*[t->max_nNodes];
        t->nodes[0]=new TreeNode();//dummy node.
        for (uint k=ROOT; k<=t->numNodes; k++) {
            t->nodes[k]=new TreeNode();
            t->nodes[k]->leaf=false;
            fscanf(fp, "%u %lf %u %u %u %lf %lf %d ",\
                   &t->nodes[k]->feature_split,&t->nodes[k]->value_split,\
                   &t->nodes[k]->nSamples,&t->nodes[k]->left_child,\
                   &t->nodes[k]->right_child,&t->nodes[k]->ini_error,\
                   &t->nodes[k]->best_error,&is_leaf);
            if (is_leaf) {
                t->nodes[k]->leaf=true;
                t->nodes[k]->pred=new REAL[n_classes];
                for (uint j=0; j<n_classes; j++) {
                    fscanf(fp, "%lf ",&t->nodes[k]->pred[j]);
                }
            }
        }
        tree[i]=t;
    }
    fclose(fp);
    return true;
}

void RF_build_trees_range(REAL** X,REAL* y,uint n_samples,BaseForest* forest,pair<uint, uint> n_trees_range,uint* oob_count,REAL* oob_prediction,REAL* importances){
#ifdef DEBUG
    if (forest->split_criterion==CRITERION_MSE) {
        assert(forest->n_classes==1);
    }
#endif
    char* str=new char[100];
    int split_criterion=forest->split_criterion;
    Tree** tree=forest->tree;
    uint n_features=forest->n_features;
    uint max_features=forest->max_features;
    uint max_depth=forest->max_depth;
    uint min_sample_leaf=forest->min_sample_leaf;
    uint find_split_algorithm=forest->find_split_algorithm;
    uint random_seed=forest->random_seed;
    uint n_classes=forest->n_classes;
    uint n_trees_beg=n_trees_range.first;
    uint n_trees_end=n_trees_range.second;
    uint verbose=forest->verbose;
    bool oob=forest->oob;
    bool compute_importance=forest->compute_importance;
    bool bootstrap=forest->bootstrap;
    
    uint i,j,k;
    REAL *oob_prediction_tmp,*importance_tmp;
    REAL** sub_X;
    REAL* sub_y;
    bool* mask;//for oob prediction
    if (oob) {
        oob_prediction_tmp=new REAL[n_samples*n_classes];
        for (i=0; i<n_samples*n_classes; i++) {
            oob_prediction[i]=0;
        }
        for (i=0; i<n_samples; i++){
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
        boost::mutex::scoped_lock rand_lock(random_number_mutex);
        uint job_random_seed=randomMT();
        rand_lock.unlock();
        if (oob) {
            for (j=0; j<n_samples; j++) {
                mask[j]=false;
            }
        }
        if (bootstrap) {
            for (j=0; j<n_samples; j++) {
                uint idx=rand_r(&job_random_seed)%n_samples;
//                uint idx=randomMT()%n_samples;
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
        //build tree
        if (split_criterion==CRITERION_MSE) {
            tree[i]=new TreeRegressor(n_features,\
                                      max_features,\
                                      min_sample_leaf,\
                                      max_depth,\
                                      find_split_algorithm,\
                                      job_random_seed,\
                                      1);
        }else{
            tree[i]=new TreeClassifier(split_criterion,\
                                       n_classes,\
                                       n_features,\
                                       max_features,\
                                       min_sample_leaf,\
                                       max_depth,\
                                       find_split_algorithm,\
                                       job_random_seed,\
                                       1);
        }
        tree[i]->build(sub_X, sub_y, n_samples);
        if (oob) {
            tree[i]->predict(X, oob_prediction_tmp, mask, n_samples, n_features);
            for (j=0; j<n_samples; j++) {
                if (mask[j]) {
                    continue;
                }
                oob_count[j]++;
                for (k=0; k<n_classes; k++) {
                    oob_prediction[j*n_classes+k]+=oob_prediction_tmp[j*n_classes+k];
                }
            }
        }
        if (compute_importance) {
            tree[i]->compute_importance(importance_tmp);
            for (j=0; j<n_features; j++) {
                importances[j]+=importance_tmp[j];
            }
            
        }
        
        //print info
        boost::mutex::scoped_lock lock(RF_mutex);
        RF_n_trees_done++;
        if (verbose) {
            if (RF_n_trees_done!=1) {
                for (j=0; j<50; j++) {
                    fprintf(stderr, "\b");
                }
            }
            str[0]='\0';
            sprintf(str, "Random forest progress: %d/%d",RF_n_trees_done,forest->n_trees);
            int str_len=strlen(str);
            fprintf(stderr, "%s",str);
            for (j=str_len; j<50; j++) {
                fprintf(stderr, " ");
            }
        }
        lock.unlock();
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
    delete []str;
}
bool BaseForest::check_parameters(){
    bool ret=true;
    if (n_classes<=0 || n_trees<=0 || n_features<=0 || max_features>n_features || n_jobs<1) {
        ret=false;
    }
    if (min_sample_leaf<=0 ||max_depth<2 || ( find_split_algorithm!=FIND_BEST && find_split_algorithm!=FIND_RANDOM) ) {
        ret=false;
    }
    return ret;
}
#ifdef DEBUG
void BaseForest::print_parameters(){
    fprintf(stderr, "n_classes %d\n",n_classes);
    fprintf(stderr, "n_jobs %d\n",n_jobs);
    fprintf(stderr, "n_trees %d\n",n_trees);
    fprintf(stderr, "n_features %d\n",n_features);
    fprintf(stderr, "max_depth %d\n",max_depth);
    fprintf(stderr, "min_sample_leaf %d\n",min_sample_leaf);
    fprintf(stderr, "max_features %d\n",max_features);
    fprintf(stderr, "find_split_algorithm %d\n",find_split_algorithm);
    fprintf(stderr, "bootstrap %d\n",bootstrap);
    fprintf(stderr, "oob %d\n",oob);
    if (oob) {
        fprintf(stderr, "oob_score %lf\n",oob_scores);
    }
    fprintf(stderr, "compute_importance %d\n",compute_importance);
    fprintf(stderr, "random_seed %d\n",random_seed);
    fprintf(stderr, "verbose %d\n",verbose);
}
#endif
#endif
