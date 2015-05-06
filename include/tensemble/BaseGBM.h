/* * * * *
 *  BaseGBM.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_BaseGBM_h
#define libTM_BaseGBM_h
#include "LossFunction.h"
#include "Tree.h"
#include "EvaluateMetric.h"
#define GBM_TREE_MAX_DEEP 15

/* define loss */
#define SQUARE_LOSS 0
#define BINOMIAL_DEVIANCE   1
#define MULTINOMIAL_DEVIANCE    2


class BaseGBM {
    /* Base gradient boosting machine,
     * Don't use it directly.
     */
    
public:
    BaseGBM();
    
    BaseGBM(int loss_function,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,REAL subsample,REAL learning_rate,bool oob,bool compute_importance,uint random_seed,uint n_jobs,int verbose);
    
    virtual ~BaseGBM();
    
    virtual int build(REAL **X, REAL *y, uint n_samples,
                      REAL** val_X=NULL,REAL* val_y=NULL,uint n_val_samples=0)=0;    
    
    virtual void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures,uint k_trees)=0;
    
    //may useful for validation
    virtual void predict_verbose(REAL** X,REAL* y,REAL* pred,uint nSamples,uint nFeatures,\
                                 uint k_trees)=0;
    
    virtual int save_model(const char* filename)=0;
    
    bool load_model(const char* filename);
    
    REAL* GetImportances();
    
protected:
    void save_model(FILE* fp=NULL);
    bool check_parameters();
#ifdef DEBUG
    void print_parameters();
#endif
    
public:
    int loss_function;
    LossFunction* loss;
    uint n_classes;
    
    Tree** tree;
    uint n_jobs;
    uint n_trees;
    uint n_features;
    uint max_depth;
    uint min_sample_leaf;
    uint max_features;
    REAL subsample;
    REAL learning_rate;
    REAL* importances;
    int oob;
    int compute_importance;
    uint random_seed;
    int verbose;
    
    /* for regresson prior_pred has size 1,and original_classes=NULL
     * for two-classes problems,prior_pred is the prior probability of positive classes.
     */
    REAL *prior_pred;
    int *original_classes;// converted classes=>original_classes
};

BaseGBM::BaseGBM(){
    tree=NULL;
    loss=NULL;
    importances=NULL;
    prior_pred=NULL;
    original_classes=NULL;
}

BaseGBM::BaseGBM(int loss_function,uint n_trees,uint n_features,uint max_depth,uint min_sample_leaf, REAL max_features_ratio,REAL subsample,REAL learning_rate,bool oob,bool compute_importance,uint random_seed,uint n_jobs,int verbose){
#ifdef DEBUG
    assert(subsample<=1.0 && subsample>0.0);
    assert(loss_function==SQUARE_LOSS || loss_function==MULTINOMIAL_DEVIANCE || loss_function==BINOMIAL_DEVIANCE);
    assert(max_features_ratio<=1.0 && max_features_ratio>0.0);
#endif
    if (max_features_ratio==MTRY_DEFAULT) {
        this->max_features=(uint)sqrt((REAL)n_features);
    }else {
        this->max_features=(uint)MAX(1.0,MIN(max_features_ratio,1.0)*n_features);
    }
    if (max_depth>GBM_TREE_MAX_DEEP) {
        max_depth=GBM_TREE_MAX_DEEP;
    }
    //subsample rate must be less than or equal to 1.0 
    if (subsample>1.0 || subsample<=0.0) {
        subsample=1.0;
    }
    int max_thread=boost::thread::hardware_concurrency();
    if (n_jobs>max_thread) {
        fprintf(stderr, "WARNNING: Number of thread is exceed the maximum number of hardware concurrency.Now set number of thread = %d\n",max_thread);  
        n_jobs=max_thread;
    }
    //reset number of threads if max_features is too small
    if (n_jobs>1 && max_features/n_jobs<THREAD_MIN_FEATURES) {
        fprintf(stderr, "WARNNING: Number of thread is %d,but max_feature=%d.Each thread at least need %d features to find best split for internal nodes in decision tree.Now set number of thread = %d\n",n_jobs,max_features,THREAD_MIN_FEATURES,MAX(1,(uint)max_features/THREAD_MIN_FEATURES));
        n_jobs=MAX(1,(uint)max_features/THREAD_MIN_FEATURES);
    }
    //if compute oob scores,must bootstrap=true
    if (oob && subsample>=1.0) {
        fprintf(stderr, "Out of bag estimation only available if subsample<1.0,auto turn out of bag estimation off.\n");
        oob=false;
    }
    this->loss_function=loss_function;
    this->n_classes=-1;
    this->tree=NULL;
    this->importances=NULL;
    this->loss=NULL;
    this->prior_pred=NULL;
    this->original_classes=NULL;
    this->n_jobs=n_jobs;
    this->n_trees=n_trees;
    this->n_features=n_features;
    this->max_depth=max_depth;
    this->min_sample_leaf=min_sample_leaf;
    this->learning_rate=learning_rate;
    this->subsample=subsample;
    this->oob=oob;
    this->compute_importance=compute_importance;
    this->random_seed=random_seed;
    this->verbose=verbose;
    if (compute_importance) {
        importances=new REAL[n_features];
    }
}
BaseGBM::~BaseGBM(){
    if (tree) {
#ifdef DEBUG
        if (loss_function!=MULTINOMIAL_DEVIANCE) {
            assert(n_classes==1);
        }
#endif
        for (uint i=0; i<n_trees*n_classes; i++) {
            delete tree[i];
            tree[i]=NULL;
        }
        delete []tree;
        tree=NULL;
    }
    if (importances) {
        delete [] importances;
        importances=NULL;
    }
    if (loss) {
        delete loss;
        loss=NULL;
    }
    if (prior_pred) {
        delete []prior_pred;
        prior_pred=NULL;
    }
    if (original_classes) {
        delete []original_classes;
        original_classes=NULL;
    }
}

void BaseGBM::save_model(FILE* fp){
    /* SAVE GBM MODEL
     * FORMAT:
     * loss_function
     * n_classes
     * prior_pred
     * original_classes(optional)
     * n_jobs
     * n_trees
     * n_features
     * max_depth
     * min_sample_leaf
     * max_features
     * subsample
     * learning_rate
     * oob
     * compute_importance
     * random_seed
     * verbose
     * tree#1
     * tree#2
     * ...
     * tree#n
     */
    if (!fp) {
        fprintf(stderr, "Invalid FILE handler for save GBM model.\n");
        return;
    }
    fprintf(fp, "loss_function %d\n",loss_function);
    fprintf(fp, "n_classes %d\n",n_classes);
    fprintf(fp, "prior_pred ");
    for (uint i=0; i<n_classes; i++) {
        fprintf(fp, " %lf",prior_pred[i]);
    }
    fprintf(fp, "\n");
    if (loss_function==BINOMIAL_DEVIANCE) {
        n_classes=2;
    }
    if (loss_function==BINOMIAL_DEVIANCE || loss_function==MULTINOMIAL_DEVIANCE) {
        fprintf(fp, "original_classes ");
        for (uint i=0; i<n_classes; i++) {
            fprintf(fp, " %d",original_classes[i]);
        }
        fprintf(fp, "\n");
    }
    if (loss_function==BINOMIAL_DEVIANCE) {
        n_classes=1;
    }
    fprintf(fp, "n_jobs %d\n",n_jobs);
    fprintf(fp, "n_trees %d\n",n_trees);
    fprintf(fp, "n_features %d\n",n_features);
    fprintf(fp, "max_depth %d\n",max_depth);
    fprintf(fp, "min_sample_leaf %d\n",min_sample_leaf);
    fprintf(fp, "max_features %d\n",max_features);
    fprintf(fp, "subsample %lf\n",subsample);
    fprintf(fp, "learning_rate %lf\n",learning_rate);
    fprintf(fp, "oob %d\n",oob);
    fprintf(fp, "compute_importance %d\n",compute_importance);
    fprintf(fp, "random_seed %d\n",random_seed);
    fprintf(fp, "verbose %d\n",verbose);
    fprintf(fp, "trees\n");
    for (uint i=0; i<n_trees; i++) {
        for (uint j=0; j<n_classes; j++) {
            Tree *t=tree[i*n_classes+j];
            fprintf(fp, "nodes %d\n",t->numNodes);
            for (uint k=ROOT; k<=t->numNodes; k++) {
                fprintf(fp, "%u %.20lg %u %u %u %lg %lg %d ",\
                        t->nodes[k]->feature_split,t->nodes[k]->value_split,\
                        t->nodes[k]->nSamples,t->nodes[k]->left_child,\
                        t->nodes[k]->right_child,t->nodes[k]->ini_error,\
                        t->nodes[k]->best_error,t->nodes[k]->leaf);
                if (t->nodes[k]->leaf) {
                    fprintf(fp, "%.20lg ",t->nodes[k]->pred[0]);
                }
            }
            fprintf(fp, "\n");
        }
    }
}
bool BaseGBM::load_model(const char *filename){
    
    int is_leaf;
    FILE* fp=fopen(filename,"r");
    if (!fp) {
        fscanf(stderr, "Error: Cannot open file %s for load GBM model.Please check your model file.\n",filename);
        return false;
    }
    int parameter_count=0;
    char cmd[100];
    fscanf(fp, "%*s");
    while (1) {
        fscanf(fp, "%100s",cmd);
        if (strcmp(cmd, "loss_function")==0) {
            parameter_count++;
            fscanf(fp, "%d",&loss_function);
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
        else if(strcmp(cmd, "learning_rate")==0){
            parameter_count++;
            fscanf(fp, "%lf",&learning_rate);
        }
        else if(strcmp(cmd, "subsample")==0){
            parameter_count++;
            fscanf(fp, "%lf",&subsample);
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
        else if(strcmp(cmd, "compute_importance")==0){
            parameter_count++;
            fscanf(fp, "%d",&compute_importance);
        }
        else if(strcmp(cmd, "original_classes")==0){
            parameter_count++;
            original_classes=new int[MAX(2,n_classes)];
            for (uint i=0; i<MAX(2,n_classes); i++) {
                fscanf(fp, "%d",&original_classes[i]);
            }
        }
        else if(strcmp(cmd, "prior_pred")==0){
            parameter_count++;
            prior_pred=new REAL[n_classes];
            for (uint i=0; i<n_classes; i++) {
                fscanf(fp, "%lf",&prior_pred[i]);
            }
        }
        else if(strcmp(cmd,"trees")==0){
            break;
        }
        else {
            fprintf(stderr, "Error: GBM read model failed,unknowed parameter %s.model file %s is corrupted.\n",cmd,filename);
            return false;
        }
    }
    if (!check_parameters()) {
        fprintf(stderr, "Error: GBM read model failed,model file %s is corrupted.\n", filename);
        return false;
    }
    //allocate trees
    tree=new Tree*[n_trees*n_classes];
    for (uint i=0; i<n_trees; i++) {
        tree[i]=NULL;
    }
    for (uint i=0; i<n_trees; i++) {
        for (uint j=0; j<n_classes; j++) {
            fscanf(fp, "%100s",cmd);
            if (strcmp(cmd, "nodes")!=0) {
                fprintf(stderr, "Error: GBM read model failed,model file %s is corrupted.\n",filename);
                return false;
            }
            Tree *t=new Tree(CRITERION_MSE, 1, n_features, max_features, min_sample_leaf, max_depth,FIND_BEST,random_seed,n_jobs);
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
                    t->nodes[k]->pred=new REAL;
                    fscanf(fp, "%lf ",&t->nodes[k]->pred[0]);
                }
            }
            tree[i*n_classes+j]=t;
        }
    }
    fclose(fp);
    return true;
}
REAL* BaseGBM::GetImportances(){
    return this->importances;
}
bool BaseGBM::check_parameters(){
    bool ret=true;
    if (n_classes<=0 || n_trees<=0 || n_features<=0 || max_features>n_features || n_jobs<1) {
        ret=false;
    }
    if (min_sample_leaf<=0 ||max_depth<2 || subsample<=0.0 || subsample>1.0 || learning_rate<=0.0 || learning_rate>=1.0) {
        ret=false;
    }
    if (loss_function!=SQUARE_LOSS && loss_function!=BINOMIAL_DEVIANCE && loss_function!=MULTINOMIAL_DEVIANCE) {
        ret=false;
    }
    return ret;
}
#ifdef DEBUG
void BaseGBM::print_parameters(){
    fprintf(stderr, "n_classes %d\n",n_classes);
    fprintf(stderr, "n_jobs %d\n",n_jobs);
    fprintf(stderr, "n_trees %d\n",n_trees);
    fprintf(stderr, "n_features %d\n",n_features);
    fprintf(stderr, "max_depth %d\n",max_depth);
    fprintf(stderr, "min_sample_leaf %d\n",min_sample_leaf);
    fprintf(stderr, "max_features %d\n",max_features);
    fprintf(stderr, "loss_function %d\n",loss_function);
    fprintf(stderr, "learning_rate %lf\n",learning_rate);
    fprintf(stderr, "subsample %f\n",subsample);
    fprintf(stderr, "compute_importance %d\n",compute_importance);
    fprintf(stderr, "random_seed %d\n",random_seed);
    fprintf(stderr, "verbose %d\n",verbose);
}
#endif
#endif
