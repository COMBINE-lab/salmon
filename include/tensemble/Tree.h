/* * * * *
 *  Tree.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_MultiThread_Tree_h
#define libTM_MultiThread_Tree_h

#include <iostream>
#include <string>
#include "TreeNode.h"
#include "FeatureData.h"
#include "Criterion.h"
#include "ClassificationCriterion.h"
#include <cmath>
#include "RandomGenerator.h"
#include "assert.h"

#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#define TREE_MAX_DEPTH  30
//using namespace boost;

enum SplitAlgorithm {
    FIND_BEST = 0,
    FIND_RANDOM = 1
};

struct MyPair {
    uint index;
    REAL value;
    bool friend operator < (const MyPair& a,const MyPair& b){
        return a.value<b.value;
    }
};

//split_args for parallel find best split,because boost:bind can handle 8 parameter at most
struct split_args {
    uint* sample_ind;
    uint s_ind_beg;
    uint s_ind_end;
    uint f_beg;
    uint f_end;
    uint nSamples;
    uint min_samples_leaf;
};
//for parallel find_best_split
void find_best_split_range(bool *skip,Criterion *criterion,REAL** X,REAL* y,split_args args,uint& feature_split,REAL& value_split,REAL& min_error);

class Tree {
    /*base decision tree*/

public:
    Criterion **p_criterion;
    int criterion_name;
    uint min_sample_leaf;
    uint n_features;
    uint n_classes;//n_classes=1 for regression
    uint max_features;
    uint max_depth;
    uint find_split_algorithm;
    uint random_seed;

    TreeNode **nodes;
    sint *nodes_deep;
    uint numNodes;
    uint max_nNodes;

    uint n_threads;

public:
    Tree(int criterion,uint n_classes,uint n_feature,uint max_features,uint min_sample_leaf,uint max_depth,uint find_split_algorithm=FIND_BEST,\
         uint random_seed=0,uint n_thread=1);

    virtual ~Tree();

    int build(REAL** X,REAL* y,uint nSamples);

    bool find_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error);

    bool find_best_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error);

    bool find_random_split(uint node_id,uint& feature_split,REAL& value_split,REAL& min_error);

    virtual void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures);

    //predict with mask
    virtual void predict(REAL** X,REAL* pred,bool* mask,uint nSamples,uint nFeatures);

    void predict_terminal_region(REAL** X,uint* region,uint nSamples,uint nFeatures);

    virtual void compute_importance(REAL* importance){};

    //for parallel find best split
    bool find_best_split_parallel(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error);

    void resize();

protected:
    void init_nodes(uint nSamples);

    void init_criterion();

    uint updateNodeSampleMap(uint node_id,uint l_node_id,uint r_node_id,\
                             uint* sample_ind,uint* sample_ind_swap,REAL** X,\
                             uint s_ind_beg,uint s_ind_end,uint feature_split,\
                             REAL value_split);
};

Tree::Tree(int criterion,uint n_classes,uint n_features,uint max_features,uint min_sample_leaf,uint max_depth,uint find_split_algorithm,uint random_seed,uint n_threads){
#ifdef DEBUG
    assert(criterion==CRITERION_MSE || criterion==CRITERION_GINI || criterion==CRITERION_ENTROPY);
    assert(find_split_algorithm==FIND_RANDOM || find_split_algorithm==FIND_BEST);
    assert(max_features<=n_features);
    assert(n_threads>=1);
    if (criterion==CRITERION_MSE) {
        assert(n_classes==1);
    }
#endif
    if (max_features>n_features || max_features<=0) {
        max_features=n_features;
    }
    int max_thread=boost::thread::hardware_concurrency();
    if (n_threads>max_thread) {
        n_threads=max_thread;
    }
    max_features=MAX(1,max_features);

    //reset number of threads if max_features is too small
    if (n_threads>1 && max_features/n_threads<THREAD_MIN_FEATURES) {
//        fprintf(stderr, "WARNNING:Number of thread is %d,but max_feature=%d.Each thread at least need %d features to find best split for internal nodes in decision tree.Now set number of thread = %d\n",n_threads,max_features,THREAD_MIN_FEATURES,MAX(1,(uint)max_features/THREAD_MIN_FEATURES));
        n_threads=MAX(1,(uint)max_features/THREAD_MIN_FEATURES);
    }
    this->n_threads=n_threads;
    this->criterion_name=criterion;
    switch (find_split_algorithm) {
        case FIND_BEST:
            this->find_split_algorithm=FIND_BEST;
            break;
        case FIND_RANDOM:
            this->find_split_algorithm=FIND_RANDOM;
            break;
        default:
            break;
    }
    this->n_features=n_features;
    this->n_classes=n_classes;
    this->max_depth=MIN(TREE_MAX_DEPTH, max_depth);
    this->max_features=max_features;
    this->min_sample_leaf=min_sample_leaf;
    this->random_seed=random_seed;
    this->numNodes=0;
    this->nodes=NULL;
    this->p_criterion=NULL;
    this->nodes_deep=NULL;
}
Tree::~Tree(){
    if (nodes) {
        for (uint i=0; i<max_nNodes; i++) {
            delete nodes[i];
            nodes[i]=NULL;
        }
        delete []nodes;
        nodes=NULL;
    }
    if (p_criterion){
        for (uint i=0; i<n_threads; i++) {
            delete p_criterion[i];
            p_criterion[i]=NULL;
        }
        delete []p_criterion;
        p_criterion=NULL;
    }
    if (nodes_deep) {
        delete []nodes_deep;
        nodes_deep=NULL;
    }
}
void Tree::resize(){
    if ((REAL)max_nNodes/numNodes<1.3) {
//        cout<<max_nNodes<<endl;
//        cout<<(REAL)max_nNodes/numNodes<<endl;
//        cout<<numNodes<<endl;
        return;
    }
    if (nodes) {
        for (uint i=numNodes+1; i<max_nNodes; i++) {
            delete nodes[i];
            nodes[i]=NULL;
        }
        max_nNodes=numNodes+1;
    }
}
void Tree::init_criterion(){
    this->p_criterion=new Criterion*[n_threads];
    switch (criterion_name) {
        case CRITERION_MSE:
            for (uint i=0; i<n_threads; i++) {
                p_criterion[i]=new MSE;
            }
            break;
        case CRITERION_GINI:
            for (uint i=0; i<n_threads; i++) {
                p_criterion[i]=new Gini(n_classes);
            }
            break;
        case CRITERION_ENTROPY:
            for (uint i=0; i<n_threads; i++) {
                p_criterion[i]=new Entropy(n_classes);
            }
            break;
        default:
            break;
    }
}
void Tree::init_nodes(uint nSamples){
    numNodes=1;
    /*===============max_nNodes must choose carefully==========================
     * too small will get poor result,too big will waste memory and may slow down the speed
     */
    max_nNodes=2*(uint)((REAL)floor((REAL)(nSamples/(1>((int)min_sample_leaf-4)?1:(min_sample_leaf-4))))+1);
    if (max_nNodes>pow((REAL)2.0,(REAL)max_depth)) {
        max_nNodes=pow((REAL)2.0, (REAL)max_depth)+1;
    }
    nodes=new TreeNode*[max_nNodes];
    nodes_deep=new sint[max_nNodes];
    for (uint i=0; i<max_nNodes; i++) {
        nodes[i]=new TreeNode;
    }
    nodes_deep[ROOT]=1;
}

uint Tree::updateNodeSampleMap(uint node_id,uint l_node_id,uint r_node_id,\
                               uint* sample_ind,uint* sample_ind_swap,REAL** X,\
                               uint s_ind_beg,uint s_ind_end,uint feature_split,\
                               REAL value_split){
    /* Map sample index in parent node to left child and right child
     * return the number of sample in left child.
     */

#ifdef DEBUG
    assert(s_ind_beg>=0 && s_ind_beg<s_ind_end);
    assert(s_ind_end>0);
#endif
    uint i,count=s_ind_beg;
    if (l_node_id>=max_nNodes || r_node_id>=max_nNodes) {
        return -1;
    }
    nodes[l_node_id]->nSamples=0;
    nodes[r_node_id]->nSamples=0;
    for (i=s_ind_beg; i<s_ind_end; i++) {
        uint idx=sample_ind[i];
        if (X[idx][feature_split]<=value_split) {
            nodes[l_node_id]->nSamples++;
            sample_ind_swap[count++]=idx;
        }
    }
    for (i=s_ind_beg; i<s_ind_end; i++) {
        uint idx=sample_ind[i];
        if (X[idx][feature_split]>value_split) {
            nodes[r_node_id]->nSamples++;
            sample_ind_swap[count++]=idx;
        }
    }
#ifdef DEBUG
    assert(count==s_ind_end);
#endif

    for (i=s_ind_beg; i<s_ind_end; i++) {
        sample_ind[i]=sample_ind_swap[i];
    }
    return nodes[l_node_id]->nSamples;
}

int Tree::build(REAL** X,REAL* y,uint nSamples){
    /* build decision tree
     * X,y:the original data
     * nSamples: number of sample in X
     */
#ifdef DEBUG
    if (criterion_name==CRITERION_MSE) {
        assert(n_classes==1);
    }else{
        assert(n_classes>=2);
    }
#endif
    uint i,l_node_id,r_node_id,work_node_id,feature_split,sample_ind_beg,sample_ind_end;
    REAL value_split,ini_error,min_error;
    uint *node_beg,*node_end,*sample_ind,*sample_ind_swap;

    //initialize criterion and nodes
    init_criterion();
    init_nodes(nSamples);
    nodes[ROOT]->nSamples=nSamples;
    node_beg=new uint[max_nNodes];
    node_end=new uint[max_nNodes];
    sample_ind=new uint[nSamples];
    sample_ind_swap=new uint[nSamples];
    node_beg[ROOT]=0;
    node_end[ROOT]=nSamples;
    for (i=0; i<nSamples; i++) {
        sample_ind[i]=i;
    }
    //main
    for (work_node_id=ROOT; work_node_id<max_nNodes; work_node_id++) {
        if (work_node_id>numNodes) {
            break;
        }
        sample_ind_beg=node_beg[work_node_id];
        sample_ind_end=node_end[work_node_id];
        //initialize threads criterion,let every thread has a criterion copy
        for (i=0; i<n_threads; i++) {
            p_criterion[i]->init(y, sample_ind,sample_ind_beg,sample_ind_end, nodes[work_node_id]->nSamples);
        }
        ini_error=p_criterion[0]->eval();
        nodes[work_node_id]->ini_error=ini_error;
        //terminal node if the number of sample in a node is less than min_sample_leaf
        if (nodes[work_node_id]->nSamples<2*min_sample_leaf || nodes_deep[work_node_id]==max_depth || numNodes>=max_nNodes-2) {
//        if (nodes[work_node_id]->nSamples<=min_sample_leaf || nodes_deep[work_node_id]==max_depth || numNodes>=max_nNodes-2) {
            nodes[work_node_id]->leaf=true;
            //allocate memory for prediction
            nodes[work_node_id]->pred=new REAL[n_classes];
            p_criterion[0]->estimate(nodes[work_node_id]->pred);
            continue;
        }
        //find split
        min_error=ini_error;
        bool flag=find_split(work_node_id, X, y,sample_ind, sample_ind_beg,sample_ind_end, nodes[work_node_id]->nSamples, feature_split, value_split,min_error);
        //terminal node
        if (!flag) {
            nodes[work_node_id]->leaf=true;
            p_criterion[0]->reset();
            //allocate memory for prediction
            nodes[work_node_id]->pred=new REAL[n_classes];
            p_criterion[0]->estimate(nodes[work_node_id]->pred);
            continue;
        }
        //found split
        nodes[work_node_id]->best_error=min_error;
        nodes[work_node_id]->feature_split=feature_split;
        nodes[work_node_id]->value_split=value_split;

        //set left,right child
        l_node_id=++numNodes;
        r_node_id=++numNodes;
        nodes[work_node_id]->left_child=l_node_id;
        nodes[work_node_id]->right_child=r_node_id;
        nodes_deep[l_node_id]=nodes_deep[r_node_id]=nodes_deep[work_node_id]+1;

        //map samples in parent node to left,right child
        uint mid=updateNodeSampleMap(work_node_id, l_node_id,r_node_id, sample_ind,sample_ind_swap, X, sample_ind_beg, sample_ind_end, feature_split, value_split);
        node_beg[l_node_id]=sample_ind_beg;
        node_end[l_node_id]=sample_ind_beg + mid;
        node_beg[r_node_id]=sample_ind_beg + mid;
        node_end[r_node_id]=sample_ind_end;
    }
    delete []sample_ind_swap;
    delete []sample_ind;
    delete []node_beg;
    delete []node_end;
    delete []nodes_deep;//no more need nodes_deep
    nodes_deep=NULL;
    //no more need criterion.
    for (uint i=0; i<n_threads; i++) {
        delete p_criterion[i];
    }
    delete []p_criterion;
    p_criterion=NULL;
    resize();
    return ENSEMBLE_SUCCESS;
}

bool Tree::find_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error){
    if (find_split_algorithm==FIND_BEST) {
        if (n_threads>1) {
            return find_best_split_parallel(node_id, X,y, sample_ind, s_ind_beg,s_ind_end,nSamples, feature_split, value_split,min_error);
        }
        return find_best_split(node_id, X,y, sample_ind, s_ind_beg,s_ind_end,nSamples, feature_split, value_split,min_error);
    }
    return -1;
}

bool Tree::find_best_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error){
    /* SINGLE-THREAD FIND SPLIT
     * find best split for internal node(node_id)
     * X,y is the original data
     * samples in internal node(node_id) is specify by spsample_ind[s_ind_beg:s_ind_end]
     */
#ifdef DEBUG
    assert(s_ind_beg>=0 && s_ind_beg<s_ind_end);
    assert((s_ind_end-s_ind_beg) ==nSamples);
#endif
    REAL ini_error=min_error;
    uint i,j,node_nSamples,f,count;
    bool same_data=true;
    //do not consider split if all data are identical.
    REAL pre_y=y[sample_ind[s_ind_beg]];
    for (i=s_ind_beg+1; i<s_ind_end; i++) {
        uint idx=sample_ind[i];
        if (y[idx]!=pre_y) {
            same_data=false;
            break;
        }
    }
    if (same_data) {
        return false;
    }

    Criterion *criterion=p_criterion[0];
    bool* skip=new bool[n_features];
    uint* skip_idx=new uint[n_features];
    node_nSamples=nodes[node_id]->nSamples;
    MyPair* x_f=new MyPair[nSamples];

    //generate random split feature
    for (i=0; i<n_features; i++) {
        skip[i]=true;
        skip_idx[i]=i;
    }
    count=n_features;
    for (i=0; i<max_features; i++) {
        j=rand_r(&this->random_seed)%count;
//        j=randomMT()%count;
        uint idx=skip_idx[j];
        skip[idx]=false;
        std::swap(skip_idx[j], skip_idx[--count]);
    }
    //find best split main loop
    for (f=0; f<n_features; f++) {
        if (skip[f]) {
            continue;
        }
        criterion->reset();
        count=0;
        for (i=s_ind_beg; i<s_ind_end; i++) {
            uint idx=sample_ind[i];
            x_f[count].index=idx;
            x_f[count].value=X[idx][f];
            count++;
        }
#ifdef DEBUG
        assert(count==nSamples);
#endif
        //sort sample with feature f
        sort(x_f, x_f+nSamples);
        if (x_f[0].value>=x_f[nSamples-1].value) {
            continue;
        }
        count=0;
        for (i=s_ind_beg; i<s_ind_end-1; i++) {
            uint loc=x_f[count].index;
            uint nLeft=criterion->update_next(y, loc);
            //do not consider split if feature value are identical.
            if (x_f[count].value==x_f[count+1].value || nLeft<min_sample_leaf || nSamples-nLeft<min_sample_leaf) {
//            if (x_f[count].value==x_f[count+1].value){
                count++;
                continue;
            }
            REAL error=criterion->eval();
            if (error<min_error) {
                min_error=error;
                feature_split=f;
                value_split=0.5*(x_f[count].value+x_f[count+1].value);
            }
            count++;
        }
    }
//    fprintf(stderr, "\t%d,%lf,%lf\n",min_error==ini_error, ini_error, min_error);
    delete []skip;
    delete []skip_idx;
    delete []x_f;
    return min_error!=ini_error;
}

bool Tree::find_best_split_parallel(uint node_id, REAL **X, REAL *y, uint *sample_ind, uint s_ind_beg, uint s_ind_end, uint nSamples, uint &feature_split, REAL &value_split, REAL &min_error){
    /* MULTI-THREAD FIND SPLIT
     * parallel find best split for internal node(node_id)
     * X,y is the original data
     * samples in internal node(node_id) is specify by spsample_ind[s_ind_beg:s_ind_end]
     */
    uint *f_split,i,j,count;
    REAL *v_split,*p_min_error;
    struct split_args args;
    f_split=new uint[n_threads];
    v_split=new REAL[n_threads];
    p_min_error=new REAL[n_threads];
    bool* skip=new bool[n_features];
    uint* skip_idx=new uint[n_features];
    args.s_ind_beg=s_ind_beg;
    args.s_ind_end=s_ind_end;
    args.sample_ind=sample_ind;
    args.nSamples=nSamples;
    args.min_samples_leaf=min_sample_leaf;
    REAL ini_error=min_error;

    //generate random split feature
    for (i=0; i<n_features; i++) {
        skip[i]=true;
        skip_idx[i]=i;
    }
    count=n_features;
    for (i=0; i<max_features; i++) {
        j=rand_r(&this->random_seed)%count;
//        j=randomMT()%count;
        uint idx=skip_idx[j];
        skip[idx]=false;
        std::swap(skip_idx[j], skip_idx[--count]);
    }
    boost::thread **thread=new boost::thread*[n_threads];
    for (i=0; i<n_threads; i++) {
        p_min_error[i]=ini_error;
        args.f_beg=(n_features/n_threads*i);
        args.f_end=(n_features/n_threads*(i+1));
        if (i==n_threads-1) {
            args.f_end=n_features;
        }
        thread[i]=new boost::thread(boost::bind(find_best_split_range,skip, p_criterion[i], X, y, args, boost::ref(f_split[i]), boost::ref(v_split[i]), boost::ref(p_min_error[i])));
    }
    for (i=0; i<n_threads; i++) {
        thread[i]->join();
        delete thread[i];
        if (p_min_error[i]<min_error) {
            min_error=p_min_error[i];
            feature_split=f_split[i];
            value_split=v_split[i];
        }
    }
    delete thread;
    delete []skip;
    delete []skip_idx;
    delete []f_split;
    delete []v_split;
    delete []p_min_error;
    return min_error!=ini_error;
}
void Tree::predict(REAL **X,REAL* pred,uint nSamples,uint nFeatures){
    /* Make Prediction
     * for classification.pred is the probability of each classes.
     */
#ifdef DEBUG
    assert(nFeatures==this->n_features);
    if (criterion_name==CRITERION_MSE) {
        assert(n_classes==1);
    }else {
        assert(n_classes>=2);
    }
#endif
    uint node_id,f_s;
    REAL v_s,normalizer=0.0;
    for (uint i=0; i<nSamples; i++) {
        node_id=ROOT;
        while (!nodes[node_id]->leaf) {
            f_s=nodes[node_id]->feature_split;
            v_s=nodes[node_id]->value_split;
            if (X[i][f_s]<=v_s) {
                node_id=nodes[node_id]->left_child;
            }else {
                node_id=nodes[node_id]->right_child;
            }
        }
        normalizer=0.0;
        for (uint j=0; j<n_classes; j++) {
            pred[n_classes*i+j]=nodes[node_id]->pred[j];
            normalizer+=nodes[node_id]->pred[j];
        }
        if (criterion_name!=CRITERION_MSE) {
            if (normalizer!=0) {
                for (uint j=0; j<n_classes; j++) {
                    pred[n_classes*i+j]/=normalizer;
                }
            }
        }
    }
}
void Tree::predict(REAL **X,REAL* pred,bool* mask,uint nSamples,uint nFeatures){
    /* Make Prediction with mask
     * for classification.pred is the probability of each classes.
     */
#ifdef DEBUG
    assert(nFeatures==this->n_features);
    if (criterion_name==CRITERION_MSE) {
        assert(n_classes==1);
    }
#endif
    uint node_id,f_s;
    REAL v_s,normalizer;
    for (uint i=0; i<nSamples; i++) {
        if (mask[i]) {
            continue;
        }
        node_id=ROOT;
        while (!nodes[node_id]->leaf) {
            f_s=nodes[node_id]->feature_split;
            v_s=nodes[node_id]->value_split;
            if (X[i][f_s]<=v_s) {
                node_id=nodes[node_id]->left_child;
            }else {
                node_id=nodes[node_id]->right_child;
            }
        }
        normalizer=0.0;
        for (uint j=0; j<n_classes; j++) {
            pred[n_classes*i+j]=nodes[node_id]->pred[j];
            normalizer+=nodes[node_id]->pred[j];
        }
        if (criterion_name!=CRITERION_MSE) {
            if (normalizer!=0) {
                for (uint j=0; j<n_classes; j++) {
                    pred[n_classes*i+j]/=normalizer;
                }
            }
        }

    }
}

void Tree::predict_terminal_region(REAL **X, uint* region,uint nSamples, uint nFeatures){
    /* GET TERMINAL REGION FOR EACH SAMPLE */
#ifdef DEBUG
    assert(nFeatures==this->n_features);
#endif
    uint node_id,f_s;
    REAL v_s;
    for (uint i=0; i<nSamples; i++) {
        node_id=ROOT;
        while (!nodes[node_id]->leaf) {
            f_s=nodes[node_id]->feature_split;
            v_s=nodes[node_id]->value_split;
            if (X[i][f_s]<=v_s) {
                node_id=nodes[node_id]->left_child;
            }else {
                node_id=nodes[node_id]->right_child;
            }
        }
        region[i]=node_id;
    }
}

void find_best_split_range(bool *skip,Criterion *criterion,REAL** X,REAL* y,split_args args,uint& feature_split,REAL& value_split,REAL& min_error){
    /* FOR MULTI-THREAD FIND SPLIT
     * find best split for internal node(node_id),feature range in [f_beg,f_end]
     * X,y is the original data
     * samples in internal node(node_id) is specify by spsample_ind[s_ind_beg:s_ind_end]
     * [f_beg,f_end] /in [0,n_features] is the feature range for each thread.
     *
     * notice：
     * criterion must be initialied before calling find_best_split_range
     * one criterion for one thread
     */
    uint s_ind_beg,s_ind_end,nSamples,f_beg,f_end,*sample_ind,min_samples_leaf;
    min_samples_leaf=args.min_samples_leaf;
    s_ind_beg=args.s_ind_beg;
    s_ind_end=args.s_ind_end;
    nSamples=args.nSamples;
    f_beg=args.f_beg;
    f_end=args.f_end;
    sample_ind=args.sample_ind;
    REAL ini_error=min_error;
#ifdef DEBUG
    assert(s_ind_beg>=0 && s_ind_beg<s_ind_end);
    assert((s_ind_end-s_ind_beg) ==nSamples);
#endif
    uint i,f,count;
    bool same_data=true;

    //do not consider split if all data are identical.
    REAL pre_y=y[sample_ind[s_ind_beg]];
    for (i=s_ind_beg+1; i<s_ind_end; i++) {
        uint idx=sample_ind[i];
        if (y[idx]!=pre_y) {
            same_data=false;
            break;
        }
    }
    if (same_data) {
        return;
    }

    MyPair* x_f=new MyPair[nSamples];
    criterion->init(y, sample_ind,s_ind_beg,s_ind_end, nSamples);
    //find best split main loop
    for (f=f_beg; f<f_end; f++) {
        if (skip[f]) {
            continue;
        }
        criterion->reset();
        count=0;
        for (i=s_ind_beg; i<s_ind_end; i++) {
            uint idx=sample_ind[i];
            x_f[count].index=idx;
            x_f[count].value=X[idx][f];
            count++;
        }
#ifdef DEBUG
        assert(count==nSamples);
#endif
        //sort sample with feature f
        sort(x_f, x_f+nSamples);
        if (x_f[0].value>=x_f[nSamples-1].value) {
            continue;
        }
        count=0;
        for (i=s_ind_beg; i<s_ind_end-1; i++) {
            uint loc=x_f[count].index;
            uint nLeft=criterion->update_next(y, loc);
            //do not consider split if feature value are identical.
            if (x_f[count].value==x_f[count+1].value || nLeft<min_samples_leaf || nSamples-nLeft<min_samples_leaf) {
//            if (x_f[count].value==x_f[count+1].value){
                count++;
                continue;
            }
            REAL error=criterion->eval();
            if (error<min_error) {
                min_error=error;
                feature_split=f;
                value_split=0.5*(x_f[count].value+x_f[count+1].value);
            }
            count++;
        }
    }
    delete []x_f;

}
#endif
