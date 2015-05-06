//
//  Tree.h
//  libTM
//
//  Created by apple on 12-10-29.
//  Copyright (c) 2012年 Rongkai Xia. All rights reserved.
//

#ifndef libTM_Tree_h
#define libTM_Tree_h
#include <iostream>
#include <string>
#include "TreeNode.h"
#include "FeatureData.h"
#include "Criterion.h"
#include <math.h>
//#include "cokus.cpp"
#include "RandomGenerator.h"
#include "assert.h"
#define LEFT(A) ((A)<<1)
#define RIGHT(A) (((A)<<1)+1)
#define MAX_DEPTH   100
#define ROOT    1

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

class Tree {
    /*decision tree*/
    
public:
    Criterion *criterion;
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

public:
    Tree(int criterion,uint n_classes,uint n_feature,uint max_features,uint min_sample_leaf,\
         uint max_depth,uint find_split_algorithm,\
         uint random_seed=0);
    
    virtual ~Tree();
    
    int build(REAL** X,REAL* y,uint nSamples);
    
    bool find_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error);
    
    bool find_best_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error);
    
    bool find_random_split(uint node_id,uint& feature_split,REAL& value_split,REAL& min_error);
    
    void predict(REAL** X,REAL* pred,uint nSamples,uint nFeatures);
    
    virtual void compute_importance(REAL* importance)=0;
protected:
    void init_nodes(uint nSamples);
    
    uint updateNodeSampleMap(uint node_id,uint l_node_id,uint r_node_id,\
                             uint* sample_ind,uint* sample_ind_swap,REAL** X,\
                             uint s_ind_beg,uint s_ind_end,uint feature_split,\
                             REAL value_split);
};
Tree::Tree(int criterion,uint n_classes,uint n_features,uint max_features,uint min_sample_leaf,\
               uint max_depth,uint find_split_algorithm,\
           uint random_seed){
#ifdef DEBUG
    assert(max_features<=n_features);
#endif
    max_features=MAX(1,max_features);
    if (criterion==CRITERION_MSE) {
        this->criterion=new MSE();
    }else{
        fprintf(stderr, "invalid criterion: %d\n",criterion);
    }
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
    this->max_depth=max_depth;
    this->max_features=max_features;
    this->min_sample_leaf=min_sample_leaf;
    this->random_seed=random_seed;
    this->numNodes=0;
}
Tree::~Tree(){
    for (uint i=0; i<max_nNodes; i++) {
        delete nodes[i];
    }
    delete []nodes;
}

void Tree::init_nodes(uint nSamples){
    numNodes=1;//ROOT node
    max_nNodes=2*(uint)((REAL)floor((REAL)(nSamples/(1>((int)min_sample_leaf-5)?1:(min_sample_leaf-5))))+1);
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
    //build decision tree
    
    uint i,l_node_id,r_node_id,work_node_id,feature_split,sample_ind_beg,sample_ind_end;
    REAL value_split,ini_error,min_error;
    uint *node_beg,*node_end,*sample_ind,*sample_ind_swap;
    
    //initialize nodes
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
    
//    cout<<"build tree: max_nNodes="<<max_nNodes<<endl;
    for (work_node_id=ROOT; work_node_id<max_nNodes; work_node_id++) {
        if (work_node_id>numNodes) {
            break;
        }
        sample_ind_beg=node_beg[work_node_id];
        sample_ind_end=node_end[work_node_id];
//        cout<<"work node="<<work_node_id<<endl;
        //initialize criterion
        criterion->init(y, sample_ind,sample_ind_beg,sample_ind_end, nodes[work_node_id]->nSamples);
        ini_error=criterion->eval();
        nodes[work_node_id]->ini_error=ini_error;
        //terminal node if the number of sample in a node is less than min_sample_leaf
        if (nodes[work_node_id]->nSamples<=min_sample_leaf || nodes_deep[work_node_id]==max_depth \
            || numNodes>=max_nNodes-2) {
            nodes[work_node_id]->leaf=true;
            //allocate memory for prediction
            nodes[work_node_id]->pred=new REAL[n_classes];
            criterion->estimate(nodes[work_node_id]->pred);
            cout<<"leaf node:"<<work_node_id<<" pred="<<nodes[work_node_id]->pred[0]<< endl;
            continue;
        }
        //find split
        bool flag=find_split(work_node_id, X, y,sample_ind, sample_ind_beg,sample_ind_end, nodes[work_node_id]->nSamples, feature_split, value_split,min_error);
//            cout<<"split flag min_error "<<min_error<<endl;
        //terminal node
        if (!flag) {
            nodes[work_node_id]->leaf=true;
//                cout<<"leaf node:"<<work_node_id<<endl;
//                criterion->init(y, mask, nodes[work_node_id]->nSamples, total_nSamples);
            criterion->reset();
            //allocate memory for prediction
            nodes[work_node_id]->pred=new REAL[n_classes];
            criterion->estimate(nodes[work_node_id]->pred);
            cout<<"leaf node:"<<work_node_id<<" pred="<<nodes[work_node_id]->pred[0]<< endl;
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
        
        fprintf(stderr, "deep=%d,(%d,%f,%f,pred=%f)\n",nodes_deep[work_node_id],feature_split+1,value_split,min_error,0.0);
        
        //map samples in parent node to left,right child
        uint mid=updateNodeSampleMap(work_node_id, l_node_id,r_node_id, sample_ind,sample_ind_swap, X, sample_ind_beg, sample_ind_end, feature_split, value_split);
        node_beg[l_node_id]=sample_ind_beg;
        node_end[l_node_id]=sample_ind_beg + mid;
        node_beg[r_node_id]=sample_ind_beg + mid;
        node_end[r_node_id]=sample_ind_end;
        //        cout<<"left_nSamples["<<l_node_id<<"]="<<nodes[l_node_id]->nSamples<<" right_nSamples["<<r_node_id<<"]="<<nodes[r_node_id]->nSamples<<endl;
//        cout<<"left_beg="<<node_beg[l_node_id]<<"left_end="<<node_end[l_node_id]<<endl;
//        cout<<"right_beg="<<node_beg[r_node_id]<<"right_end="<<node_end[r_node_id]<<endl;
#ifdef PRINT_DEBUG
        fprintf(stderr, "deep=%d,node_id=%d,feature_split=%d,value_split=%f,min error=%f\n",\
                deep,work_node_id,feature_split,value_split,min_error);
#endif
    }
//    delete []mask;
    //release NodeSampleMap;
    delete []sample_ind_swap;
    delete []sample_ind;
    delete []node_beg;
    delete []node_end;
    delete []nodes_deep;//no more need nodes_deep
    delete criterion; //no more need Criterion
    return SUCCESS;
}

bool Tree::find_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error){
    if (find_split_algorithm==FIND_BEST) {
        return find_best_split(node_id, X,y, sample_ind, s_ind_beg,s_ind_end,nSamples, feature_split, value_split,min_error);
    }
    return -1;
}

bool Tree::find_best_split(uint node_id,REAL** X,REAL* y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples,uint& feature_split,REAL& value_split,REAL& min_error){
    
#ifdef DEBUG
    assert(s_ind_beg>=0 && s_ind_beg<s_ind_end);
    assert((s_ind_end-s_ind_beg) <=nSamples);
#endif
    
    uint i,j,node_nSamples,f,count;
    bool* skip=new bool[n_features];
    uint* skip_idx=new uint[n_features];
    min_error=HUGE_VAL;
    node_nSamples=nodes[node_id]->nSamples;
    MyPair* x_f=new MyPair[nSamples];
    
    //generate random split feature
    for (i=0; i<n_features; i++) {
        skip[i]=true;
        skip_idx[i]=i;
    }
    count=n_features;
//    cout<<"======================"<<endl;
    for (i=0; i<max_features; i++) {
//        j=rand()%count;
        j=randomMT()%count;
        uint idx=skip_idx[j];
        skip[idx]=false;
        std::swap(skip_idx[j], skip_idx[--count]);
//        cout<<"skip="<<idx<<endl;
    }
//    criterion->init(y, mask, node_nSamples, total_nSamples);
    //main loop
//    cout<<"find_best_split"<<endl;
    for (f=0; f<n_features; f++) {
//        cout<<"feature="<<f<<endl;
        if (skip[f]) {
            continue;
        }
        criterion->reset();
//        criterion->init(y, mask, node_nSamples, total_nSamples);
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
//            cout<<"nLeft="<<nLeft<<endl;
            //do not consider split if feature value are identical.
            if (x_f[count].value==x_f[count+1].value) {
                count++;
                continue;
            }
            REAL error=criterion->eval();
            cout<<"f="<<f<<' '<<error<<endl;
            if (error<min_error) {
                fprintf(stderr, "f=%d,pre_error=%.8f,now_error=%.8f\n",f,min_error,error);
                min_error=error;
                feature_split=f;
                value_split=0.5*(x_f[count].value+x_f[count+1].value);
            }
//            cout<<"error="<<error<<endl;
            count++;
        }
//        fprintf(stderr, "node_id=%d,nSamples=%d,f=%d,min_error=%.8lf\n",node_id,node_nSamples,f,min_error);
    }
    delete []skip;
    delete []skip_idx;
    delete []x_f;
    return min_error!=HUGE_VAL;
}
void Tree::predict(REAL **X,REAL* pred,uint nSamples,uint nFeatures){
#ifdef DEBUG
    assert(nFeatures==this->n_features);
#endif
    uint node_id,f_s;
    REAL v_s;
    for (uint i=0; i<nSamples; i++) {
        node_id=ROOT;
//        fprintf(stderr, "sample %d:",i);
//        for (uint j=0; j<n_features; j++) {
//            cout<<"X["<<j<<"]="<<X[j][i]<<endl;
//        }
        while (!nodes[node_id]->leaf) {
//            cout<<node_id<<' '<<nodes[node_id]->leaf<<endl;
//            fprintf(stderr, "(%d,%d,%f)\n",node_id,nodes[node_id]->feature_split,\
                    nodes[node_id]->value_split);
            f_s=nodes[node_id]->feature_split;
            v_s=nodes[node_id]->value_split;
            if (X[i][f_s]<=v_s) {
//                cout<<X[f_s][i]<<":"<<v_s<<" left"<<endl;
//                node_id=LEFT(node_id);
                node_id=nodes[node_id]->left_child;
            }else {
//                cout<<X[f_s][i]<<":"<<v_s<<" right"<<endl;
//                node_id=RIGHT(node_id);
                node_id=nodes[node_id]->right_child;
            }
        }
//        cout<<node_id<<" pred="<<nodes[node_id]->pred[0]<<endl;
        for (uint j=0; j<n_classes; j++) {
            pred[i]=nodes[node_id]->pred[j];
        }
    }
}
#endif
