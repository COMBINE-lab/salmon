/* * * * *
 *  ClassificationCriterion.h 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_ClassificationCriterion_h
#define libTM_ClassificationCriterion_h
#include "Criterion.h"
#include "cmath"
class ClassificationCriterion:public Criterion{
    //Abstract criterion for classification.
public:
    ClassificationCriterion(uint n_classes);
    virtual ~ClassificationCriterion();
    virtual void init(REAL *y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples);
    virtual void estimate(REAL* pred);
    virtual void reset();
    virtual REAL eval()=0;
    virtual int update (REAL *y,uint* loc,uint a,uint b);
    virtual int update_next(REAL *y,uint a);
public:
    uint *label_count_init;
    uint *label_count_left;
    uint *label_count_right;
};
ClassificationCriterion::ClassificationCriterion(uint n_classes){
#ifdef DEBUG
    assert(n_classes>=2);
#endif
    this->n_classes=n_classes;
    label_count_init=new uint[n_classes];
    label_count_left=new uint[n_classes];
    label_count_right=new uint[n_classes];
}
ClassificationCriterion::~ClassificationCriterion(){
    delete []label_count_left;
    delete []label_count_right;
    delete []label_count_init;
}
void ClassificationCriterion::init(REAL *y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples){
    /* Initialise the criterion class; assume all samples
     * are in the right branch.
     * put all samples y[idx] which idx belong to (sample_ind[s_ind_beg],sample_ind[s_ind_end]) 
     * to right right branch 
     * sample_ind stored the real index in y
     * y[i] must be non-negative and less than n_classes
     */
    this->nSamples=nSamples;
    this->nLeft=0;
    for (uint i=0; i<n_classes; i++) {
        label_count_left[i]=0;
        label_count_right[i]=0;
        label_count_init[i]=0;
    }
#ifdef DEBUG
    uint count=0;
#endif
    for (uint i=s_ind_beg; i<s_ind_end; i++) {
        uint idx=sample_ind[i];
#ifdef DEBUG
            assert(y[i]>=0 && y[i]<n_classes);
            count++;
#endif
        uint category=(uint)y[idx];
        label_count_init[category]++;
    }
#ifdef DEBUG
    assert(count==nSamples);
#endif
    reset();
}
void ClassificationCriterion::reset(){
    nLeft=0;
    for (uint i=0; i<n_classes; i++) {
        label_count_left[i]=0;
        label_count_right[i]=label_count_init[i];
    }
}
void ClassificationCriterion::estimate(REAL* pred){
    for (uint i=0; i<n_classes; i++) {
        pred[i]=label_count_init[i];
    }
}

//==================FIX ME===========================
int ClassificationCriterion::update(REAL *y, uint *loc, uint a, uint b){
    /*Update the criteria for each value in interval [a,b) and [b,end],
     [a,b) in left side,[b,end] in right side(where a and b
     are indices in loc).
     */
#ifdef DEBUG
    assert(a>=0);
    assert(b<=nSamples);
#endif
    for (uint i=a; i<b; i++) {
        uint idx=loc[i];
        uint category=y[idx];
        label_count_left[category]++;
        label_count_right[category]--;
        nLeft++;
    }
    return nLeft;

}
//===================================================


int ClassificationCriterion::update_next(REAL *y, uint a){
    /* put sample[a] into left side and update the criteria*/
#ifdef DEBUG
    assert(a>=0);
#endif
    uint category=y[a];
#ifdef DEBUG
    assert(category>=0 && category<n_classes);
#endif
    label_count_left[category]++;
    label_count_right[category]--;
    nLeft++;
    return nLeft;

}

class Entropy:public ClassificationCriterion {
    //Cross Entropy splitting criteria.
    
public:
    Entropy(uint n_classes);
    
    REAL eval();
};

Entropy::Entropy(uint n_classes):ClassificationCriterion(n_classes){
    
}

REAL Entropy::eval(){
    REAL l_entropy,r_entropy;
    l_entropy=r_entropy=0;
    uint nRight=nSamples-nLeft;
    for (uint i=0; i<n_classes; i++) {
        if (nLeft>0 && label_count_left[i]>0) {
            l_entropy-=(1.0*label_count_left[i]/nLeft)*log((REAL)1.0*label_count_left[i]/nLeft);
        }
        if (nRight>0 && label_count_right[i]>0) {
            r_entropy-=(1.0*label_count_right[i]/nRight)*log((REAL)1.0*label_count_right[i]/nRight);
        }
    }
    l_entropy*=1.0*nLeft/nSamples;
    r_entropy*=1.0*nRight/nSamples;
    return l_entropy+r_entropy;
}

class Gini:public ClassificationCriterion {
    //Gini spliting criteria
    
public:
    Gini(uint n_classes);
    REAL eval();
};

Gini::Gini(uint n_classes):ClassificationCriterion(n_classes){
    
}
REAL Gini::eval(){
    REAL l_score,r_score;
    uint nRight=nSamples-nLeft;
    l_score=(REAL)nLeft;
    r_score=(REAL)nRight;
    for (uint i=0; i<n_classes; i++) {
        //Caution:avoid overflow.
        if (label_count_left[i]>0) {
            l_score-=(REAL)label_count_left[i]/nLeft*label_count_left[i];
        }
        if (label_count_right[i]>0) {
            r_score-=(REAL)label_count_right[i]/nRight*label_count_right[i];
        }
    }
#ifdef DEBUG
    assert(l_score>=0);
    assert(r_score>=0);
#endif
//    if (nLeft!=0) {
//        l_score/=nLeft;
//    }
//    if (nRight!=0) {
//        r_score/=nRight;
//    }    
    l_score/=nSamples;
    r_score/=nSamples;
    return l_score+r_score;
}

#endif
