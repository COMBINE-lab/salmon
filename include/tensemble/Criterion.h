/* * * * *
 *  Criterion.h 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_Criterion_h
#define libTM_Criterion_h
#include "TypeDef.h"
#include "assert.h"
#define CRITERION_MSE   1
#define CRITERION_ENTROPY 2
#define CRITERION_GINI 3

class Criterion {
    //splitting criteria (regression and classification)
public:
//    Criterion(uint n_classes=1){
//    }
    virtual ~Criterion(){};
public:
    virtual void init(REAL *y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples)=0;
    virtual void estimate(REAL* pred)=0;
    virtual void reset()=0;
    virtual REAL eval()=0;
    virtual int update (REAL *y,uint* loc,uint a,uint b)=0;
    virtual int update_next(REAL *y,uint a)=0;
public:
    //interval [left:mid] is in the left side,[mid+1:end] in right side
    uint nSamples;
    uint nLeft;
    uint n_classes;
};


class RegressionCriterion:public Criterion {
    //Abstract criterion for regression
    
public:
    virtual void init(REAL *y,uint* sample_ind,uint s_ind_beg,uint s_ind_end,uint nSamples);
    virtual void estimate(REAL* pred);
    virtual void reset();
    virtual REAL eval()=0;
    virtual int update (REAL *y,uint* loc,uint a,uint b);
    virtual int update_next(REAL *y,uint a);
public:
    REAL mean_all;
    REAL mean_left;
    REAL mean_right;
    REAL sq_sum_all;
    REAL sq_sum_left;
    REAL sq_sum_right;
    REAL var_left;
    REAL var_right;
};
void RegressionCriterion::init(REAL *y,uint* sample_ind,\
                               uint s_ind_beg,uint s_ind_end,uint nSamples){
    /* Initialise the criterion class; assume all samples
     * are in the right branch.
     * put all samples y[idx] which idx belong to (sample_ind[s_ind_beg],sample_ind[s_ind_end]) 
     * to right right branch 
     * sample_ind stored the real index in y
     */
    this->n_classes=1;
    this->nSamples=nSamples;
    this->nLeft=0;
    mean_all=0.0;
    mean_left=0.0;
    mean_right=0.0;
    sq_sum_all=0.0;
    sq_sum_left=0.0;
    sq_sum_right=0.0;
    var_left=0.0;
    var_right=0.0;
#ifdef DEBUG
    uint count=0;
#endif
    for (uint i=s_ind_beg; i<s_ind_end; i++) {
        uint idx=sample_ind[i];
        mean_all+=y[idx];
        sq_sum_all+=y[idx]*y[idx];
#ifdef DEBUG
        count++;
#endif
    }
    
#ifdef DEBUG
    assert(count==nSamples);
#endif
    mean_all/=nSamples;
    reset();
}
void RegressionCriterion::reset(){
    nLeft=0;
    mean_left=0.0;
    mean_right=mean_all;
    sq_sum_left=0.0;
    sq_sum_right=sq_sum_all;
    var_left=0.0;
    var_right=sq_sum_all-nSamples*(mean_right*mean_right);
}
void RegressionCriterion::estimate(REAL* pred){
    //regression only need pred[0]
    pred[0]=this->mean_all;
}

//=====================FIX ME==========================================
int RegressionCriterion::update(REAL* y,uint* loc,uint a,uint b){
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
        sq_sum_left+=y[idx]*y[idx];
        sq_sum_right-=y[idx]*y[idx];
        mean_left=(nLeft*mean_left+y[idx])/(nLeft+1);
        mean_right=((nSamples-nLeft)*mean_right-y[idx])/(nSamples-nLeft-1);
        nLeft++;
    }
    return nLeft;
}
//======================================================================

int RegressionCriterion::update_next(REAL* y,uint a){
     /* put sample y[a] into left side and update the criteria*/
#ifdef DEBUG
    assert(a>=0);
#endif
    sq_sum_left+=y[a]*y[a];
    sq_sum_right-=y[a]*y[a];
    mean_left=(nLeft*mean_left+y[a])/(nLeft+1);
    if (nLeft+1==nSamples) {
        mean_right=0.0;
    }else{
        mean_right=(1.0*(nSamples-nLeft)*mean_right-y[a])/(nSamples-nLeft-1);
    }
    nLeft++;
    var_left=sq_sum_left-nLeft*(mean_left*mean_left);
    var_right=sq_sum_right-(nSamples-nLeft)*(mean_right*mean_right);
    return nLeft;
}

class MSE:public RegressionCriterion {
    /*Mean squared error impurity criterion.*/
    
public:
    REAL eval();
};
REAL MSE::eval(){
    return this->var_left+this->var_right;
}
#endif
