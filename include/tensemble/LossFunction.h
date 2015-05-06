/* * * * *
 *  LossFunction.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_LossFunction_h
#define libTM_LossFunction_h
#include "Estimator.h"
#include "Tree.h"

class LossFunction {
    /* ABSTRACT LOSS FUNCTION
     * parameter n_classes and k in the following function is for classification
     * get_init_estimate: return the prior estimate(initialize predction)
     * loss: calculate the loss
     * compute_residual: return the residual(or pseudoresponse).For multi-class problems,return the
     * k-th residual(or pseudoresponse)
     */
    
public:
    REAL alpha;//just for QuantileLossFunction and HuberLossFunction
public:
    virtual ~LossFunction(){};
    
    virtual void get_init_estimate(REAL* y,REAL* pred,REAL* prior,uint nSamples,uint n_classes=1)=0;
    
    virtual REAL loss(REAL* y,REAL* pred,bool* mask,bool invert_mask,uint nSamples,uint n_classes=1)=0;
    
    virtual void compute_residual(REAL* residual, REAL* y,REAL* pred,uint nSamples,uint n_classes=1,uint k=0)=0;
    
    virtual void update(Tree* tree, REAL** X,bool* mask,REAL* y,REAL* pred,REAL* residual,\
                        uint nSamples,uint nFeatures,REAL learning_rate,uint n_classes=1,uint k=0);
    
    virtual void update_terminal_region(Tree* tree, uint* region, bool* mask,REAL* y, \
                                        REAL* residual,uint nSamples,uint n_classes)=0;
};

void LossFunction::update(Tree* tree, REAL** X,bool* mask,REAL* y,REAL* pred,REAL* residual,\
                          uint nSamples,uint nFeature,REAL learning_rate,uint n_classes,uint k){
    /* UPDATE FUNCTION FOR GRADIENT BOOSTING
     * update k-th tree for gradient boosting.
     * k=0,n_classes=1 for regression and two-classes classification
     */
    uint *region=new uint[nSamples];
    tree->predict_terminal_region(X, region, nSamples, nFeature);
    update_terminal_region(tree, region, mask, y, residual, nSamples,n_classes);
    for (uint i=0; i<nSamples; i++) {
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        uint node_id=region[i];
        pred[i]+=learning_rate*tree->nodes[node_id]->pred[0];
//#ifdef DEBUG
//        assert(pred[i]>=0);
//#endif
    }
    delete []region;
}

class SquareLoss:public LossFunction {
    //Square loss
    
public:
    void get_init_estimate(REAL* y,REAL* pred,REAL* prior,uint nSamples,uint n_classes);
    
    REAL loss(REAL* y, REAL* pred,bool* mask,bool invert_mask,uint nSamples,uint n_classes);
    
    void compute_residual(REAL* residual, REAL* y,REAL* pred,uint nSamples,uint n_classes,uint k);
    
    void update(Tree* tree, REAL** X,bool* mask,REAL* y,REAL* pred,REAL* residual,\
                        uint nSamples,uint nFeatures,REAL learning_rate,uint n_classes,uint k);
    
    void update_terminal_region(Tree* tree, uint* region, bool* mask,REAL* y, \
                                        REAL* residual,uint nSamples,uint n_classes);
};

void SquareLoss::get_init_estimate(REAL *y, REAL *pred, REAL* prior, uint nSamples,uint n_classes){
    return MeanEstimator(y, pred, prior,nSamples);
}

REAL SquareLoss::loss(REAL *y, REAL *pred,bool* mask,bool invert_mask,uint nSamples,uint n_classes){
    /* compute the loss
     * if invert_mask=true,return train score,else return OOB score.
     */
    REAL ret=0.0;
    uint count=0;
    for (uint i=0; i<nSamples; i++) {
        if ( (mask[i] && !invert_mask) || (!mask[i] && invert_mask) ) {
            continue;
        }
        ret+=(y[i]-pred[i])*(y[i]-pred[i]);
        count++;
    }
    return ret/count;
}

void SquareLoss::compute_residual(REAL *residual, REAL *y, REAL *pred, uint nSamples,uint n_classes,uint k){
    /* compute the working residual, equals negative gradient. 
     */
//    uint count=0;
    for (uint i=0; i<nSamples; i++) {
//        if (!mask[i]) {
//            continue;
//        }
        residual[i]=y[i]-pred[i];
    }
}

void SquareLoss::update(Tree *tree, REAL **X, bool *mask, REAL *y, REAL *pred, REAL *residual, uint nSamples, uint nFeatures, REAL learning_rate,uint n_classes, uint k){
    /* Least squares does not need to update terminal regions.*/
#ifdef DEBUG
    assert(n_classes==1);
    assert(k==0);
#endif
    REAL* r_pred=new REAL[nSamples];
    tree->predict(X, r_pred, nSamples, nFeatures);
    for (uint i=0; i<nSamples; i++) {
        pred[i]+=learning_rate*r_pred[i];
    }
    delete []r_pred;
}

void SquareLoss::update_terminal_region(Tree *tree, uint *region, bool *mask, REAL *y, REAL *residual, uint nSamples,uint n_classes){
    /* Least squares does not need to update terminal regions.*/
    
}

class BinomialDeviance:public LossFunction {
    //negative binomial log-likelihood loss
    
public:
    void get_init_estimate(REAL* y,REAL* pred,REAL* prior,uint nSamples,uint n_classes);
    
    REAL loss(REAL* y, REAL* pred,bool* mask,bool invert_mask,uint nSamples,uint n_classes);
    
    void compute_residual(REAL* residual, REAL* y,REAL* pred,uint nSamples,uint n_classes,uint k);
    
    void update_terminal_region(Tree* tree, uint* region, bool* mask,REAL* y, \
                                REAL* residual,uint nSamples,uint n_classes);
};

void BinomialDeviance::get_init_estimate(REAL *y, REAL *pred, REAL* prior, uint nSamples,uint n_classes){
    return LogOddsEstimator(y, pred, prior, nSamples);
}

REAL BinomialDeviance::loss(REAL *y, REAL *pred,bool* mask,bool invert_mask,uint nSamples,uint n_classes){
    /* compute the loss
     * if invert_mask=true,return train score,else return OOB score.
     */
    REAL deviance=0;
    uint count=0;
    for (uint i=0; i<nSamples; i++) {
        if ( (mask[i] && !invert_mask) || (!mask[i] && invert_mask) ) {
            continue;
        }
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        deviance-=y[i]*pred[i]-log(1.0+exp(pred[i]));
        count++;
    }
    return deviance/count;
}

void BinomialDeviance::compute_residual(REAL *residual, REAL *y, REAL *pred, uint nSamples,uint n_classes,uint k){
    /* compute the working residual, equals negative gradient. 
     * y has zero-one values.e.g y[i]=1 if sample i belong positive class,otherwise y[i]=0.
     */

//    uint count=0;
    for (uint i=0; i<nSamples; i++) {
//        if (!mask[i]) {
//            continue;
//        }
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        residual[i]=y[i]-1.0/(1.0+exp(-pred[i]));
    }
}

void BinomialDeviance::update_terminal_region(Tree* tree, uint* region, bool* mask, REAL* y,\
                            REAL* residual,uint nSamples,uint n_classes){
    /* UPDATE TERMINAL REGION FOR BinomialDeviance
     * Make a single Newton-Raphson step
     * sample i is out of bag sample if mask[i]=false
     */
    uint i,k,numNodes=tree->numNodes;
    REAL *numerator,*denominator;
    numerator=new REAL[numNodes+1];
    denominator=new REAL[numNodes+1];
    for (i=ROOT; i<=numNodes; i++) {
        numerator[i]=denominator[i]=0;
    }
//    uint count=0;
    for (i=0; i<nSamples; i++) {
        if (!mask[i]) {
            continue;
        }
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        k=region[i];
        numerator[k]+=residual[i];
        denominator[k]+=(y[i]-residual[i])*(1.0-y[i]+residual[i]);
//        count++;
    }
    for (i=ROOT; i<=numNodes; i++) {
        //update terminal nodes
        if (tree->nodes[i]->leaf) {
            //using pred[0] to store the GBM prediction
            if (denominator[i]==0) {
                tree->nodes[i]->pred[0]=0;
            }else{
                tree->nodes[i]->pred[0]=numerator[i]/denominator[i];
            }
        }
    }
    delete []numerator;
    delete []denominator;
}

class MultinomialDeviance:public LossFunction {
    //multi-class logistic loss
    
public:
    void get_init_estimate(REAL* y,REAL* pred,REAL* prior,uint nSamples,uint n_classes);
    
    REAL loss(REAL* y, REAL* pred,bool* mask,bool invert_mask,uint nSamples,uint n_classes);
    
    void compute_residual(REAL* residual, REAL* y,REAL* pred,uint nSamples,uint n_classes,uint k);
    
    void update_terminal_region(Tree* tree, uint* region, bool* mask,REAL* y, \
                                REAL* residual,uint nSamples,uint n_classes);
};

void MultinomialDeviance::get_init_estimate(REAL *y, REAL *pred,REAL* prior, uint nSamples, uint n_classes){
    return MulticlassPriorEstimator(y, pred, prior, nSamples, n_classes);
}

REAL MultinomialDeviance::loss(REAL *y, REAL *pred,bool* mask,bool invert_mask, uint nSamples,uint n_classes){
    /* multinomial deviance loss
     * here y belong to [0,1,...,n_classes-1],not zero-one values
     * if invert_mask=true,return train score,else return OOB score.
     */
    
    uint *zero_one_y=new uint[nSamples*n_classes];
    uint i,j,k;
    REAL deviance=0.0,sum_exp;
    uint count=0;
    memset(zero_one_y, 0, nSamples*n_classes*sizeof(uint));
    for (i=0; i<nSamples; i++) {
#ifdef DEBUG
        assert(y[i]>=0 && y[i]<n_classes);
#endif
        k=y[i];
        zero_one_y[k*nSamples+i]=1;
    }
    
    for (i=0; i<nSamples; i++) {
        if ( (mask[i] && !invert_mask) || (!mask[i] && invert_mask) ) {
            continue;
        }
        sum_exp=0.0;
        for (j=0; j<n_classes; j++) {
            sum_exp+=exp(pred[j*nSamples+i]);
        }
#ifdef DEBUG
        assert(sum_exp!=0);
#endif
        for (j=0; j<n_classes; j++) {
            deviance-=((REAL)zero_one_y[j*nSamples+i]*pred[j*nSamples+i]-log(sum_exp));
        }
        count++;
    }
    delete []zero_one_y;
    return deviance/count;
}

void MultinomialDeviance::compute_residual(REAL *residual, REAL *y, REAL *pred, uint nSamples,uint n_classes, uint k){
    /* Compute residual(or pseudoresponse,negative gradient) for the k-th class. 
     * y has zero-one values.e.g y[i]=1 if sample i belong k-th class,otherwise y[i]=0.
     * pred has size nSamples*n_classes,that is pred=[pred_1,pred_2,...,pred_n_classes]
     * pred_k=[f_k(X_1),f_k(X_2),...,f_k(X_nSamples)]
     */
    uint i,j;
    REAL sum_exp,prob_k;
    for (i=0; i<nSamples; i++) {
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        sum_exp=0.0;
        for (j=0; j<n_classes; j++) {
            sum_exp+=exp(pred[j*nSamples+i]);
        }
        if (sum_exp==0.0) {
            prob_k=0.0;
        }else{
            prob_k=exp(pred[k*nSamples+i])/sum_exp;
        }
        residual[i]=y[i]-prob_k;
    }
}

void MultinomialDeviance::update_terminal_region(Tree *tree, uint *region, bool *mask, REAL *y, REAL *residual, uint nSamples,uint n_classes){
    /* UPDATE TERMINAL REGION FOR MultinomialDeviance
     * Make a single Newton-Raphson step
     * sample i is out of bag sample if mask[i]=false
     */
    uint i,k,numNodes=tree->numNodes;
    REAL *numerator,*denominator;
    numerator=new REAL[numNodes+1];
    denominator=new REAL[numNodes+1];
    for (i=ROOT; i<=numNodes; i++) {
        numerator[i]=denominator[i]=0.0;
    }
//    uint count=0;
    for (i=0; i<nSamples; i++) {
        if (!mask[i]) {
            continue;
        }
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        k=region[i];
        numerator[k]+=residual[i];
        denominator[k]+=(y[i]-residual[i])*(1-y[i]+residual[i]);
    }
    for (i=ROOT; i<=numNodes; i++) {
        //update terminal nodes
        if (tree->nodes[i]->leaf) {
            //using pred[0] to store the GBM prediction
            if (denominator[i]==0) {
                tree->nodes[i]->pred[0]=0.0;
            }else{
                tree->nodes[i]->pred[0]=((REAL)(n_classes-1.0)/n_classes)*numerator[i]/denominator[i];
            }
        }
    }
    delete []numerator;
    delete []denominator;
}
#endif
