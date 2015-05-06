/* * * * *
 *  Estimator.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_Estimator_h
#define libTM_Estimator_h
#include "TypeDef.h"
#include "assert.h"
#include <map>

void MeanEstimator(REAL* y,REAL* pred,REAL* prior,uint nSamples){
    REAL mean=0;
    for (uint i=0; i<nSamples; i++) {
        mean+=y[i];
    }
    mean/=nSamples;
//    mean=0;
    for (uint i=0; i<nSamples; i++) {
        pred[i]=mean;
    }
    *prior=mean;
}

void LogOddsEstimator(REAL* y,REAL* pred,REAL* prior,uint nSamples){
    /* for two class problem,y must be 0 or 1 */
    
    uint n_pos=0;
    for (uint i=0; i<nSamples; i++) {
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        n_pos+=y[i];
    }
    REAL prob=log((REAL)n_pos/(nSamples-n_pos));
    for (uint i=0; i<nSamples; i++) {
        pred[i]=prob;
    }
    *prior=prob;
}

void MulticlassPriorEstimator(REAL* y,REAL* pred,REAL* prior,uint nSamples,uint n_classes){
    /* for multi-class problems,y must be 0,1,...,n_classes-1 
     * pred is a nSamples*n_classes vector
     */
    
#ifdef DEBUG
    assert(n_classes>2);
#endif
    map<uint,uint> count;
    uint i,j,offset;
    for (i=0; i<n_classes; i++) {
        count[i]=0;
    }
    for (i=0; i<nSamples; i++) {
#ifdef DEBUG
        assert(y[i]>=0 && y[i]<n_classes);
#endif
        count[y[i]]++;
    }
    for (j=0; j<n_classes; j++) {
        offset=j*nSamples;
        for (i=0; i<nSamples; i++) {
            pred[offset+i]=count[j]/nSamples;
        }
    }
    //return prior
    for (j=0; j<n_classes; j++) {
        prior[j]=count[j]/nSamples;
    }
}
#endif
