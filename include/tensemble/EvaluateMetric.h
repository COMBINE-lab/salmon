/* * * * *
 *  EvaluateMetric.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_EvaluateMetric_h
#define libTM_EvaluateMetric_h
#include "TypeDef.h"

REAL mse(const REAL* preds, const REAL* label,uint nSamples) {
    REAL r = 0;
    int i, N = nSamples;
    for (i=0; i<N; i++)
        r += (label[i] - preds[i])*(label[i] - preds[i]); 
    return 1.0 / N * r;
}

REAL rmse(const REAL* preds, const REAL* label,uint nSamples) {
    REAL r = 0;
    int i, N = nSamples;
    for (i=0; i<N; i++)
        r += (label[i] - preds[i])*(label[i] - preds[i]); 
    return sqrt(1.0 / N * r);
}

REAL R2(const REAL* preds, const REAL* label,uint nSamples) {
    REAL r1 = 0,r2=0,r3=0;
    REAL mean_data=0,mean_preds=0;
    int i, N = nSamples;
    for (i=0; i<N; i++){
        mean_data+=label[i];
        mean_preds+=preds[i];
    }
    mean_data/=N;
    mean_preds/=N;
    for( i=0 ; i<N ; i++ )
    {
        r1+=(label[i]-mean_data)*(preds[i]-mean_preds);
        r2+=(label[i]-mean_data)*(label[i]-mean_data);
        r3+=(preds[i]-mean_preds)*(preds[i]-mean_preds);
    }
    if (r2==0.0 || r3==0.0) {
        if (r1==0.0) {
            return 1.0;
        }
        return 0.0;
    }
    return r1*r1/(r2*r3);
}

REAL BinomialDevianceLoss(REAL *y, REAL *pred,uint nSamples){
    /* binomial deviance loss
     * here y belong to [0,1]
     */    
    REAL deviance=0;
    uint count=0;
    for (uint i=0; i<nSamples; i++) {
#ifdef DEBUG
        assert(y[i]==0 || y[i]==1);
#endif
        deviance+=y[i]*pred[i]-log(1.0+exp(pred[i]));
        count++;
    }
    return deviance/count;
}

REAL Accuracy(REAL* y,REAL *pred,uint nSamples){
    REAL acc=0;
    for (uint i=0; i<nSamples; i++) {
        if (y[i]==pred[i]) {
            acc++;
        }
    }
    return acc/nSamples;
}
#endif
