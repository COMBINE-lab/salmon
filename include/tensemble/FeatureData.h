/* * * * *
 *  FeatureData.h 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_FeatureData_h
#define libTM_FeatureData_h
#include <iostream>
#include <algorithm>
#include "TypeDef.h"
using namespace std;

class FeatureData { // represents a training data set distributed among processors feature-wise
public:
    // constructor/destructor
//    FeatureData();
    FeatureData(uint nSamples, uint nFeaures);
    ~FeatureData();
    
	// reading and initialization
    void sort();
	
	// manage sample_in_node
	void reset();
	
    // queries
    uint getNumSamples();
    uint getNumFeatures();
    int getSampleInNode(uint i);//get node of sample[i]
    uint_vec getNodeSamples(uint i);//get all sample from node i
    void setNode(uint i_sample, uint i_node);//set i_sample to i_node

    REAL getFeature(uint f, uint i);
    REAL getSortedFeature(uint f, uint i);
    uint getSortedIndex(uint f, uint i);
    
    uint whoHasFeature(uint f);
    bool isLocalFeature(uint f);
    uint localFeatureIndex(uint gf);
    uint globalFeatureIndex(uint lf);
    
    // prediction
    void updatePred(uint i, double p);
    void updateResiduals();
    
    
public:
	// dataset descriptors
	uint nSamples; // number of data instances
    uint nFeatures; // number of features stored on this processor
    
	// static attributes
    REAL** rawfeatures; // feature values ordered by instance
    REAL** sortedfeatures; // feature values ordered by value
    uint** sortedindices; // indices of original instances for each sorted feature
    REAL* label; // target label value of each instance
    
	// level-specific attributes
	int* sample_in_node; // indicated sample in which nodes.
	
};

//FeatureData::FeatureData(){
//    this->nSamples=this->nFeatures=-1;
//    rawfeatures=NULL;
//    sortedfeatures=NULL;
//    sortedindices=NULL;
//    label=NULL;
//    sample_in_node=NULL;
//}

FeatureData::FeatureData(uint nSamples, uint nFeaures) {
    this->nSamples=nSamples;
    this->nFeatures = nFeaures;

    // rawfeatures: initialized to minimum value (for missing values)
    rawfeatures = new REAL*[nSamples];
	for (int i=0; i<nSamples; i++) {
        rawfeatures[i] = new REAL[nFeaures];
		for (int j=0; j<nFeaures; j++)
            rawfeatures[i][j] = 0.f; //-9999999.f; // TODO replace with better minimum value -- min float?
	}
	
	// sortedfeatures, sortedindices: limited init, completed during sort()
    sortedfeatures = new REAL*[nFeaures];
    sortedindices = new uint*[nFeaures];
    for (int i=0; i<nFeaures; i++) {
        sortedfeatures[i] = new REAL[nSamples];
        sortedindices[i] = new uint[nSamples];
	}
	
	// label: limited init, read from file
    label = new REAL[nSamples];
	
	// node: initialized to 0
    sample_in_node = new int[nSamples];
    for (int i=0; i<nSamples; i++)
		sample_in_node[i] = -1;

}

FeatureData::~FeatureData() {
    // delete all 1-d arrays: qid, label, node, pred, residual, idealdcg
    if (label) {
        delete [] label;
    }
    if (sample_in_node) {
        delete [] sample_in_node;
    }
    label=NULL;
    sample_in_node=NULL;
	
	// delete all 2-d arrays: rawfeatures, sortedfeatures, sortedindices
    for (int i=0; i<nSamples; i++) {
        if(rawfeatures[i])
            delete [] rawfeatures[i];
        rawfeatures[i] = NULL;
    }
	for (int i=0; i<nFeatures; i++) {
        if(sortedfeatures[i])
            delete [] sortedfeatures[i];
        sortedfeatures[i] = NULL;
        if(sortedindices[i])
            delete [] sortedindices[i];
        sortedindices[i] = NULL;
    }
    delete[] rawfeatures;
    delete[] sortedfeatures;
    delete[] sortedindices;
}


void FeatureData::reset() {
    // clear nodes before next tree
	for (int i=0; i<nSamples; i++) {
        sample_in_node[i]=-1;
	}
}

class FeatureValuePair {
public:
    uint index;
	REAL value;
};

struct CompareFeatureValuePairs {
    bool operator() (FeatureValuePair* fv1, FeatureValuePair* fv2) {
        return (fv1->value < fv2->value);
    }
};

void FeatureData::sort() {
    // initialize FeatureValue array
    CompareFeatureValuePairs cfvp;
    FeatureValuePair** pairs = new FeatureValuePair*[nSamples];
    for (int i=0; i<nSamples; i++)
        pairs[i] = new FeatureValuePair();
    
    // sort each feature
	for (int f=0; f<nFeatures; f++) {
	    // load feature into pairs array
        for (int i=0; i<nSamples; i++) {
            pairs[i]->index = i;
            pairs[i]->value = rawfeatures[f][i];
        }
	    
	    // sort pairs array
		std::sort(pairs, pairs + nSamples, cfvp);
		
		// write feature to sortedfeatures, sortedindices
        for (int i=0; i<nSamples; i++) {
            sortedfeatures[f][i] = pairs[i]->value;
            sortedindices[f][i] = pairs[i]->index;
        }
	}
	
	// delete FeatureValue array
    for (int i=0; i<nSamples; i++) {
        delete pairs[i];
        pairs[i] = NULL;
    }
    delete [] pairs;
}

uint FeatureData::getNumSamples() {
    return nSamples;
}

int FeatureData::getSampleInNode(uint i) {
    return sample_in_node[i];
}

uint_vec FeatureData::getNodeSamples(uint i){
    uint_vec ret;
    for (uint k=0; k<nSamples; k++) {
        if (sample_in_node[k]==i) {
            ret.push_back(k);
        }
    }
    return ret;
}

void FeatureData::setNode(uint i_sample, uint i_node) {
    sample_in_node[i_sample]=i_node;
}


REAL FeatureData::getFeature(uint f, uint i) {
    return rawfeatures[f][i];
}

REAL FeatureData::getSortedFeature(uint f, uint i) {
    return sortedfeatures[f][i];
}

uint FeatureData::getSortedIndex(uint f, uint i) {
    return sortedindices[f][i];
}

#endif
