/* * * * *
 *  ReadData.h 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_ReadData_h
#define libTM_ReadData_h
#include "FeatureData.h"
#include "fstream"
#include <cstdlib>
using namespace std;
#define PRINT_DEBUG

struct Data {
public:
    Data(){
        X=NULL;
        y=NULL;
        n_samples=n_features=-1;
    }
    ~Data(){
        if (X) {
            for (int i=0; i<n_samples; i++) {
                delete []X[i];
            }
            delete []X;
            X=NULL;
        }
        if (y) {
            delete []y;
            y=NULL;
        }
    }
    void set_size(int n_samples,int n_features){
        this->n_samples=n_samples;
        this->n_features=n_features;
        X=new REAL*[n_samples];
        y=new REAL[n_samples];
        for (int i=0; i<n_samples; i++) {
            X[i]=new REAL[n_features];
        }
    }
public:
    REAL** X;
    REAL* y;
    int n_samples;
    int n_features;
};

bool readData(Data& data, const char* x_file,const char* y_file=NULL){
    char* endptr,*idx,*val;
    ifstream x_input(x_file);
    ifstream y_input;
    if (y_file) {
        y_input.open(y_file);
    }
    int x_linenum,y_linenum;
    x_linenum=y_linenum=0;
    int n_samples=0;
    string strline;
	if (x_input.fail()) {
	    fprintf(stderr, "Error: unable to open data point file \"%s\".\n", x_file);
	    return false;
    }
    if (y_file && y_input.fail()) {
        fprintf(stderr, "Error: unable to open target file \"%s\".\n",y_file);
        return false;
    }
    //check file correctness
    getline(x_input,strline);
    while (!x_input.eof()) {
        x_linenum++;
        getline(x_input,strline);
    }
    if (y_file) {
        getline(y_input,strline);
        while (!y_input.eof()) {
            y_linenum++;
            getline(y_input,strline);
        }
        if (x_linenum!=y_linenum+1) {
            fprintf(stderr, "Error: number of samples in data point file(%d samples) and target file(%d samples) are not consistent.Please check your input file.\n",x_linenum-1,y_linenum );
            return false;
        }
        y_input.clear();
        y_input.seekg(0,ios_base::beg);
    }
    //reset
    x_input.clear();
    x_input.seekg(0,ios_base::beg);
    
    n_samples=x_linenum-1;
    getline(x_input, strline);
    int n_features=(int)strtol(strline.c_str(), &endptr, 10);
    if (endptr == strline.c_str() || (*endptr != '\0' && !isspace(*endptr))){
        fprintf(stderr, "Error: wrong input format at line %d in file %s\n",1,x_file);
        return false;
    }
    if (n_features<=0) {
        fprintf(stderr, "Error: number of features is 0 in training file.Please check your input file.\n");
        return false;
    }
    //read data
    data.set_size(n_samples, n_features);
    for (uint i=0; i<n_samples; i++) {
        if (y_file) {
            getline(y_input,strline);
            data.y[i]=strtod(strline.c_str(), &endptr);        
            if (endptr == strline.c_str() || (*endptr != '\0' && !isspace(*endptr))){
                fprintf(stderr, "Error: wrong input format at line %d in file %s\n",i+1,y_file);
                return false;
            }
        }
        
        getline(x_input, strline);
        char* line = strdup(strline.c_str()); 
//        cout<<line<<endl;
        idx = NULL;
        idx=strtok(line, ":");
        if (idx==NULL) {
            fprintf(stderr, "Error: Empty line at line %d in file %s\n",i+2,x_file);
            free(line);
            return false;
        }
        while (idx!=NULL) { // tok is feature:value
            val=strtok(NULL, " \t");
            if (val==NULL) {
                break;
            }
            int index=(int)strtol(idx, &endptr, 10);
            // check for errors
            if (endptr == idx || (*endptr != '\0' && !isspace(*endptr))){
                fprintf(stderr, "Error: wrong input format at line %d in file %s\n",i+2,x_file);
                free(line);
                return false;
            }
            if (index==0 || index>n_features) {
                fprintf(stderr, "Error: wrong input format at line %d in file %s.index must from 1 to n_features(number of features).Please check your input file.\n",i+2,x_file);
                free(line);
                return false;
            }
            REAL value=strtod(val, &endptr);
            if (endptr == val || (*endptr != '\0' && !isspace(*endptr))){
                fprintf(stderr, "Error: wrong input format at line %d in file %s\n",i+2,x_file);
                free(line);
                return false;
            }
            data.X[i][index-1]=value;
            idx=strtok(NULL, ":");
        }
        free(line);
    }
    return true;
}

#endif
