/* * * * *
 *  TreeNode.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */


#ifndef libTM_TreeNode_h
#define libTM_TreeNode_h
#include    "TypeDef.h"

class TreeNode {
public:
    uint feature_split;
    REAL value_split;
    bool leaf;
    uint nSamples;
    uint left_child;
    uint right_child;
    REAL* pred;//only leaf node need allocate memory
    
    //for compute feature importance
    REAL ini_error;
    REAL best_error;
    TreeNode():feature_split(0),value_split(0.0),leaf(false),nSamples(0),left_child(0),right_child(0),pred(NULL),ini_error(0.0),best_error(0.0){
    }
    ~TreeNode(){
        if (pred) {
            delete pred;
            pred=NULL;
        }
    }
};

#endif
