/* * * * *
 *  RandomGenerator.h
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_RandomGenerator_h
#define libTM_RandomGenerator_h
//#include    "cokus.h"
#define MAX_UINT_COKUS 4294967295
#define uint32 unsigned long
#define SMALL_INT char

#ifdef MATLAB
#define SMALL_INT_CLASS mxCHAR_CLASS //will be used to allocate memory t
#endif


//this two function are implemented in cokus.cpp
void seedMT(uint32 seed);
uint32 randomMT(void);

REAL unif_rand(){
    return (((REAL)randomMT())/((REAL)MAX_UINT_COKUS));
}

#endif
