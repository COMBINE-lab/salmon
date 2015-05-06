/* * * * *
 *  TypeDef.h 
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Copyright (C) 2012, Rongkai Xia, shinekaixia@gmail.com
 *
 * * * * */

#ifndef libTM_TypeDef_h
#define libTM_TypeDef_h
#include <vector>
#include "cstdio"
#define ENSEMBLE_SUCCESS 0
#define ENSEMBLE_FAIL    -1

#define EPS 1e-18
#define LEFT(A) ((A)<<1)
#define RIGHT(A) (((A)<<1)+1)

#define MAX_DEPTH   100
#define ROOT    1
#define THREAD_MIN_FEATURES 5
#define MTRY_DEFAULT    0

typedef double REAL;
typedef unsigned int uint;
typedef short int sint;
typedef std::vector<uint> uint_vec;


#define MAX(A,B)    (((A)>(B))?(A):(B))
#define MIN(A,B)    (((A)>(B))?(B):(A))


#endif
