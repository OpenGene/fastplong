#ifndef POLY_X_H
#define POLY_X_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "filterresult.h"
#include "options.h"
#include "read.h"

using namespace std;

class PolyX{
public:
    PolyX();
    ~PolyX();
    static void trimPolyX(Read* r1, FilterResult* fr, int compareReq);
    static bool test();


};


#endif