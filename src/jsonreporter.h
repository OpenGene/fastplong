#ifndef JSON_REPORTER_H
#define JSON_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "stats.h"
#include "filterresult.h"
#include <fstream>
#include <atomic>
#include "common.h"
#include "util.h"

using namespace std;

class JsonReporter{
public:
    JsonReporter(Options* opt);
    ~JsonReporter();

    void report(FilterResult* result, Stats* preStats1, Stats* postStats1);

private:
    Options* mOptions;
};


#endif
