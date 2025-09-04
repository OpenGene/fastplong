#ifndef FILTER_H
#define FILTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "options.h"
#include "read.h"

using namespace std;

class Filter{
public:
    Filter(Options* opt);
    ~Filter();
    int passFilter(Read* r);
    bool passLowComplexityFilter(Read* r);
    Read* trimAndCut(Read* r, int front, int tail, int& frontTrimmed);
    //return all (start, end) pairs (0-based, inclusive) whose mean quality < quality in the read
    vector<pair<int, int>> detectLowQualityRegions(Read* r, int windowSize, int quality);

private:
    bool match(vector<string>& list, string target, int threshold);

private:
    Options* mOptions;
};


#endif
