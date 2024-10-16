#ifndef ADAPTER_TRIMMER_H
#define ADAPTER_TRIMMER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "filterresult.h"
#include "options.h"
#include "read.h"

using namespace std;

class AdapterTrimmer{
public:
    AdapterTrimmer();
    ~AdapterTrimmer();

    static int trimBySequenceStart(Read* r1, FilterResult* fr, string& adapter, double edMax = 0.3);
    static int trimBySequenceEnd(Read* r1, FilterResult* fr, string& adapter, double edMax = 0.3);
    static int trimByMultiSequences(Read* r1, FilterResult* fr, vector<string>& adapterList, double edMax = 0.3);
    static bool findMiddleAdapters(Read* r, string& startAdater, string& endAdapter, int& start, int& len, double edMax = 0.3);
    static int searchMiddleAdapter(string* read, string& adapter, double edMax = 0.3);
    static bool test();


};


#endif