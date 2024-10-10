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
    static int trimByMultiSequences(Read* r1, FilterResult* fr, vector<string>& adapterList);
    static bool test();


};


#endif