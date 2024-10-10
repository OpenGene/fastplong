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

    static bool trimBySequenceStart(Read* r1, FilterResult* fr, string& adapter, double edMax = 0.2);
    static bool trimBySequenceEnd(Read* r1, FilterResult* fr, string& adapter, double edMax = 0.2);
    static bool trimByMultiSequences(Read* r1, FilterResult* fr, vector<string>& adapterList, bool incTrimmedCounter = true);
    static bool test();


};


#endif