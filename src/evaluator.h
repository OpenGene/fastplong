#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "util.h"
#include "read.h"

using namespace std;

class Evaluator{
public:
    Evaluator(Options* opt);
    ~Evaluator();
    // evaluate how many reads are stored in the input file
    void evaluateReadNum(long& readNum);
    void evalAdapterAndReadNum(Options* opt, long& readNum);
    void evaluateSeqLenAndCheckRNA();
    string extendKeyToAdapter(int key, unsigned int* counts, unsigned long* positionAcc, int keylen, bool isRNA = false, bool leftFirst = true);
    int getTopKey(unsigned int* counts, int keylen);

    static string matchKnownAdapter(string seq);
    string int2seq(unsigned int val, int seqlen, bool isRNA = false);
    int seq2int(string* seq, int pos, int seqlen, int lastVal = -1);
    int seq2int(string& seq, int pos, int seqlen, int lastVal = -1);
private:
    Options* mOptions;
    string getAdapterWithSeed(int seed, Read** loadedReads, long records, int keylen);
};


#endif
