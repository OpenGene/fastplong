#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "read.h"
#include "options.h"

using namespace std;

class Stats{
public:
    // this @guessedCycles parameter should be calculated using the first several records
    Stats(Options* opt, int guessedCycles = 0, int bufferMargin = 1024);
    ~Stats();
    int getCycles();
    long getReads();
    long getBases();
    long getQ20();
    long getQ30();
    long getGCNumber();
    // by default the qualified qual score is Q20 ('5')
    void statRead(Read* r);

    static Stats* merge(vector<Stats*>& list);
    void print();
    void summarize(bool forced = false);
    // a port of JSON report
    void reportJson(ofstream& ofs, string padding);
    // a port of HTML report
    void reporHtmlMedianQualLengthDensity(ofstream& ofs, string filteringType);
    void reporHtmlMedianQualHist(ofstream& ofs, string filteringType);
    void reportHtmlBasicInfo(ofstream& ofs, string filteringType);
    void reportHtmlQuality(ofstream& ofs, string filteringType);
    void reportHtmlContents(ofstream& ofs, string filteringType);
    void reportHtmlKMER(ofstream& ofs, string filteringType);
    bool isLongRead();
    int getMeanLength();
    void calcLengthHistogram();

public:
    static string list2string(double* list, long size, long* coords);
    template<class T>
    static string list2string(T* list, long size);
    static int base2val(char base);

private:
    void extendBuffer(int newBufLen);
    string makeKmerTD(int i, int j);
    string kmer3(int val);
    string kmer2(int val);

private:
    Options* mOptions;
    long mReads;
    int mEvaluatedSeqLen;
    /* 
    why we use 8 here?
    map A/T/C/G/N to 0~7 by their ASCII % 8:
    'A' % 8 = 1
    'T' % 8 = 4
    'C' % 8 = 3
    'G' % 8 = 7
    'U' % 8 = 5
    'N' % 8 = 6
    */
    long *mCycleQ30Bases[8];
    long *mCycleQ20Bases[8];
    long *mCycleBaseContents[8];
    long *mCycleBaseQual[8];
    long *mCycleTotalBase;
    long *mCycleTotalQual;
    long *mKmer;
    long mBaseQualHistogram[128];
    long mMedianReadQualHistogram[128];
    long mMedianReadQualBases[128];

    map<string, double*> mQualityCurves;
    map<string, double*> mContentCurves;
    map<string, long> mOverRepSeq;
    map<string, long*> mOverRepSeqDist;

    //qual-length distribution
    map<char, vector<int>> mQualLength;

    int mCycles;
    int mBufLen;
    long mBases;
    long mQ20Bases[8];
    long mQ30Bases[8];
    long mBaseContents[8];
    long mQ5Total;
    long mQ7Total;
    long mQ10Total;
    long mQ15Total;
    long mQ20Total;
    long mQ30Total;
    long mQ40Total;
    bool summarized;
    long mKmerMax;
    long mKmerMin;
    int mKmerBufLen;
    long mLengthSum;
    vector<int> mLengthVec;
    int mMinLen;
    int mMaxLen;
    int mMedianLen;
    int mN50Len;
    map<int, int> mLengthHist;
    bool mNeedCalcLength;
};

#endif