#ifndef THREAD_CONFIG_H
#define THREAD_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "stats.h"
#include "writer.h"
#include "options.h"
#include "filterresult.h"
#include "singleproducersingleconsumerlist.h"

using namespace std;

class ThreadConfig{
public:
    ThreadConfig(Options* opt, int threadId, bool paired = false);
    ~ThreadConfig();
    inline Stats* getPreStats1() {return mPreStats1;}
    inline Stats* getPostStats1() {return mPostStats1;}
    inline Writer* getWriter1() {return mWriter1;}
    inline FilterResult* getFilterResult() {return mFilterResult;}

    void initWriter(string filename1);
    void initWriter(string filename1, string filename2);

    void addFilterResult(int result, int readNum);

    int getThreadId() {return mThreadId;}
    // for splitting output
    // increase mCurrentSplitReads by readNum, and check it with options->split.size;
    void markProcessed(long readNum);
    void initWriterForSplit();
    bool canBeStopped();
    void cleanup();

    // input list
    void setInputList(SingleProducerSingleConsumerList<ReadPack*>* list);
    SingleProducerSingleConsumerList<ReadPack*>* getInput(){return mInputList;}

private:
    void deleteWriter();
    void writeEmptyFilesForSplitting();

private:
    Stats* mPreStats1;
    Stats* mPostStats1;
    Writer* mWriter1;
    Options* mOptions;
    FilterResult* mFilterResult;
    SingleProducerSingleConsumerList<ReadPack*>* mInputList;

    // for spliting output
    int mThreadId;
    int mWorkingSplit;
    long mCurrentSplitReads;
    bool mCanBeStopped;
};

#endif