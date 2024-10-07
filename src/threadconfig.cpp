#include "threadconfig.h"
#include "util.h"

ThreadConfig::ThreadConfig(Options* opt, int threadId, bool paired){
    mOptions = opt;
    mThreadId = threadId;
    mWorkingSplit = threadId;
    mCurrentSplitReads = 0;
    mPreStats1 = new Stats(mOptions);
    mPostStats1 = new Stats(mOptions);

    mWriter1 = NULL;

    mFilterResult = new FilterResult(opt, paired);
    mCanBeStopped = false;
    mInputList = NULL;
}

ThreadConfig::~ThreadConfig() {
    cleanup();
}

void ThreadConfig::cleanup() {
    if(mOptions->split.enabled && mOptions->split.byFileNumber)
        writeEmptyFilesForSplitting();
    deleteWriter();
    if(mInputList) {
        delete mInputList;
        mInputList = NULL;
    }
    if(mPreStats1) {
        delete mPreStats1;
        mPreStats1 = NULL;
    }
    if(mPostStats1) {
        delete mPostStats1;
        mPostStats1 = NULL;
    }
    if(mFilterResult) {
        delete mFilterResult;
        mFilterResult = NULL;
    }
}


void ThreadConfig::setInputList(SingleProducerSingleConsumerList<ReadPack*>* list) {
    mInputList = list;
}

void ThreadConfig::deleteWriter() {
    if(mWriter1 != NULL) {
        delete mWriter1;
        mWriter1 = NULL;
    }
}

void ThreadConfig::initWriter(string filename1) {
    deleteWriter();
    mWriter1 = new Writer(mOptions, filename1, mOptions->compression);
}

void ThreadConfig::initWriter(string filename1, string filename2) {
    deleteWriter();
    mWriter1 = new Writer(mOptions, filename1, mOptions->compression);
}

void ThreadConfig::addFilterResult(int result, int readNum) {
    mFilterResult->addFilterResult(result, readNum);
}


void ThreadConfig::initWriterForSplit() {
    if(mOptions->out.empty())
        return ;

    // use 1-based naming
    string num = to_string(mWorkingSplit + 1);
    // padding for digits like 0001
    if(mOptions->split.digits > 0){
        while(num.size() < mOptions->split.digits)
            num = "0" + num;
    }

    string filename1 = joinpath(dirname(mOptions->out), num + "." + basename(mOptions->out));
    initWriter(filename1);
}

void ThreadConfig::markProcessed(long readNum) {
    mCurrentSplitReads += readNum;
    if(!mOptions->split.enabled)
        return ;
    // if splitting is enabled, check whether current file is full
    if(mCurrentSplitReads >= mOptions->split.size) {
        // if it's splitting by file number, totally we cannot exceed split.number
        // if it's splitting by file lines, then we don't need to check
        if(mOptions->split.byFileLines || mWorkingSplit + mOptions->thread < mOptions->split.number ){
            mWorkingSplit += mOptions->thread;
            initWriterForSplit();
            mCurrentSplitReads = 0;
        } else {
            // this thread can be stoped now since all its tasks are done
            // only a part of threads have to deal with the remaining reads
            if(mOptions->split.number % mOptions->thread >0 
                && mThreadId >= mOptions->split.number % mOptions->thread)
                mCanBeStopped = true;
        }
    }
}

// if a task of writting N files is assigned to this thread, but the input file doesn't have so many reads to input
// write some empty files so it will not break following pipelines
void ThreadConfig::writeEmptyFilesForSplitting() {
    while(mWorkingSplit + mOptions->thread < mOptions->split.number) {
        mWorkingSplit += mOptions->thread;
            initWriterForSplit();
            mCurrentSplitReads = 0;
    }
}

bool ThreadConfig::canBeStopped() {
    return mCanBeStopped;
}