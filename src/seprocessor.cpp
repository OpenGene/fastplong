#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"

SingleEndProcessor::SingleEndProcessor(Options* opt){
    mOptions = opt;
    mReaderFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mLeftWriter =  NULL;
    mFailedWriter = NULL;

    mPackReadCounter = 0;
    mPackProcessedCounter = 0;

    mReadPool = new ReadPool(mOptions);
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
    if(mReadPool) {
        delete mReadPool;
        mReadPool = NULL;
    }
    delete[] mInputLists;
}

void SingleEndProcessor::initOutput() {
    if(!mOptions->failedOut.empty())
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);
    if(mOptions->out.empty() && !mOptions->outputToSTDOUT)
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out, mOptions->outputToSTDOUT);
}

void SingleEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mFailedWriter) {
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out.empty())
        return;

    if(mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}

bool SingleEndProcessor::process(){
    if(!mOptions->split.enabled)
        initOutput();

    mInputLists = new SingleProducerSingleConsumerList<ReadPack*>*[mOptions->thread];

    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        mInputLists[t] = new SingleProducerSingleConsumerList<ReadPack*>();
        configs[t] = new ThreadConfig(mOptions, t, false);
        configs[t]->setInputList(mInputLists[t]);
        initConfig(configs[t]);
    }

    std::thread readerThread(std::bind(&SingleEndProcessor::readerTask, this));

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::processorTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writerTask, this, mLeftWriter));
    if(mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writerTask, this, mFailedWriter));

    readerThread.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(!mOptions->split.enabled) {
        if(leftWriterThread)
            leftWriterThread->join();
        if(failedWriterThread)
            failedWriterThread->join();
    }

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    finalPreStats->calcLengthHistogram();
    Stats* finalPostStats = Stats::merge(postStats);
    finalPostStats->calcLengthHistogram();
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    // read filter results to the first thread's
    for(int t=1; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
    }

    cerr << "Before filtering:"<<endl;
    finalPreStats->print();
    cerr << endl;
    cerr << "After filtering:"<<endl;
    finalPostStats->print();

    cerr << endl;
    cerr << "Filtering result:"<<endl;
    finalFilterResult->print();


    // make JSON report
    JsonReporter jr(mOptions);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;
    if(failedWriterThread)
        delete failedWriterThread;

    if(!mOptions->split.enabled)
        closeOutput();

    return true;
}

void SingleEndProcessor::recycleToPool(int tid, Read* r) {
    // failed to recycle, then delete it
    if(!mReadPool->input(tid, r))
        delete r;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){
    string* outstr = new string();
    string* failedOut = new string();
    int tid = config->getThreadId();

    int readPassed = 0;
    for(int p=0;p<pack->count;p++){

        // original read1
        Read* or1 = pack->data[p];

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);

        int frontTrimmed = 0;
        // trim in head and tail, and apply quality cut in sliding window
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front, mOptions->trim.tail, frontTrimmed);

        if(r1 != NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        vector<Read*> outReads;

        if(r1 != NULL && mOptions->adapter.enabled){
            int trimmed = 0;
            if(!mOptions->adapter.sequenceStart.empty())
                trimmed += AdapterTrimmer::trimBySequenceStart(r1, config->getFilterResult(), mOptions->adapter.sequenceStart, mOptions->adapter.edMax, mOptions->adapter.trimmingExtension);
            if(!mOptions->adapter.sequenceEnd.empty())
                trimmed += AdapterTrimmer::trimBySequenceEnd(r1, config->getFilterResult(), mOptions->adapter.sequenceEnd,mOptions->adapter.edMax, mOptions->adapter.trimmingExtension);
            if(mOptions->adapter.hasFasta) {
                trimmed += AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, mOptions->adapter.edMax,mOptions->adapter.trimmingExtension);
            }
            if(trimmed > 0) {
                config->getFilterResult()->addReadTrimmed(trimmed);
            }

            //search for middle adapter
            int start = -1;
            int len = 0;
            bool foundMiddleAdapter = AdapterTrimmer::findMiddleAdapters(r1, mOptions->adapter.sequenceStart, mOptions->adapter.sequenceEnd, start, len, mOptions->adapter.edMax,mOptions->adapter.trimmingExtension);
            if(foundMiddleAdapter) {
                //break the read
                outReads = r1->breakByGap(start, len);
                //cerr << "break at " << start << ", " << len << ", read len: " << r1->length() << endl;
                //cerr << r1->mSeq->substr(start, len) << endl;
            } else {
                outReads.push_back(r1);
            }
        } else if(r1 != NULL) {
            outReads.push_back(r1);
        }

        //break by low quality regions
        if(mOptions->breakOpt.enabled && outReads.size() > 0) {
            vector<Read*> tmpReads;
            for(int i=0; i<outReads.size(); i++) {
                Read* rr = outReads[i];
                vector<pair<int, int>> regions = mFilter->detectLowQualityRegions(rr, mOptions->breakOpt.windowSize, mOptions->breakOpt.quality);
                if(regions.size() > 0) {
                    vector<Read*> brs = rr->breakByRegions(regions);
                    for(int j=0; j<brs.size(); j++)
                        tmpReads.push_back(brs[j]);
                    // release the read if it's created by breaking gap
                    if(rr != or1 && rr != r1)
                        recycleToPool(tid, rr);
                } else {
                    tmpReads.push_back(rr);
                }
            }
            outReads = tmpReads;
        }
        // mask by low quality regions
        if(mOptions->mask.enabled && outReads.size() > 0) {
            for(int i=0; i<outReads.size(); i++) {
                Read* rr = outReads[i];
                vector<pair<int, int>> regions = mFilter->detectLowQualityRegions(rr, mOptions->mask.windowSize, mOptions->mask.quality);
                for(int j=0; j<regions.size(); j++) {
                    rr->maskRegionWithN(regions[j].first, regions[j].second - regions[j].first + 1);
                }
            }
        }

        bool passed = false;
        for(int i=0; i<outReads.size(); i++) {

            Read* outr = outReads[i];
            int result = mFilter->passFilter(outr);

            config->addFilterResult(result, 1);

            if( outr != NULL &&  result == PASS_FILTER) {
                outr->appendToString(outstr);

                passed = true;
                // stats the read after filtering
                config->getPostStats1()->statRead(outr);
            } else if(mFailedWriter && outReads.size() == 1) {
                or1->appendToStringWithTag(failedOut, FAILED_TYPES[result]);
            }

            // release the read if it's created by breaking gap
            if(outr != or1 && outr!= r1 && outr != NULL)
                recycleToPool(tid, outr);
        }

        if(passed)
            readPassed++;

        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            recycleToPool(tid, r1);

        recycleToPool(tid, or1);
    }

    if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out.empty())
            config->getWriter1()->writeString(outstr);
    } 

    if(mLeftWriter) {
        mLeftWriter->input(tid, outstr);
        outstr = NULL;
    }
    if(mFailedWriter) {
        // write failed data
        mFailedWriter->input(tid, failedOut);
        failedOut = NULL;
    }

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    if(outstr)
        delete outstr;
    if(failedOut)
        delete failedOut;

    delete pack->data;
    delete pack;

    mPackProcessedCounter++;

    return true;
}

void SingleEndProcessor::readerTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader(mOptions->in, true);
    reader.setReadPool(mReadPool);
    int count=0;
    bool needToBreak = false;
    while(true){
        Read* read = reader.read();
        if(!read || needToBreak){
            // the last pack
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            mInputLists[mPackReadCounter % mOptions->thread]->produce(pack);
            mPackReadCounter++;
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M reads";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            mInputLists[mPackReadCounter % mOptions->thread]->produce(pack);
            mPackReadCounter++;
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
            // if the processor is far behind this reader, sleep and wait to limit memory usage
            while( mPackReadCounter - mPackProcessedCounter > PACK_IN_MEM_LIMIT){
                //cerr<<"sleep"<<endl;
                slept++;
                usleep(100);
            }
            readNum += count;
            // if the writer threads are far behind this reader, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while(mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) {
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
            // re-evaluate split size
            // TODO: following codes are commented since it may cause threading related conflicts in some systems
            /*if(mOptions->split.needEvaluation && !splitSizeReEvaluated && readNum >= mOptions->split.size) {
                splitSizeReEvaluated = true;
                // greater than the initial evaluation
                if(readNum >= 1024*1024) {
                    size_t bytesRead;
                    size_t bytesTotal;
                    reader.getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }*/
        }
    }

    for(int t=0; t<mOptions->thread; t++)
        mInputLists[t]->setProducerFinished();

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mReaderFinished = true;
    if(mOptions->verbose) {
        loginfo("Loading completed with " + to_string(mPackReadCounter) + " packs");
    }
    //lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete[] data;
}

void SingleEndProcessor::processorTask(ThreadConfig* config)
{
    SingleProducerSingleConsumerList<ReadPack*>* input = config->getInput();
    while(true) {
        if(config->canBeStopped()){
            break;
        }
        while(input->canBeConsumed()) {
            ReadPack* data = input->consume();
            processSingleEnd(data, config);
        }
        if(input->isProducerFinished()) {
            if(!input->canBeConsumed()) {
                if(mOptions->verbose) {
                    string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                    loginfo(msg);
                }
                break;
            }
        } else {
            usleep(100);
        }
    }
    input->setConsumerFinished();        

    mFinishedThreads++;
    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
    }

    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::writerTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
