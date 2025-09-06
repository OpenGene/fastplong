#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>

using namespace std;

#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

class MaskOptions {
public:
    MaskOptions() {
        enabled = false;
        windowSize = 20;
        quality = 15;
    }
public:
    bool enabled;
    int windowSize;
    int quality;
};

class BreakOptions {
public:
    BreakOptions() {
        enabled = false;
        windowSize = 20;
        quality = 15;
    }
public:
    bool enabled;
    int windowSize;
    int quality;
};

class LowComplexityFilterOptions {
public:
    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
public:
    bool enabled;
    double threshold;
};


class PolyXTrimmerOptions {
public:
    PolyXTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class QualityCutOptions {
public:
    QualityCutOptions() {
        enabledFront = false;
        enabledTail = false;
        windowSizeShared = 4;
        qualityShared = 20;
        windowSizeFront = windowSizeShared;
        qualityFront = qualityShared;
        windowSizeTail = windowSizeShared;
        qualityTail = qualityShared;
    }
public:
    // enable 5' cutting by quality
    bool enabledFront;
    // enable 3' cutting by quality
    bool enabledTail;
    // the sliding window size
    int windowSizeShared;
    // the mean quality requirement
    int qualityShared;
    // the sliding window size for cutting by quality in 5'
    int windowSizeFront;
    // the mean quality requirement for cutting by quality in 5'
    int qualityFront;
    // the sliding window size for cutting by quality in 3'
    int windowSizeTail;
    // the mean quality requirement for cutting by quality in 3'
    int qualityTail;
};

class SplitOptions {
public:
    SplitOptions() {
        enabled = false;
        needEvaluation = false;
        number = 0;
        size = 0;
        digits = 4;
        byFileNumber = false;
        byFileLines = false;
    }
public:
    bool enabled;
    // number of files
    int number;
    // lines of each file
    long size;
    // digits number of file name prefix, for example 0001 means 4 digits
    int digits;
    // need evaluation?
    bool needEvaluation;
    bool byFileNumber;
    bool byFileLines;
};

class AdapterOptions {
public:
    AdapterOptions() {
        enabled = true;
        detected = false;
        trimmingExtension = 10;
        edMax = 0.25;
    }
public:
    bool enabled;
    string sequenceStart;
    string sequenceEnd;
    vector<string> seqsInFasta;
    string fastaFile;
    bool detected;
    bool hasFasta;
    // extend to make a cleaner trimming
    int trimmingExtension;
    // the threshold of edit_distance/match_length, suggest (0.1 ~ 0.3)
    double edMax;
};

class TrimmingOptions {
public:
    TrimmingOptions() {
        front = 0;
        tail = 0;
    }
public:
    // trimming first cycles for read
    int front;
    // trimming last cycles for read
    int tail;
    // max length of read
    int maxLen;
};

class QualityFilteringOptions {
public:
    QualityFilteringOptions() {
        enabled = true;
        // '0' = Q15
        qualifiedQual = '0';
        unqualifiedPercentLimit = 40;
        nBasePercentLimit = 10;
    }
public:
    // quality filter enabled
    bool enabled;
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    char qualifiedQual;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedPercentLimit;
    // if n_base_number > nBaseLimit, then discard this read
    int nBaseLimit;
    // if n_base_percent > nBasePercentLimit, then discard this read
    int nBasePercentLimit;
    // if average qual score < avgQualReq, then discard this read
    int avgQualReq;
};

class ReadLengthFilteringOptions {
public:
    ReadLengthFilteringOptions() {
        enabled = false;
        requiredLength = 15;
        maxLength = 0;
    }
public:
    // length filter enabled
    bool enabled;
    // if read_length < requiredLength, then this read is discard
    int requiredLength;
    // length limit, 0 for no limitation
    int maxLength;
};

class Options{
public:
    Options();
    void init();
    bool validate();
    bool adapterCuttingEnabled();
    bool polyXTrimmingEnabled();
    string getReadStartAdapter();
    string getReadEndAdapter();
    vector<string> makeListFromFileByLine(string filename);
    bool shallDetectAdapter();
    void loadFastaAdapters();

public:
    // file name of read input
    string in;
    // file name of read output
    string out;
    // file name of failed reads output
    string failedOut;
    // json file
    string jsonFile;
    // html file
    string htmlFile;
    // html report title
    string reportTitle;
    // compression level
    int compression;
    // do not rewrite existing files
    bool dontOverwrite;
    // read STDIN
    bool inputFromSTDIN;
    // write STDOUT
    bool outputToSTDOUT;
    // only process first N reads
    int readsToProcess;
    // worker thread number
    int thread;
    // trimming options
    int seqLen;
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;
    // length filtering options
    ReadLengthFilteringOptions lengthFilter;
    // adapter options
    AdapterOptions adapter;
    // multiple file splitting options
    SplitOptions split;
    // options for quality cutting
    QualityCutOptions qualityCut;
    // 3' end polyX trimming
    PolyXTrimmerOptions polyXTrim;
    // low complexity filtering
    LowComplexityFilterOptions complexityFilter;
    // N masking by quality
    MaskOptions mask;
    // break reads into high-quality fragments, and discard low-quality fragments
    BreakOptions breakOpt;
    // output debug information
    bool verbose;
    // the buffer size for writer
    size_t writerBufferSize;
    // is direct RNA sequencing data, which contains U not T
    bool isRNA;

};

#endif
