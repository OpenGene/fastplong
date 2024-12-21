#include "stats.h"
#include <memory.h>
#include <sstream>
#include "util.h"
#include "htmlreporter.h"

#define KMER_LEN 5

Stats::Stats(Options* opt, int guessedCycles, int bufferMargin){
    mOptions = opt;
    mReads = 0;
    mLengthSum = 0;

    mEvaluatedSeqLen = mOptions->seqLen;

    if(guessedCycles == 0) {
        guessedCycles = mEvaluatedSeqLen;
    }

    mCycles = guessedCycles;
    mBases = 0;
    mQ5Total = 0;
    mQ7Total = 0;
    mQ10Total = 0;
    mQ15Total = 0;
    mQ20Total = 0;
    mQ30Total = 0;
    mQ40Total = 0;
    summarized = false;
    mKmerMin = 0;
    mKmerMax = 0;


    mMinLen = 0;
    mMaxLen = 0;
    mMedianLen = 0;
    mN50Len = 0;
    mNeedCalcLength = true;

    // extend the buffer to make sure it's long enough
    mBufLen = guessedCycles + bufferMargin;

    for(int i=0; i<8; i++){
        mQ20Bases[i] = 0;
        mQ30Bases[i] = 0;
        mBaseContents[i] = 0;

        mCycleQ30Bases[i] = new long[mBufLen];
        memset(mCycleQ30Bases[i], 0, sizeof(long) * mBufLen);

        mCycleQ20Bases[i] = new long[mBufLen];
        memset(mCycleQ20Bases[i], 0, sizeof(long) * mBufLen);

        mCycleBaseContents[i] = new long[mBufLen];
        memset(mCycleBaseContents[i], 0, sizeof(long) * mBufLen);

        mCycleBaseQual[i] = new long[mBufLen];
        memset(mCycleBaseQual[i], 0, sizeof(long) * mBufLen);
    }
    mCycleTotalBase = new long[mBufLen];
    memset(mCycleTotalBase, 0, sizeof(long)*mBufLen);

    mCycleTotalQual = new long[mBufLen];
    memset(mCycleTotalQual, 0, sizeof(long)*mBufLen);

    mKmerBufLen = 2<<(KMER_LEN * 2);
    mKmer = new long[mKmerBufLen];
    memset(mKmer, 0, sizeof(long)*mKmerBufLen);

    memset(mBaseQualHistogram, 0, sizeof(long)*128);
    memset(mMedianReadQualHistogram, 0, sizeof(long)*128);
    memset(mMedianReadQualBases, 0, sizeof(long)*128);
}

void Stats::extendBuffer(int newBufLen){
    if(newBufLen <= mBufLen)
        return ;

    long* newBuf = NULL;

    for(int i=0; i<8; i++){
        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleQ30Bases[i], sizeof(long) * mBufLen);
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleQ20Bases[i], sizeof(long) * mBufLen);
        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleBaseContents[i], sizeof(long) * mBufLen);
        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = newBuf;

        newBuf = new long[newBufLen];
        memset(newBuf, 0, sizeof(long)*newBufLen);
        memcpy(newBuf, mCycleBaseQual[i], sizeof(long) * mBufLen);
        delete mCycleBaseQual[i];
        mCycleBaseQual[i] = newBuf;
    }
    newBuf = new long[newBufLen];
    memset(newBuf, 0, sizeof(long)*newBufLen);
    memcpy(newBuf, mCycleTotalBase, sizeof(long)*mBufLen);
    delete mCycleTotalBase;
    mCycleTotalBase = newBuf;

    newBuf = new long[newBufLen];
    memset(newBuf, 0, sizeof(long)*newBufLen);
    memcpy(newBuf, mCycleTotalQual, sizeof(long)*mBufLen);
    delete mCycleTotalQual;
    mCycleTotalQual = newBuf;

    mBufLen = newBufLen;
}

Stats::~Stats() {
    for(int i=0; i<8; i++){
        delete mCycleQ30Bases[i];
        mCycleQ30Bases[i] = NULL;

        delete mCycleQ20Bases[i];
        mCycleQ20Bases[i] = NULL;

        delete mCycleBaseContents[i];
        mCycleBaseContents[i] = NULL;

        delete mCycleBaseQual[i];
        mCycleBaseQual[i] = NULL;
    }

    delete mCycleTotalBase;
    delete mCycleTotalQual;

    // delete memory of curves
    map<string, double*>::iterator iter;
    for(iter = mQualityCurves.begin(); iter != mQualityCurves.end(); iter++) {
        delete iter->second;
    }
    for(iter = mContentCurves.begin(); iter != mContentCurves.end(); iter++) {
        delete iter->second;
    }
    delete mKmer;
}

void Stats::summarize(bool forced) {
    if(summarized && !forced)
        return;

    // first get the cycle and count total bases
    for(int c=0; c<mBufLen; c++) {
        mBases += mCycleTotalBase[c];
        if (mCycleTotalBase[c] == 0){
            mCycles = c;
            break;
        }
    }
    if(mCycleTotalBase[mBufLen-1]>0)
        mCycles = mBufLen;

    // Q20, Q30, base content
    for(int i=0; i<8; i++) {
        for(int c=0; c<mCycles; c++) {
            mQ20Bases[i] += mCycleQ20Bases[i][c];
            mQ30Bases[i] += mCycleQ30Bases[i][c];
            mBaseContents[i] += mCycleBaseContents[i][c];
        }
        mQ20Total += mQ20Bases[i];
        mQ30Total += mQ30Bases[i];
    }

    for(char c=40; c<127-33; c++) {
        mQ40Total += mBaseQualHistogram[c+33];
    }
    mQ30Total = mQ40Total;
    for(char c=30; c<40; c++) {
        mQ30Total += mBaseQualHistogram[c+33];
    }
    mQ20Total = mQ30Total;
    for(char c=20; c<30; c++) {
        mQ20Total += mBaseQualHistogram[c+33];
    }
    mQ15Total = mQ20Total;
    for(char c=15; c<20; c++) {
        mQ15Total += mBaseQualHistogram[c+33];
    }
    mQ10Total = mQ15Total;
    for(char c=10; c<15; c++) {
        mQ10Total += mBaseQualHistogram[c+33];
    }
    mQ7Total = mQ10Total;
    for(char c=7; c<10; c++) {
        mQ7Total += mBaseQualHistogram[c+33];
    }
    mQ5Total = mQ7Total;
    for(char c=5; c<7; c++) {
        mQ5Total += mBaseQualHistogram[c+33];
    }


    // quality curve for mean qual
    double* meanQualCurve = new double[mCycles];
    memset(meanQualCurve, 0, sizeof(double)*mCycles);
    for(int c=0; c<mCycles; c++) {
        meanQualCurve[c] = (double)mCycleTotalQual[c] / (double)mCycleTotalBase[c];
    }
    mQualityCurves["mean"] = meanQualCurve;

    // quality curves and base content curves for different nucleotides
    char alphabets[5] = {'A', 'T', 'C', 'G', 'N'};
    if(mOptions->isRNA)
        alphabets[1]='U';
    for(int i=0; i<5; i++) {
        char base = alphabets[i];
        // get last 3 bits
        char b = base & 0x07;
        double* qualCurve = new double[mCycles];
        memset(qualCurve, 0, sizeof(double)*mCycles);
        double* contentCurve = new double[mCycles];
        memset(contentCurve, 0, sizeof(double)*mCycles);
        for(int c=0; c<mCycles; c++) {
            if(mCycleBaseContents[b][c] == 0)
                qualCurve[c] = meanQualCurve[c];
            else
                qualCurve[c] = (double)mCycleBaseQual[b][c] / (double)mCycleBaseContents[b][c];
            contentCurve[c] = (double)mCycleBaseContents[b][c] / (double)mCycleTotalBase[c];
        }
        mQualityCurves[string(1, base)] = qualCurve;
        mContentCurves[string(1, base)] = contentCurve;
    }

    // GC content curve
    double* gcContentCurve = new double[mCycles];
    memset(gcContentCurve, 0, sizeof(double)*mCycles);
    char gBase = 'G' & 0x07;
    char cBase = 'C' & 0x07;
    for(int c=0; c<mCycles; c++) {
        gcContentCurve[c] = (double)(mCycleBaseContents[gBase][c] + mCycleBaseContents[cBase][c]) / (double)mCycleTotalBase[c];
    }
    mContentCurves["GC"] = gcContentCurve;

    mKmerMin = mKmer[0];
    mKmerMax = mKmer[0];
    for(int i=0; i<mKmerBufLen; i++) {
        if(mKmer[i] > mKmerMax)
            mKmerMax = mKmer[i];
        if(mKmer[i] < mKmerMin)
            mKmerMin = mKmer[i];
    }

    summarized = true;
}

int Stats::getMeanLength() {
    if(mReads == 0)
        return 0.0;
    else
        return mLengthSum/mReads;
}

void Stats::statRead(Read* r) {
    int len = r->length();

    mLengthVec.push_back(len);
    mNeedCalcLength = true;

    mLengthSum += len;

    if(mBufLen < len) {
        extendBuffer(max(len + 100, (int)(len * 1.5)));
    }
    const char* seqstr = r->mSeq->c_str();
    const char* qualstr = r->mQuality->c_str();

    int* qualHist = new int[128];
    memset(qualHist, 0, sizeof(int)*128);

    int kmer = 0;
    bool needFullCompute = true;
    for(int i=0; i<len; i++) {
        char base = seqstr[i];
        char qual = qualstr[i];
        // get last 3 bits
        char b = base & 0x07;

        const char q20 = '5';
        const char q30 = '?';

        mBaseQualHistogram[qual]++;
        qualHist[qual]++;

        if(qual >= q30) {
            mCycleQ30Bases[b][i]++;
            mCycleQ20Bases[b][i]++;
        } else if(qual >= q20) {
            mCycleQ20Bases[b][i]++;
        }

        mCycleBaseContents[b][i]++;
        mCycleBaseQual[b][i] += (qual-33);

        mCycleTotalBase[i]++;
        mCycleTotalQual[i] += (qual-33);

        if(base == 'N'){
            needFullCompute = true;
            continue;
        }

        // 5 bases required for kmer computing
        if(i<4)
            continue;

        // calc 5 KMER
        // 0x3FC == 0011 1111 1100
        if(!needFullCompute){
            int val = base2val(base);
            if(val < 0){
                needFullCompute = true;
                continue;
            } else {
                kmer = ((kmer<<2) & 0x3FC ) | val;
                mKmer[kmer]++;
            }
        } else {
            bool valid = true;
            kmer = 0;
            for(int k=0; k<5; k++) {
                int val = base2val(seqstr[i - 4 + k]);
                if(val < 0) {
                    valid = false;
                    break;
                }
                kmer = ((kmer<<2) & 0x3FC ) | val;
            }
            if(!valid) {
                needFullCompute = true;
                continue;
            } else {
                mKmer[kmer]++;
                needFullCompute = false;
            }
        }

    }

    //calculate the median
    if(len > 0) {
        int total = 0;
        char median = 0;
        int half = len>>1;
        while(true) {
            total += qualHist[median];
            if(total > half)
                break;
            median++;
        }
        mMedianReadQualHistogram[median]++;
        mMedianReadQualBases[median]+=len;

        map<char, vector<int>>::iterator iter = mQualLength.find(median);
        if(iter == mQualLength.end()) {
            mQualLength[median] = vector<int>();
        }
        mQualLength[median].push_back(len);
    }

    delete[] qualHist;

    mReads++;
}

void Stats::calcLengthHistogram() {
    for(size_t i=0; i< mLengthVec.size(); i++) {
        int len = mLengthVec[i];
        if(mLengthHist.count(len) > 0) {
            mLengthHist[len] ++;
        } else {
            mLengthHist[len] = 1;
        }
    }

    map<int, int>::iterator iter;
    long totalBase = 0;
    int readnum = 0;
    mN50Len = 0;
    mMedianLen = 0;
    if(!mLengthHist.empty()) {
        mMinLen = mLengthHist.begin()->first;
        mMaxLen = mLengthHist.rbegin()->first;
    }
    for(iter = mLengthHist.begin(); iter != mLengthHist.end(); iter++) {
        totalBase += iter->first * iter->second;
        if(mN50Len == 0 && totalBase > mLengthSum / 2)
            mN50Len = iter->first;
        readnum += iter->second;
        if(mMedianLen == 0 &&readnum > mLengthVec.size() / 2) 
            mMedianLen = iter->first;

        if(mMedianLen>0 && mN50Len>0)
            break;
    }

    mNeedCalcLength = false;
}

int Stats::base2val(char base) {
    switch(base){
        case 'A':
            return 0;
        case 'T':
        case 'U':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 3;
        default:
            return -1;
    }
}

int Stats::getCycles() {
    if(!summarized)
        summarize();
    return mCycles;
}

long Stats::getReads() {
    if(!summarized)
        summarize();
    return mReads;
}

long Stats::getBases() {
    if(!summarized)
        summarize();
    return mBases;
}

long Stats::getQ20() {
    if(!summarized)
        summarize();
    return mQ20Total;
}

long Stats::getQ30() {
    if(!summarized)
        summarize();
    return mQ30Total;
}

long Stats::getGCNumber() {
    if(!summarized)
        summarize();
    return mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07];
}

void Stats::print() {
    if(!summarized) {
        summarize();
    }
    cerr << "total reads: " << mReads << endl;
    cerr << "total bases: " << mBases << endl;
    cerr << "Q20 bases: " << mQ20Total << "(" << (mQ20Total*100.0)/mBases << "%)" << endl;
    cerr << "Q30 bases: " << mQ30Total << "(" << (mQ30Total*100.0)/mBases << "%)" << endl;
}

void Stats::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"total_reads\": " << mReads << "," << endl;
    ofs << padding << "\t" << "\"total_bases\": " << mBases << "," << endl;
    ofs << padding << "\t" << "\"q20_bases\": " << mQ20Total << "," << endl;
    ofs << padding << "\t" << "\"q30_bases\": " << mQ30Total << "," << endl;
    ofs << padding << "\t" << "\"total_cycles\": " << mCycles << "," << endl;

    // quality curves
    string qualNames[5] = {"A", "T", "C", "G", "mean"};
    if(mOptions->isRNA)
        qualNames[1]="U";
    ofs << padding << "\t" << "\"quality_curves\": {" << endl;
    for(int i=0 ;i<5; i++) {
        string name=qualNames[i];
        double* curve = mQualityCurves[name];
        ofs << padding << "\t\t" << "\"" << name << "\":[";
        for(int c = 0; c<mCycles; c++) {
            ofs << curve[c];
            // not the end
            if(c != mCycles - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(i != 5-1)
            ofs << ",";
        ofs << endl; 
    }
    ofs << padding << "\t" << "}," << endl;

    // content curves
    string contentNames[6] = {"A", "T", "C", "G", "N", "GC"};
    if(mOptions->isRNA)
        contentNames[1]="U";
    ofs << padding << "\t" << "\"content_curves\": {" << endl;
    for(int i=0 ;i<6; i++) {
        string name=contentNames[i];
        double* curve = mContentCurves[name];
        ofs << padding << "\t\t" << "\"" << name << "\":[";
        for(int c = 0; c<mCycles; c++) {
            ofs << curve[c];
            // not the end
            if(c != mCycles - 1)
                ofs << ",";
        }
        ofs << "]";
        // not the end;
        if(i != 6-1)
            ofs << ",";
        ofs << endl; 
    }
    ofs << padding << "\t" << "}," << endl;

    // KMER counting
    ofs << padding << "\t" << "\"kmer_count\": {" << endl;
    for(int i=0; i<64; i++) {
        string first = kmer3(i);
        for(int j=0; j<16; j++) {
            int target = (i<<4) + j;
            long count = mKmer[target];
            string last = kmer2(j);
            ofs << padding << "\t\t\"" << first << last << "\":" << count;
            if(j != 16-1)
                ofs << ",";
        }
        if(i != 64-1)
            ofs << "," << endl;
        else
            ofs << endl;
    }
    ofs << padding << "\t" << "}" << endl;

    ofs << padding << "}," << endl;
}

string Stats::list2string(double* list, long size, long* coords) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        // coords is 1,2,3,...
        long start = 0;
        if(i>0)
            start = coords[i-1];
        long end = coords[i];

        double total = 0.0;
        for(int k=start; k<end; k++)
            total += list[k];

        // get average
        if(end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

template<class T>
string Stats::list2string(T* list, long size) {
    stringstream ss;
    for(long i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

bool Stats::isLongRead() {
    return mCycles > 300;
}

void Stats::reporHtmlMedianQualHist(ofstream& ofs, string filteringType) {
    string subsection = filteringType + ": Read median quality statistics";
    string divName = replace(subsection, " ", "_");

    ofs << "<div class='subsection_title'>" + subsection + "</div>\n";

    int minVal = 0;
    int maxVal = 0;
    for(int i=0; i<127-33; i++) {
        if(mMedianReadQualBases[i+33] == 0)
            minVal++;
        else
            break;
    }
    for(int i=127-33; i>=0; i--) {
        if(mMedianReadQualBases[i+33] > 0) {
            maxVal = i;
            break;
        }
    }

    int offset = max(0, minVal-1); 
    int total = min(127-33, maxVal - minVal + 2);
    total = max(0, total);
    long *x = new long[total];
    for(int i=0; i<total; i++) {
        x[i] = i + offset;
    }

    double* percentReads = new double[total];
    memset(percentReads, 0, sizeof(double)*total);
    double* percentBases = new double[total];
    memset(percentBases, 0, sizeof(double)*total);
    for(int i=0; i<total; i++) {
        percentReads[i] = (double)mMedianReadQualHistogram[i+offset+33] * 100.0 / (double)mReads;
        percentBases[i] = (double)mMedianReadQualBases[i+offset+33] * 100.0 / (double)mBases;
    }

    ofs << "<div id='mean_qual_length_histogram_figure'>\n";
    ofs << "<div class='figure' id='plot_median_qual_hist_" + divName + "' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var readNum=";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percentReads, total) + "],";
    json_str += "name: '% reads',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += ";\n";


    json_str += "var baseNum=";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percentBases, total) + "],";
    json_str += "name: '% accumulated bases',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,128,0,1.0)', width:1}\n";
    json_str += "}";

    json_str += ";\n";

    json_str += "var data = [readNum, baseNum];;\n";

    json_str += "var layout={legend: {x: 0, y: 1.0},title:'Read median quality distribution', xaxis:{title:'read median quality score'}, yaxis:{title:'Percent (%)'}};\n";
    json_str += "Plotly.newPlot('plot_median_qual_hist_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percentReads;
    delete[] percentBases;
}

void Stats::reporHtmlMedianQualLengthDensity(ofstream& ofs, string filteringType) {
    string subsection = filteringType + ": Density plot of read median quality and read length";
    string divName = replace(subsection, " ", "_");

    ofs << "<div class='subsection_title'>" + subsection + "</div>\n";

    short *x = new short[mReads];
    int *y = new int[mReads];
    long cur = 0;
    map<char, vector<int>>::iterator iter;
    for(iter = mQualLength.begin(); iter != mQualLength.end(); iter++) {
        for(size_t l=0; l<iter->second.size(); l++) {
            x[cur] = iter->first - 33;
            y[cur] = iter->second[l];
            cur++;
        }
    }

    ofs << "<div id='mean_qual_length_density_figure'>\n";
    ofs << "<div class='figure' id='plot_median_qual_length_density_" + divName + "' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var density=";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, mReads) + "],";
    json_str += "y:[" + Stats::list2string(y, mReads) + "],";
    json_str += "name: '% reads',";
    json_str += "type:'histogram2dcontour',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += ";\n";


    json_str += "var data = [density];\n";

    json_str += "var layout={legend: {x: 0, y: 1.0},title:' Density plot of read median quality and read length', xaxis:{title:'read median quality score'}, yaxis:{title:'Read length', type:'log'}};\n";
    json_str += "Plotly.newPlot('plot_median_qual_length_density_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] y;
}


void Stats::reportHtmlBasicInfo(ofstream& ofs, string filteringType) {
    // KMER
    string subsection = filteringType + ": Basic statistics";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");

    ofs << "<div class='subsection_title'>" + subsection + "</div>\n";
    ofs << "<table>\n";
    HtmlReporter::outputRow(ofs, "total reads:", HtmlReporter::formatNumber(mReads));
    HtmlReporter::outputRow(ofs, "total bases:", HtmlReporter::formatNumber(mBases));
    HtmlReporter::outputRow(ofs, "minimum length:", HtmlReporter::formatNumber(mMinLen));
    HtmlReporter::outputRow(ofs, "maximum length:", HtmlReporter::formatNumber(mMaxLen));
    HtmlReporter::outputRow(ofs, "median length:", HtmlReporter::formatNumber(mMedianLen));
    HtmlReporter::outputRow(ofs, "mean length:", HtmlReporter::formatNumber(getMeanLength()));
    HtmlReporter::outputRow(ofs, "N50 length:", HtmlReporter::formatNumber(mN50Len));
    HtmlReporter::outputRow(ofs, "GC content:", HtmlReporter::getPercents(getGCNumber(),mBases) + "%");
    HtmlReporter::outputRow(ofs, "Q5 bases:", HtmlReporter::formatNumber(mQ5Total) + " (" + HtmlReporter::getPercents(mQ5Total,mBases) + "%)");
    HtmlReporter::outputRow(ofs, "Q7 bases:", HtmlReporter::formatNumber(mQ7Total) + " (" + HtmlReporter::getPercents(mQ7Total, mBases) + "%)");
    HtmlReporter::outputRow(ofs, "Q10 bases:", HtmlReporter::formatNumber(mQ10Total) + " (" + HtmlReporter::getPercents(mQ10Total,mBases) + "%)");
    HtmlReporter::outputRow(ofs, "Q15 bases:", HtmlReporter::formatNumber(mQ15Total) + " (" + HtmlReporter::getPercents(mQ15Total, mBases) + "%)");
    HtmlReporter::outputRow(ofs, "Q20 bases:", HtmlReporter::formatNumber(mQ20Total) + " (" + HtmlReporter::getPercents(mQ20Total,mBases) + "%)");
    HtmlReporter::outputRow(ofs, "Q30 bases:", HtmlReporter::formatNumber(mQ30Total) + " (" + HtmlReporter::getPercents(mQ30Total, mBases) + "%)");
    HtmlReporter::outputRow(ofs, "Q40 bases:", HtmlReporter::formatNumber(mQ40Total) + " (" + HtmlReporter::getPercents(mQ40Total, mBases) + "%)");
    ofs << "</table>\n";
}

void Stats::reportHtmlKMER(ofstream& ofs, string filteringType) {

    // KMER
    string subsection = filteringType + ": KMER counting";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'>" + subsection + "</div>\n";
    ofs << "<div  id='" << divName << "'>\n";
    ofs << "<div class='sub_section_tips'>Darker background means larger counts. The count will be shown on mouse over.</div>\n";
    ofs << "<table class='kmer_table' style='width:680px;'>\n";
    ofs << "<tr>";
    ofs << "<td></td>";
    // the heading row
    for(int h=0; h<16; h++) 
        ofs << "<td style='color:#333333'>" << kmer2(h) << "</td>";
    ofs << "</tr>\n";
    // content
    for(int i=0; i<64; i++) {
        ofs << "<tr>";

        ofs << "<td style='color:#333333'>" << kmer3(i) << "</td>";
        for(int j=0; j<16; j++) {
            ofs << makeKmerTD(i,j) ;
        }
        ofs << "</tr>\n";
    }
    ofs << "</table>\n";
    ofs << "</div>\n";
}

string Stats::makeKmerTD(int i, int j) {
    int target = (i<<4) + j;
    long val = mKmer[target];
    // 3bp + 2bp = 5bp
    string first = kmer3(i);
    string last = kmer2(j);
    string kmer = first+last;
    double meanBases = (double)(mBases+1) / mKmerBufLen;
    double prop = val / meanBases;
    /*double frac = 0.5;
    if(prop > 2.0) 
        frac = (prop-2.0)/20.0 + 0.5;
    else if(prop< 0.5)
        frac = prop;

    frac = max(0.01, min(1.0, frac));*/
    int r =0;
    int g = 0;
    int b = 0;
    if(prop <= 0.3){
        double frac = prop*2.0;
        b = 255 - 256*frac;
        g = 255*frac;
        r = b*frac;
    } else if(prop>3.0){
        double frac = 2.0/prop;
        r = 255  - 128*frac;
        g = 128*frac;
        b = r*frac;
    } else {
        r = g = b = 196;
    }
    stringstream ss;
    ss << "<td style='background:#"; 
    if(r<16)
        ss << "0";
    ss<<hex<<r;
    if(g<16)
        ss << "0";
    ss<<hex<<g;
    if(b<16)
        ss << "0";
    ss<<hex<<b;
    ss << dec << "' title='"<< kmer << ": " << val << "\n" << prop << " times as mean value'>";
    ss << kmer << "</td>";
    return ss.str();
}

string Stats::kmer3(int val) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    if(mOptions->isRNA)
        bases[1]='U';
    string ret(3, ' ');
    ret[0] = bases[(val & 0x30) >> 4];
    ret[1] = bases[(val & 0x0C) >> 2];
    ret[2] = bases[(val & 0x03)];
    return ret;
}

string Stats::kmer2(int val) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    if(mOptions->isRNA)
        bases[1]='U';
    string ret(2, ' ');
    ret[0] = bases[(val & 0x0C) >> 2];
    ret[1] = bases[(val & 0x03)];
    return ret;
}

void Stats::reportHtmlQuality(ofstream& ofs, string filteringType) {

    // quality
    string subsection = filteringType  + ": quality";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'>" + subsection + "</div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    string alphabets[5] = {"A", "T", "C", "G", "mean"};
    if(mOptions->isRNA)
        alphabets[1]="U";
    string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<5; b++) {
        string base = alphabets[b];
        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mQualityCurves[base], total, x) + "],";
        json_str += "name: '" + base + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'quality'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

void Stats::reportHtmlContents(ofstream& ofs, string filteringType) {

    // content
    string subsection = filteringType + ": base contents";
    string divName = replace(subsection, " ", "_");
    divName = replace(divName, ":", "_");
    string title = "";

    ofs << "<div class='subsection_title'>" + subsection + "</div>\n";
    ofs << "<div id='" + divName + "'>\n";
    ofs << "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n";
    ofs << "<div class='figure' id='plot_" + divName + "'></div>\n";
    ofs << "</div>\n";
    
    string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
    if(mOptions->isRNA)
        alphabets[1]="U";
    string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(255, 0, 0, 1.0)", "rgba(20,20,20,1.0)"};
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    long *x = new long[mCycles];
    int total = 0;
    if(!isLongRead()) {
        for(int i=0; i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
    } else {
        const int fullSampling = 40;
        for(int i=0; i<fullSampling && i<mCycles; i++){
            x[total] = i+1;
            total++;
        }
        // down sampling if it's too long
        if(mCycles>fullSampling) {
            double pos = fullSampling;
            while(true){
                pos *= 1.05;
                if(pos >= mCycles)
                    break;
                x[total] = (int)pos;
                total++;
            }
            // make sure lsat one is contained
            if(x[total-1] != mCycles){
                x[total] = mCycles;
                total++;
            }
        }
    }
    // four bases
    for (int b = 0; b<6; b++) {
        string base = alphabets[b];
        long count = 0;
        if(base.size()==1) {
            char b = base[0] & 0x07;
            count = mBaseContents[b];
        } else {
            count = mBaseContents['G' & 0x07] + mBaseContents['C' & 0x07] ;
        }
        string percentage = to_string((double)count * 100.0 / mBases);
        if(percentage.length()>5)
            percentage = percentage.substr(0,5);
        string name = base + "(" + percentage + "%)"; 

        json_str += "{";
        json_str += "x:[" + list2string(x, total) + "],";
        json_str += "y:[" + list2string(mContentCurves[base], total, x) + "],";
        json_str += "name: '" + name + "',";
        json_str += "mode:'lines',";
        json_str += "line:{color:'" + colors[b] + "', width:1}\n";
        json_str += "},";
    }
    json_str += "];\n";
    json_str += "var layout={title:'" + title + "', xaxis:{title:'position'";
    // use log plot if it's too long
    if(isLongRead()) {
        json_str += ",type:'log'";
    }
    json_str += "}, yaxis:{title:'base content ratios'}};\n";
    json_str += "Plotly.newPlot('plot_" + divName + "', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
}

Stats* Stats::merge(vector<Stats*>& list) {
    if(list.size() == 0)
        return NULL;

    //get the most long cycles
    int cycles = 0;
    for(int t=0; t<list.size(); t++) {
        list[t]->summarize();
        cycles = max(cycles, list[t]->getCycles());
    }

    Stats* s = new Stats(list[0]->mOptions, cycles, 0);

    // init overrepresented seq maps
    map<string, long>::iterator iter;

    for(int t=0; t<list.size(); t++) {
        int curCycles =  list[t]->getCycles();
        // merge read number
        s->mReads += list[t]->mReads;
        s->mLengthSum += list[t]->mLengthSum;

        // merge length vector
        s->mLengthVec.reserve(s->mLengthVec.size() + list[t]->mLengthVec.size());
        s->mLengthVec.insert(s->mLengthVec.end(), list[t]->mLengthVec.begin(), list[t]->mLengthVec.end());

        // merge per cycle counting for different bases
        for(int i=0; i<8; i++){
            for(int j=0; j<cycles && j<curCycles; j++) {
                s->mCycleQ30Bases[i][j] += list[t]->mCycleQ30Bases[i][j];
                s->mCycleQ20Bases[i][j] += list[t]->mCycleQ20Bases[i][j];
                s->mCycleBaseContents[i][j] += list[t]->mCycleBaseContents[i][j];
                s->mCycleBaseQual[i][j] += list[t]->mCycleBaseQual[i][j];
            }
        }

        // merge per cycle counting for all bases
        for(int j=0; j<cycles && j<curCycles; j++) {
            s->mCycleTotalBase[j] += list[t]->mCycleTotalBase[j];
            s->mCycleTotalQual[j] += list[t]->mCycleTotalQual[j];
        }

        // merge kMer
        for(int i=0; i<s->mKmerBufLen; i++) {
            s->mKmer[i] += list[t]->mKmer[i];
        }

        // merge base/read qual histogram
        for(int i=0; i<128; i++) {
            s->mBaseQualHistogram[i] += list[t]->mBaseQualHistogram[i];
            s->mMedianReadQualHistogram[i] += list[t]->mMedianReadQualHistogram[i];
            s->mMedianReadQualBases[i] += list[t]->mMedianReadQualBases[i];
        }

        // merge qual-length distribution
        map<char, vector<int>>::iterator iter;
        map<char, vector<int>>::iterator overallIter;
        for(iter = list[t]->mQualLength.begin(); iter != list[t]->mQualLength.end(); iter++) {
            overallIter = s->mQualLength.find(iter->first);
            if(overallIter == s->mQualLength.end()) {
                s->mQualLength[iter->first] = vector<int>();
            }
            s->mQualLength[iter->first].insert(s->mQualLength[iter->first].end(), iter->second.begin(), iter->second.end());
        }
    }

    s->summarize();

    return s;
}
