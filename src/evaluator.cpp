#include "evaluator.h"
#include "fastqreader.h"
#include <map>
#include <memory.h>
#include "nucleotidetree.h"
#include "knownadapters.h"

Evaluator::Evaluator(Options* opt){
    mOptions = opt;
}


Evaluator::~Evaluator(){
}

void Evaluator::evaluateSeqLenAndCheckRNA() {
    if(mOptions->in.empty()) {
        return;
    }

    FastqReader reader(mOptions->in);

    long records = 0;
    bool reachedEOF = false;
    bool isRNA = false;

    // get seqlen
    int seqlen=0;
    long numT = 0;
    long numU = 0;
    while(records < 100) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        const char* seq = r->mSeq->data();
        int rlen = r->length();
        if(rlen > seqlen)
            seqlen = rlen;
        for (int c=0; c<rlen; c++) {
            const char base = seq[c];
            if(base == 'T')
                numT++;
            else if(base == 'U')
                numU++;
        }
        records ++;
        delete r;
    }
    if(numT > 0 && numU >0) {
        error_exit("This data contains both U and T");
    } else if (numU > 0) {
        isRNA = true;
        cerr << "RNA direct sequencing data" << endl;
    }

    mOptions->seqLen = seqlen;
    mOptions->isRNA = isRNA;

}

void Evaluator::evaluateReadNum(long& readNum) {
    FastqReader reader(mOptions->in);

    const long READ_LIMIT = 512*1024;
    const long BASE_LIMIT = 151 * 512*1024;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    bool reachedEOF = false;
    bool first = true;
    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        records++;
        bases += r->length();
        delete r;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }
}

void Evaluator::evalAdapterAndReadNum(Options* opt, long& readNum) {
    string filename = mOptions->in;

    FastqReader reader(filename);
    // stat up to 64K reads or 512M bases
    const long READ_LIMIT = 64*1024;
    const long BASE_LIMIT = 8192 * READ_LIMIT;
    long records = 0;
    long bases = 0;
    size_t firstReadPos = 0;

    size_t bytesRead;
    size_t bytesTotal;

    Read** loadedReads = new Read*[READ_LIMIT];
    memset(loadedReads, 0, sizeof(Read*)*READ_LIMIT);
    bool reachedEOF = false;
    bool first = true;

    while(records < READ_LIMIT && bases < BASE_LIMIT) {
        Read* r = reader.read();
        if(!r) {
            reachedEOF = true;
            break;
        }
        if(first) {
            reader.getBytes(bytesRead, bytesTotal);
            firstReadPos = bytesRead;
            first = false;
        }
        int rlen = r->length();
        bases += rlen;
        loadedReads[records] = r;
        records++;
    }

    readNum = 0;
    if(reachedEOF){
        readNum = records;
    } else if(records>0) {
        // by the way, update readNum so we don't need to evaluate it if splitting output is enabled
        reader.getBytes(bytesRead, bytesTotal);
        double bytesPerRead = (double)(bytesRead - firstReadPos) / (double) records;
        // increase it by 1% since the evaluation is usually a bit lower due to bad quality causes lower compression rate
        readNum = (long) (bytesTotal*1.01 / bytesPerRead);
    }

    // we need at least 100 valid records to evaluate
    if(records < 100) {
        for(int r=0; r<records; r++) {
            delete loadedReads[r];
            loadedReads[r] = nullptr;
        }
        delete[] loadedReads;
        return ;
    }

    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail);

    // why we add trim_tail here? since the last cycle are usually with low quality and should be trimmed
    const double FOLD_THRESHOLD = 100.0;
    const int keylen = 10;
    int size = 1 << (keylen*2 );
    unsigned int* counts = new unsigned int[size];
    unsigned long* positionAcc = new unsigned long[size];

    // read end adapter
    if(opt->adapter.sequenceStart == "auto") {
        cerr << "Trying to detect adapter sequence at read start"<<endl;
        long total = 0;
        int totalKey = 0;
        memset(counts, 0, sizeof(unsigned int)*size);
        memset(positionAcc, 0, sizeof(unsigned long)*size);
        for(int i=0; i<records; i++) {
            Read* r = loadedReads[i];
            const char* data = r->mSeq->c_str();
            int key = -1;
            for(int pos = 0; pos <= r->length()-keylen-shiftTail && pos<128; pos++) {
                key = seq2int(r->mSeq, pos, keylen, key);
                if(key >= 0) {
                    counts[key]++;
                    positionAcc[key]+=pos;
                    total++;
                }
            }
        }
        for(int k=0; k<size; k++) {
        if(counts[k] >0)
            totalKey++;
        }
        // set AAAAAAAAAA = 0;
        counts[0] = 0;

        int key = getTopKey(counts, keylen);
        long count = counts[key];
        if(count>10 && count*totalKey > total * FOLD_THRESHOLD) {
            string adapter = extendKeyToAdapter(key, counts, positionAcc, keylen, false);
            if(adapter.length() > 16){
                cerr << "Detected: " << adapter << endl;
                mOptions->adapter.sequenceStart = adapter;
            } else {
                cerr << "Found possible adapter sequence, but it's too short: " << adapter << ", specify -s " << adapter << " to force trimming using this adapter"  << endl;
            }
        } else {
            cerr << "Not detected" << endl;
        }
    }

    // read start adapter
    if(opt->adapter.sequenceEnd == "auto") {
        cerr << "Trying to detect adapter sequence at read end"<<endl;
        long total = 0;
        int totalKey = 0;
        memset(counts, 0, sizeof(unsigned int)*size);
        memset(positionAcc, 0, sizeof(unsigned long)*size);
        for(int i=0; i<records; i++) {
            Read* r = loadedReads[i];
            const char* data = r->mSeq->c_str();
            int key = -1;
            int startpos = max(0, r->length()-keylen-shiftTail - 128);
            for(int pos = startpos; pos <= r->length()-keylen-shiftTail; pos++) {
                key = seq2int(r->mSeq, pos, keylen, key);
                if(key >= 0) {
                    counts[key]++;
                    positionAcc[key]+=r->length() - pos;
                    total++;
                }
            }
        }
        for(int k=0; k<size; k++) {
        if(counts[k] >0)
            totalKey++;
        }
        // set AAAAAAAAAA = 0;
        counts[0] = 0;

        int key = getTopKey(counts, keylen);
        long count = counts[key];
        if(count>10 && count*totalKey > total * FOLD_THRESHOLD) {
            string adapter = extendKeyToAdapter(key, counts, positionAcc, keylen, mOptions->isRNA, true);
            if(adapter.length() > 16){
                cerr << "Detected: " << adapter << endl;
                mOptions->adapter.sequenceEnd = adapter;
            } else {
                cerr << "Found possible adapter sequence, but it's too short: " << adapter << ", specify -e " << adapter << " to force trimming using this adapter"  << endl;
            }
        } else {
            cerr << "Not detected" << endl;
        }
    }

    delete[] counts;
    delete[] positionAcc;
    for(int r=0; r<records; r++) {
        delete loadedReads[r];
        loadedReads[r] = nullptr;
    }
    delete[] loadedReads;

}

int Evaluator::getTopKey(unsigned int* counts, int keylen) {
    // get the top N
    int size = 1 << (keylen*2 );
    int topkey = -1;
    unsigned int topCount = 0;
    for(int k=0; k<size; k++) {
        unsigned int val = counts[k];

        // skip low complexity seq
        int atcg[4] = {0};
        for(int i=0; i<keylen; i++) {
            int baseOfBit = (k >> (i*2)) & 0x03;
            atcg[baseOfBit]++;
        }
        bool lowComplexity = false;
        int zeroNum = 0;
        for(int b=0; b<4; b++) {
            if(atcg[b] >= keylen-4)
                lowComplexity=true;
            if(atcg[b] == 0)
                zeroNum++;
        }
        if(zeroNum>=2)
            lowComplexity=true;
        //repeative
        if((k >> keylen) == (k & ((0x01<<keylen)-1)) )
            lowComplexity = true;
        int diff = 0;
        for(int s=0; s<keylen - 1; s++) {
            int curBase = ( val >> ((keylen -s)*2) ) & 0x03;
            int lastBase = ( val >> ((keylen -s - 1)*2) ) & 0x03;
            if(curBase != lastBase)
                diff++;
        }
        if(diff <3){
            continue;
        }
        if(lowComplexity)
            continue;
        // too many GC
        if(atcg[2] + atcg[3] >= keylen-2)
            continue;

        // starts with GGGG
        if( (k>>12) == 0xff)
            continue;
        if(k == 0)
            continue;
        
        if(val > topCount) {
            topCount = val;
            topkey = k;
        }
    }
    return topkey;
}

string Evaluator::extendKeyToAdapter(int key, unsigned int* counts, unsigned long* positionAcc, int keylen, bool isRNA, bool leftFirst) {
    string adapter = int2seq(key, keylen, isRNA);
    const int mask = (1 << (keylen*2 )) - 1;
    const int MAX_LEN = 64;
    char bases[4] = {'A', 'T', 'C', 'G'};
    if(isRNA)
        bases[1] = 'U';

    bool leftFinished = false;
    bool rightFinished = false;
    bool extendingLeft = true;
    if(!leftFirst)
        extendingLeft = false;
    while(true) {
        int curkey = key;
        while(adapter.length() < MAX_LEN) {
            int totalCount = 0;
            bool extended = false;
            for(int b=0; b<4; b++) {
                int newkey;
                if(extendingLeft)
                    newkey = (b<<((keylen-1)*2)) | (curkey >> 2);
                else
                    newkey = b | (mask & (curkey << 2));
                totalCount += counts[newkey];
            }

            for(int b=0; b<4; b++) {
                int newkey;
                if(extendingLeft)
                    newkey = (b<<((keylen-1)*2)) | (curkey >> 2);
                else
                    newkey = b | (mask & (curkey << 2));
                if(counts[newkey] == 0)
                    continue;
                double offset = (double)positionAcc[newkey]/counts[newkey] - (double)positionAcc[curkey]/counts[curkey];
                //cerr << (double)counts[newkey] / (double)totalCount << ", " << (double)counts[newkey] / (double)counts[key];
                //cerr << ", offset: " << offset << endl;
                if((double)counts[newkey] / (double)totalCount < 0.7)
                    continue;
                if( (double)counts[newkey] / (double)counts[key] < 0.5) {
                    continue;
                }
                // offset should be near -1.0
                if(offset > 2  || offset < -4)
                    continue;
                //successful to extend
                curkey = newkey;
                extended = true;
                if(extendingLeft) {
                    adapter.insert(adapter.begin(), bases[b]);
                } else {
                    adapter.insert(adapter.end(), bases[b]);
                }
                break;
            }

            if(!extended) {
                if(extendingLeft)
                    leftFinished = true;
                else
                    rightFinished = true;
                break;
            }

            if(adapter.length() == MAX_LEN) {
                leftFinished = true;
                rightFinished = true;
                break;
            }
        }

        //finish one side, go the other side
        if(extendingLeft)
            extendingLeft = false;
        else
            extendingLeft = true;

        if(leftFinished && rightFinished)
            break;
    }

    return adapter;

}

string Evaluator::getAdapterWithSeed(int seed, Read** loadedReads, long records, int keylen) {
    // we have to shift last cycle for evaluation since it is so noisy, especially for Illumina data
    const int shiftTail = max(1, mOptions->trim.tail);
    NucleotideTree forwardTree(mOptions);
    // forward search
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq->c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail; pos++) {
            key = seq2int(r->mSeq, pos, keylen, key);
            if(key == seed) {
                forwardTree.addSeq(r->mSeq->substr(pos+keylen, r->length()-keylen-shiftTail-pos));
            }
        }
    }
    bool reachedLeaf = true;
    string forwardPath = forwardTree.getDominantPath(reachedLeaf);

    NucleotideTree backwardTree(mOptions);
    // backward search
    for(int i=0; i<records; i++) {
        Read* r = loadedReads[i];
        const char* data = r->mSeq->c_str();
        int key = -1;
        for(int pos = 20; pos <= r->length()-keylen-shiftTail; pos++) {
            key = seq2int(r->mSeq, pos, keylen, key);
            if(key == seed) {
                string seq =  r->mSeq->substr(0, pos);
                string rcseq = reverse(seq);
                backwardTree.addSeq(rcseq);
            }
        }
    }
    string backwardPath = backwardTree.getDominantPath(reachedLeaf);

    string adapter = reverse(backwardPath) + int2seq(seed, keylen) + forwardPath;
    if(adapter.length()>60)
        adapter.resize(60);

    string matchedAdapter = matchKnownAdapter(adapter);
    if(!matchedAdapter.empty()) {
        map<string, string> knownAdapters = getKnownAdapter();
        cerr << knownAdapters[matchedAdapter] << endl << matchedAdapter << endl;
        return matchedAdapter;
    } else {
        if(reachedLeaf) {
            cerr << adapter << endl;
            return adapter;
        } else {
            return "";
        }
    }
}

string Evaluator::matchKnownAdapter(string seq) {
    map<string, string> knownAdapters = getKnownAdapter();
    map<string, string>::iterator iter;
    for(iter = knownAdapters.begin(); iter != knownAdapters.end(); iter++) {
        string adapter = iter->first;
        string desc = iter->second;
        if(seq.length()<adapter.length()) {
            continue;
        }
        int diff = 0;
        for(int i=0; i<adapter.length() && i<seq.length(); i++) {
            if(adapter[i] != seq[i])
                diff++;
        }
        if(diff == 0)
            return adapter;
    }
    return "";
}

string Evaluator::int2seq(unsigned int val, int seqlen, bool isRNA) {
    char bases[4] = {'A', 'T', 'C', 'G'};
    if(isRNA)
        bases[1] = 'U';
    string ret(seqlen, 'N');
    int done = 0;
    while(done < seqlen) {
        ret[seqlen - done - 1] = bases[val & 0x03];
        val = (val >> 2);
        done++;
    }
    return ret;
}

int Evaluator::seq2int(string* seq, int pos, int keylen, int lastVal) {
    return seq2int(*seq, pos, keylen, lastVal);
}

int Evaluator::seq2int(string& seq, int pos, int keylen, int lastVal) {
    int rlen = seq.length();
    if(lastVal >= 0) {
        const int mask = (1 << (keylen*2 )) - 1;
        int key = (lastVal<<2) & mask;
        char base = seq[pos + keylen - 1];
        switch (base) {
            case 'A':
                key += 0;
                break;
            case 'T':
            case 'U':
                key += 1;
                break;
            case 'C':
                key += 2;
                break;
            case 'G':
                key += 3;
                break;
            default:
                // N or anything else
                return -1;
        }
        return key;
    } else {
        int key = 0;
        for(int i=pos; i<keylen+pos; i++) {
            key = (key << 2);
            char base = seq[i];
            switch (base) {
                case 'A':
                    key += 0;
                    break;
                case 'T':
                case 'U':
                    key += 1;
                    break;
                case 'C':
                    key += 2;
                    break;
                case 'G':
                    key += 3;
                    break;
                default:
                    // N or anything else
                    return -1;
            }
        }
        return key;
    }
}

