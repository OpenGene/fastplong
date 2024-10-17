#include "adaptertrimmer.h"
#include "editdistance.h"
#include <math.h>

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::findMiddleAdapters(Read* r, string& startAdater, string& endAdapter, int& start, int& len, double edMax) {
    len = -1;

    int startAdaterPos = searchAdapter(r->mSeq, startAdater, edMax);
    int endAdapterPos = searchAdapter(r->mSeq, endAdapter, edMax);

     // extend it to make a cleaner cut
    const int EXTEND = 20;

    if(startAdaterPos >=0 && endAdapterPos>=0) {
        start = min(startAdaterPos, endAdapterPos);
        int end = max(startAdaterPos + startAdater.length(), endAdapterPos + endAdapter.length());

        start = max(0, start-EXTEND);
        end = min(r->length(), end+EXTEND);
        len = end - start;
        return true;
    } if(startAdaterPos >=0){
        int end = min(r->length(), startAdaterPos + (int)startAdater.length() +EXTEND); 
        start = max(0, startAdaterPos-EXTEND);
        len = end - start;
        return true;
    } if(endAdapterPos >=0){
        int end = min(r->length(), endAdapterPos + (int)endAdapter.length() +EXTEND); 
        start = max(0, endAdapterPos-EXTEND);
        len = end - start;
        return true;
    }

    return false;
}

int AdapterTrimmer::trimByMultiSequences(Read* r, FilterResult* fr, vector<string>& adapterList, double edMax) {
    int matchReq = 4;
    if(adapterList.size() > 16)
        matchReq = 5;
    if(adapterList.size() > 256)
        matchReq = 6;
    int trimmed = 0;

    string* originalSeq = r->mSeq;
    for(int i=0; i<adapterList.size(); i++) {
        trimmed += trimBySequenceStart(r, fr, adapterList[i], edMax);
        trimmed += trimBySequenceEnd(r, fr, adapterList[i], edMax);
    }

    return trimmed;
}

int AdapterTrimmer::searchAdapter(string* read, string& adapter, double edMax, int searchStart, int searchLen) {
    int minMismatch = 99999; // initialized with a large mismatch
    int pos = -1;
    // for the best match
    const char* adata = adapter.c_str();
    const char* rdata = read->c_str();
    int rlen = read->length();
    int alen = adapter.length();

    int threshold = round(edMax * alen);

    int searchEnd = rlen;
    if(searchLen > 0) {
        searchEnd = min(rlen, searchLen + searchStart);
    }

    for(int p = searchStart; p < searchEnd - alen; p++) {
        int mismatch = 0;
        for(int i=0; i<alen; i++) {
            if(rdata[p+i] != adata[i])
                mismatch++;
        }
        if(mismatch < minMismatch ) {
            minMismatch = mismatch;
            pos = p;
        }
    }

    if(pos >= 0 ) {
        int ed = edit_distance(rdata+pos, alen, adata, alen);
        if(ed < threshold)
            return pos;
        else
            return -1;
    } else
        return -1;
    
}

int AdapterTrimmer::trimBySequenceStart(Read* r, FilterResult* fr, string& adapterseq, double edMax) {
    const int WINDOW = 200;
    const int PATTERN_LEN = 16;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(rlen < PATTERN_LEN)
        return 0;

    int plen = min(PATTERN_LEN, alen);
    int pos = -1;

    // search by full match
    pos = searchAdapter(r->mSeq, adapterseq, edMax, 0, WINDOW);
    if(pos >= 0) {
        int cmplen = min(pos+plen, alen);
        if(fr) {
            fr->addAdapterTrimmed(adapterseq.substr(alen - cmplen, cmplen));
            r->trimFront(pos + plen);
            return pos+plen;
        }
    }

    // adapter not found by above full match
    // search part adapter

    int mined = -1;
    //from tail to front, search by partly match of edit distance
    for(int p=0; p<rlen-plen && p<WINDOW - plen; p++) {
        int ed = edit_distance(rdata + p, plen, adata + alen - plen, plen);
        if(ed< round(edMax * plen)) {
            if(pos < 0) {
                pos = p;
                mined = ed;
            } else if(ed > mined) // last one is best
                break;
            else {
                pos = p;
                mined = ed;
            }
        }
    }

    if(pos > 0) {
        // extend to compare the whole adapter
        int cmplen = min(pos+plen, alen);
        if(edit_distance(rdata + pos + plen - cmplen, cmplen, adata + alen - cmplen, cmplen) < round(edMax * cmplen) ){
            if(fr)
                fr->addAdapterTrimmed(adapterseq.substr(alen - cmplen, cmplen));
            r->trimFront(pos + plen);
            return pos+plen;
        }
    }

    return 0;
}

int AdapterTrimmer::trimBySequenceEnd(Read* r, FilterResult* fr, string& adapterseq, double edMax) {
    const int WINDOW = 200;
    const int PATTERN_LEN = 16;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(rlen < PATTERN_LEN)
        return false;

    int plen = min(PATTERN_LEN, alen);
    int pos = -1;

    // search by full match
    int searchStart = max(0, rlen - WINDOW);
    pos = searchAdapter(r->mSeq, adapterseq, edMax, searchStart, WINDOW);
    if(pos >= 0) {
        int cmplen = min(pos+plen, alen);
        if(fr)
            fr->addAdapterTrimmed(adapterseq.substr(0, cmplen));
        r->resize(rlen - plen -pos);
        return pos + plen;
    }

    // adapter not found by above full match
    // search part adapter

    int mined = -1;
    //from tail to front
    for(int p=0; p<rlen-plen && p<WINDOW - plen; p++) {
        int ed = edit_distance(rdata + rlen - plen -p, plen, adata, plen);
        if(ed< round(edMax * plen)) {
            if(pos < 0) {
                pos = p;
                mined = ed;
            } else if(ed > mined) // last one is best
                break;
            else {
                pos = p;
                mined = ed;
            }
        }
    }

    if(pos > 0) {
        // extend to compare the whole adapter
        int cmplen = min(pos+plen, alen);
        if(edit_distance(rdata + rlen -plen - pos, cmplen, adata, cmplen) < round(edMax * cmplen) ){
            if(fr)
                fr->addAdapterTrimmed(adapterseq.substr(0, cmplen));
            r->resize(rlen - plen -pos);
            return pos + plen;
        }
    }
    return 0;
}

bool AdapterTrimmer::test() {
    Read r("@name",
        "AGGTGCTGCGCATACTTTTCCACGGGGATACTACTGGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEE///EEEEEEEEEEEEEEEEEEEEE");
    string adapter = "GCGCATACTTTTCCACGGGGATACTACTG";
    int trimmed = AdapterTrimmer::trimBySequenceStart(&r, NULL, adapter);
    if (*r.mSeq != "GGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT") {
        cerr << "expect GGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT" << endl;
        cerr << "return " << *r.mSeq  << endl;
        return false;
    }

    Read r2("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAAGCGCATACTTTTCCACGGGGA",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEET");
    trimmed = AdapterTrimmer::trimBySequenceEnd(&r2, NULL, adapter);
    if (*r2.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA") {
        cerr << "expect TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA" << endl;
        cerr << "return " << *r.mSeq  << endl;
        return false;
    }

    /*Read read("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    vector<string> adapterList;
    adapterList.push_back("GCTAGCTAGCTAGCTA");
    adapterList.push_back("AAATTTCCCGGGAAATTTCCCGGG");
    adapterList.push_back("ATCGATCGATCGATCG");
    adapterList.push_back("AATTCCGGAATTCCGG");
    trimmed = AdapterTrimmer::trimByMultiSequences(&read, NULL, adapterList);
    if (*read.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG") {
        cerr << read.mSeq << endl;
        return false;
    }*/

    return true;
}