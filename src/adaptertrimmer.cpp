#include "adaptertrimmer.h"
#include "editdistance.h"

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}


bool AdapterTrimmer::trimByMultiSequences(Read* r, FilterResult* fr, vector<string>& adapterList, bool incTrimmedCounter) {
    int matchReq = 4;
    if(adapterList.size() > 16)
        matchReq = 5;
    if(adapterList.size() > 256)
        matchReq = 6;
    bool trimmed = false;

    string* originalSeq = r->mSeq;
    for(int i=0; i<adapterList.size(); i++) {
        trimmed |= trimBySequenceStart(r, NULL, adapterList[i], matchReq);
        trimmed |= trimBySequenceEnd(r, NULL, adapterList[i], matchReq);
    }

    if(trimmed) {
        string adapter = originalSeq->substr(r->length(), originalSeq->length() - r->length());
        if(fr)
            fr->addAdapterTrimmed(adapter, incTrimmedCounter);
        else
            cerr << adapter << endl;
    }

    return trimmed;
}

bool AdapterTrimmer::trimBySequenceStart(Read* r, FilterResult* fr, string& adapterseq, double edMax) {
    const int WINDOW = 200;
    const int PATTERN_LEN = 16;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(rlen < PATTERN_LEN)
        return false;

    int plen = min(PATTERN_LEN, alen);

    bool found = false;
    int pos = -1;
    for(int p=0; p<rlen-plen && p<WINDOW - plen; p++) {
        if(edit_distance(rdata + p, plen, adata + alen - plen, plen) < round(edMax * plen)) {
            // extend to compare the whole adapter
            int cmplen = min(p+plen, alen);
            if(edit_distance(rdata + p + plen - cmplen, cmplen, adata + alen - cmplen, cmplen) < round(edMax * cmplen) ){
                fr->addAdapterTrimmed(adapterseq.substr(alen - cmplen, cmplen));
                r->trimFront(p + plen);
                return true;
            }
        }
    }
    return false;
}

bool AdapterTrimmer::trimBySequenceEnd(Read* r, FilterResult* fr, string& adapterseq, double edMax) {
    const int WINDOW = 200;
    const int PATTERN_LEN = 16;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(rlen < PATTERN_LEN)
        return false;

    int plen = min(PATTERN_LEN, alen);

    bool found = false;
    int pos = -1;
    //from tail to front
    for(int p=0; p<rlen-plen && p<WINDOW - plen; p++) {
        if(edit_distance(rdata + rlen - plen -p, plen, adata, plen) < round(edMax * plen)) {
            // extend to compare the whole adapter
            int cmplen = min(p+plen, alen);
            if(edit_distance(rdata + rlen -plen - p, cmplen, adata, cmplen) < round(edMax * cmplen) ){
                fr->addAdapterTrimmed(adapterseq.substr(0, cmplen));
                r->resize(rlen - plen -p);
                return true;
            }
        }
    }
    return false;
}

bool AdapterTrimmer::test() {
    Read r("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////E");
    string adapter = "TTTTCCACGGGGATACTACTG";
    bool trimmed = AdapterTrimmer::trimBySequenceEnd(&r, NULL, adapter);
    if (*r.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA")
        return false;

    Read read("@name",
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
    }

    return true;
}