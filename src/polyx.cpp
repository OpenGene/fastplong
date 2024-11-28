#include "polyx.h"
#include "common.h"

PolyX::PolyX(){
}


PolyX::~PolyX(){
}

void PolyX::trimPolyX(Read* r, FilterResult* fr, int compareReq) {
    const int allowOneMismatchForEach = 8;
    const int maxMismatch = 5;

    const char* data = r->mSeq->c_str();

    int rlen = r->length();


    int atcgNumbers[4] = {0, 0, 0, 0};
    int pos = 0;
    for(pos=0; pos<rlen; pos++) {
        switch(data[rlen - pos - 1]) {
            case 'A':
                atcgNumbers[0]++;
                break;
            case 'T':
                atcgNumbers[1]++;
                break;
            case 'C':
                atcgNumbers[2]++;
                break;
            case 'G':
                atcgNumbers[3]++;
                break;
            case 'N':
                atcgNumbers[0]++;
                atcgNumbers[1]++;
                atcgNumbers[2]++;
                atcgNumbers[3]++;
                break;
            default:
                break;
        }

        int cmp = (pos+1);
        int allowedMismatch = min(maxMismatch, cmp/allowOneMismatchForEach);

        bool needToBreak = true;
        for(int b=0; b<4; b++) {
            if(cmp - atcgNumbers[b] <= allowedMismatch)
                needToBreak = false;
        }
        if(needToBreak && (pos >= allowOneMismatchForEach || pos+1 >= compareReq-1)) {
            break;
        }
    }

    // has polyX
    if(pos+1 >= compareReq) {
        // find the poly
        int poly;
        int maxCount = -1;
        for(int b=0; b<4; b++) {
            if(atcgNumbers[b] > maxCount){
                maxCount = atcgNumbers[b];
                poly = b;
            }
        }
        char polyBase = ATCG_BASES[poly];
        while(data[rlen - pos - 1] != polyBase && pos>=0)
            pos--;

        r->resize(rlen - pos - 1);
        if(fr)
          fr->addPolyXTrimmed(poly, pos + 1);
    }
}

