#ifndef ___EDITDISTANCE__H__
#define ___EDITDISTANCE__H__

#include <stdint.h>
#include <string>

// struct PatternMap {
//     uint64_t p_[256][4];
//     unsigned int tmax_;
//     unsigned int tlen_;
// };

using namespace std;

unsigned int edit_distance(const char *a, const unsigned int asize, const char *b, const unsigned int bsize);
// void create_patternmap(struct PatternMap *pm, const int64_t *a, const unsigned int size);
// unsigned int edit_distance_by_patternmap(struct PatternMap *mp, const int64_t *b, const unsigned int size);

unsigned int edit_distance(string a, string b);

unsigned int hamming_distance(const char *a, const unsigned int asize, const char *b, const unsigned int bsize);

bool editdistance_test();

#endif
