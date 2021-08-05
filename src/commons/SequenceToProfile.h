//
// Created by Charlotte Tumescheit on 2021/06/08.
//

#ifndef TEST_KERAS_SEQUENCETOPROFILE_H
#define TEST_KERAS_SEQUENCETOPROFILE_H

#include <vector>
#include <fdeep/fdeep.hpp>

class SequenceToProfile {
public:
    SequenceToProfile(int seqType, size_t maxLen);
    ~SequenceToProfile();

    // Use NN to generate profile for each 13-mer
    void sequenceToProfile(const char *seq, unsigned int seqLen);
    // turn middle position into a kmer, add padding if needed
    std::vector<float> determineKmer(const char *seq, unsigned int seqLen, int i);

    float *profile;

private:
    int seqType;
    int kmerSize;
    fdeep::model model;
    std::vector<float> fkmer;

};



#endif //TEST_KERAS_SEQUENCETOPROFILE_H