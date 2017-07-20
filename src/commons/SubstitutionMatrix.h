#ifndef SUBSTITUTION_MATRIX_H
#define SUBSTITUTION_MATRIX_H

// Written by Maria Hauser mhauser@genzentrum.lmu.de
//  
// Represents a simple amino acid substitution matrix. 
// Can read a straight-forward substitution matrix file (.mat ending of the file) and probabilities file (.out ending of the file).
// If a probabilities file is given, it calculates a biased matrix (produces shorter = more precise alignments).
//

#include <cstddef>
#include "BaseMatrix.h"

class SubstitutionMatrix: public BaseMatrix {

    public:
        SubstitutionMatrix(const char *scoringMatrixFileName_, float bitFactor, float scoreBias);

        virtual ~SubstitutionMatrix();

        virtual float getBitFactor() {return bitFactor; }
    
        virtual double getBackgroundProb(size_t aa_index) { return pBack[aa_index]; }

        static void calcLocalAaBiasCorrection(const BaseMatrix *m ,const int *int_sequence, const int N, float *compositionBias);

        static void calcGlobalAaBiasCorrection( short *profileScores,
                                               const size_t profileAASize,
                                               const int N);
        bool estimateLambdaAndBackground(const double ** mat, int alphabetSize, double * pBack, double & lambda);


        void setupLetterMapping();


        struct FastMatrix{
            const char ** matrix;
            const char * matrixData;
            const size_t asciiStart;
            FastMatrix(const char ** matrix, const char * matrixData, const size_t asciiStart):
                    matrix(matrix), matrixData(matrixData), asciiStart(asciiStart)
            {}
        };

        // build matrix from ~ (=0) to ~(=122)
        static FastMatrix createAsciiSubMat(SubstitutionMatrix & submat){
            const size_t asciiStart = 0;
            const size_t asciiEnd = 'z';
            const size_t range = asciiEnd-asciiStart;
            char ** matrix = new char *[range];
            char * matrixData = new char[range*range];
            for(size_t i = 0; i < range; i++) {
                matrix[i] = matrixData+(i*range);
                int curr_i = submat.aa2int[asciiStart+(int)i];
                for (size_t j = 0; j < range; j++) {
                    int curr_j = submat.aa2int[asciiStart+(int)j];
                    matrix[i][j] = submat.subMatrix[curr_i][curr_j];
                }
            }
            return FastMatrix((const char**) matrix,
                              (const char*) matrixData,
                              asciiStart);
        }

private:

        const char* scoringMatrixFileName;

        int parseAlphabet(char * word, char * int2aa, int * aa2int);

        int readProbMatrix(const std::string &matrixData);

        float bitFactor;
};

#endif
