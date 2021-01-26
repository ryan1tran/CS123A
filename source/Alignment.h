/* CLASS FOR PERFORMING ALIGNMENTS */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "AlignmentPair.h"
#include "IO.h"
#include <random>

class Alignment
{
    private:
        std::vector<Sequence*> sequences;

    public:
        explicit Alignment(const std::string& filename);
        std::vector<Sequence*> getSequences();
        void setSequences(std::vector<Sequence*> sequences);
        void MSA();
        void Alignment::smithWatermanAlignment(std::string sequence1, std::string sequence2, int match, int mismatch, int singleGap, int extendedGap, AlignmentPair* ap);
};

#endif
