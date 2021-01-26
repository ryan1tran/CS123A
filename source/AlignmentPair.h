/* CLASS FOR HOLDING TWO SEQUENCES AND RELEVANT DATA */

#ifndef ALIGNMENTPAIR_H
#define ALIGNMENTPAIR_H

#include "Sequence.h"
#include <ctime>

class AlignmentPair
{
    private:
        std::string s1;
        std::string s2;
        std::string name;
        int score;
        std::string alignedSequence1;
        std::string alignedSequence2;

    public:
        AlignmentPair(Sequence* s1, Sequence* s2);
        AlignmentPair(Sequence* seq, std::string s);
        const std::string &getS1() const;
        void setS1(const std::string &s1);
        const std::string &getS2() const;
        void setS2(const std::string &s2);
        const std::string &getName() const;
        void setName(const std::string &name);
        int getScore() const;
        void setScore(int score);
        const std::string &getAlignedSequence1() const;
        void setAlignedSequence1(const std::string &alignedSequence1);
        const std::string &getAlignedSequence2() const;
        void setAlignedSequence2(const std::string &alignedSequence2);
};


#endif
