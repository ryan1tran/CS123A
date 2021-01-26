#include "AlignmentPair.h"

// constructor for when two Sequence*'s are passed in
AlignmentPair::AlignmentPair(Sequence* s1, Sequence* s2)
{
    this->s1 = s1->getSequence();
    this->s2 = s2->getSequence();
    this->name = s1->getName();
    score = 0;
    alignedSequence1 = "";
    alignedSequence2 = "";
}

// constructor for when a Sequence* and a string are passed in
AlignmentPair::AlignmentPair(Sequence* seq, std::string s)
{
    this->s1 = seq->getSequence();
    this->s2 = std::move(s);
    this->name = seq->getName();
    score = 0;
    alignedSequence1 = "";
    alignedSequence2 = "";
}

/* standard getters and setters below */

const std::string &AlignmentPair::getS1() const
{
    return s1;
}

void AlignmentPair::setS1(const std::string &s1)
{
    this->s1 = s1;
}

const std::string &AlignmentPair::getS2() const
{
    return s2;
}

void AlignmentPair::setS2(const std::string &s2)
{
    this->s2 = s2;
}

const std::string &AlignmentPair::getName() const
{
    return name;
}

void AlignmentPair::setName(const std::string &name)
{
    this->name = name;
}

int AlignmentPair::getScore() const {
    return score;
}

void AlignmentPair::setScore(int score) {
    AlignmentPair::score = score;
}

const std::string &AlignmentPair::getAlignedSequence1() const
{
    return alignedSequence1;
}

void AlignmentPair::setAlignedSequence1(const std::string &alignedSequence1)
{
    this->alignedSequence1 = alignedSequence1;
}

const std::string &AlignmentPair::getAlignedSequence2() const
{
    return alignedSequence2;
}

void AlignmentPair::setAlignedSequence2(const std::string &alignedSequence2)
{
    this->alignedSequence2 = alignedSequence2;
}
