/* CLASS FOR HOLDING FASTA DATA */

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <iostream>

class Sequence
{
    private:
        std::string name;
        std::string sequence;

    public:
        Sequence(std::string name, std::string sequence);
        std::string getName();
        std::string getSequence();
};

#endif
