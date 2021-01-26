#include "Sequence.h"

// constructor with both parameters provided since this is the only expected case
Sequence::Sequence(std::string name, std::string sequence)
{
    this->name = std::move(name);
    this->sequence = std::move(sequence);
}

/* standard getters and setters below */

std::string Sequence::getName()
{
    return this->name;
}

std::string Sequence::getSequence()
{
    return this->sequence;
}
