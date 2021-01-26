/* HEADER FILE FOR INPUT/OUTPUT FUNCTIONS */

#ifndef IO_H
#define IO_H

#include "AlignmentPair.h"
#include <fstream>
#include <iomanip>

std::vector<Sequence*> readInput(const std::string& filename);
std::string writeToOutput(std::string& file_name, const std::vector<AlignmentPair*>& pairs, std::string consensusSequence, int numSequences, int match, int mismatch, int singleGap, int extendedGap);
#endif
