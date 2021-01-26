#include "IO.h"

// reads the file and returns a vector of Sequence*'s
std::vector<Sequence*> readInput(const std::string& filename)
{
    // vector of Sequence*'s that will hold the inputs and be returned at the end of the function
    std::vector<Sequence*> inputs;

    // open the input file
    std::ifstream source;
    source.open(filename);

    if (!source) // if the file does not open successfully, exit with error
        exit(1);

    std::vector<std::string> lines;
    std::string line;
    std::vector<int> nameLines;
    int index = 0;

    // read through the whole file
    while (std::getline(source, line))
    {
        if (line.at(0) == '>')
            nameLines.push_back(index); // keep track of lines that start with '>'

        lines.push_back(line); // store all lines
        index++;
    }

    std::string name, sequence;

    // separate the lines into names and sequences, then store them as a Sequence* in the return vector
    for (int i = 0; i < nameLines.size(); i++)
    {
        int pos = nameLines.at(i);
        name = lines.at(pos);

        if (i == nameLines.size() - 1) // last FASTA record in the file
        {
            for (int j = pos + 1; j < lines.size(); j++)
                sequence += lines.at(j);

            inputs.push_back(new Sequence(name, sequence));

        } else { // anything but the last FASTA record in the file

            for (int j = pos + 1; j < nameLines.at(i + 1); j++)
                sequence += lines.at(j);

            inputs.push_back(new Sequence(name, sequence));
            name.clear();
            sequence.clear();
        }
    }

    // return vector of input Sequence*'s
    return inputs;
}

// writes results of MSA to an output file; returns the name of the file if successful
std::string writeToOutput(std::string& file_name, const std::vector<AlignmentPair*>& pairs, std::string consensusSequence, int numSequences, int match, int mismatch, int singleGap, int extendedGap)
{
    // get the current date
    time_t now = std::time(nullptr);
    tm *ltm = localtime(&now);
    std::string month = std::to_string(1 + ltm->tm_mon);
    std::string day = std::to_string(ltm->tm_mday);
    std::string year = std::to_string(1900 + ltm->tm_year);

    // add the date, a scalar, and ".txt" to the end of the file name (to avoid overwriting output files)
    int scalar = 1;
    file_name += "_" + month + "-" + day + "-" + year + "_" + std::to_string(scalar) + ".txt";

    // while the file name already exists
    while (std::ifstream(file_name))
    {
        // edit the name by incrementing the scalar number that comes before ".txt"
        int fileLength = file_name.length();
        int intLength = std::to_string(scalar).length();
        file_name.erase(fileLength - intLength - 4, intLength + 4);

        file_name += std::to_string(scalar++) + ".txt";
    }

    // open the file
    std::ofstream new_file(file_name);

    if (!new_file) // if the file does not open, return an empty string
        return "";


    // writes sequence and input information
    new_file << "DETAILS & PARAMETERS" << std::endl;
    new_file << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    new_file << "Number of FASTA sequences: " << numSequences << std::endl;
    new_file << "Match score: " << match
           << "\nMismatch penalty: " << mismatch
           << "\nSingle-gap penalty: " << singleGap
           << "\nExtended-gap penalty: " << extendedGap << std::endl << std::endl << std::endl;


    // find the max alignment length and max FASTA name length (accession number only) for future formatting
    int maxAlignmentLength = 0;
    int maxNameLength = 0;
    int maxLength;

    for (auto & pair : pairs)
        if(pair->getAlignedSequence1().length() > maxAlignmentLength)
            maxAlignmentLength = pair->getAlignedSequence1().length();

    for (auto & pair : pairs)
    {
        int space = pair->getName().find(' ');
        pair->setName(pair->getName().substr(1, space)); // keep only the accession numbers to save space

        if (pair->getName().length() > maxNameLength)
            maxNameLength = pair->getName().length();
    }

    maxLength = maxAlignmentLength + maxNameLength;

    // writes headers for the alignment output; also includes a tip regarding how the output is displayed
    new_file << std::left << std::setw(maxNameLength + 1) << "ACCESSION NUMBER" << "ALIGNMENT (Turn off word wrap if the alignment does not fit in the window.)" << std::endl;

    // underlines for more clarity
    for (int i = 0; i < maxLength + 2; i++)
        new_file << "~";
    new_file << std::endl;

    // writes the names of the sequences and their respective alignments
    for (auto & pair : pairs)
        new_file << std::left << std::setw(maxNameLength) << pair->getName() << " " << pair->getAlignedSequence1() << std::endl;

    // line to separate sequence alignments from the consensus sequence
    for (int i = 0; i < maxLength + 2; i++)
        new_file << "=";

    // writes the consensus sequence
    new_file << std::endl << std::left << std::setw(maxNameLength + 1) << "" << consensusSequence << std::endl;
    new_file << std::left << std::setw(maxNameLength + 1) << "";

    // writes an asterisk if the nucleotide at a particular position is the same for all sequences; writes a space if not all matching
    bool matching;
    char nucleotide;

    for (int i = 0; i < consensusSequence.length(); i++)
    {
        // reset matching to true at the beginning of every loop iteration
        matching = true;
        nucleotide = consensusSequence.at(i);

        for (auto &pair : pairs)
        {
            // if any nucleotide does not match
            if (pair->getAlignedSequence1().at(i) != nucleotide)
            {

                matching = false; // set matching as false
                break;            // and immediately break
            }
        }

        // write the respective icon depending on the value of matching
        if (matching)
            new_file << "*";
        else
            new_file << " ";
    }

    // close the file
    new_file.close();

    // return the file name to be used later and to signify a successful write
    return file_name;
}
