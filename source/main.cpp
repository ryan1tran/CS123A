/* MAIN DRIVER CODE */

#include "Sequence.h"
#include "Alignment.h"

// main function that takes in command line inputs
int main(int argc, char * argv[])
{
    if (argc < 2) // if there are less than two arguments, exit without error
    {
        std::cout << "Not enough arguments." << std::endl;
        exit(0);
    } else if (argc > 2) { // if there are too many arguments, exit without error
        std::cout << "Too many arguments." << std::endl;
        exit(0);
    } else { // if exactly two arguments, continue

        // create and run a multiple sequence alignment using Alignment and its .MSA() function
        Alignment alignment(argv[1]);
        alignment.MSA();
    }

    return 0;
}
