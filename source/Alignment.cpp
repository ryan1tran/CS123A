#include "Alignment.h"

// takes in file name and reads in the FASTA records
Alignment::Alignment(const std::string& filename)
{
    // parses the file and stores the returned list of Sequence*'s
    setSequences(readInput(filename));
}

std::vector<Sequence*> Alignment::getSequences()
{
    return this->sequences;
}

void Alignment::setSequences(std::vector<Sequence*> sequences)
{
    this->sequences = std::move(sequences);
}

// performs the multiple sequence alignment
void Alignment::MSA()
{
    // create two vectors with the Sequence*'s in them and another vector for AlignmentPair's
    std::vector<Sequence*> sequenceList = getSequences();
    std::vector<Sequence*> sequencesStored = getSequences();
    std::vector<AlignmentPair*> pairs;

    // number of sequences total; to be used much later
    int numSequences = sequenceList.size();

    // cool logo for fun

    std::cout << "\n"
                 " 888b     d888          888b     d888  .d8888b.        d8888 \n"
                 " 8888b   d8888          8888b   d8888 d88P  Y88b      d88888 \n"
                 " 88888b.d88888          88888b.d88888 Y88b.          d88P888 \n"
                 " 888Y88888P888 888  888 888Y88888P888  \"Y888b.      d88P 888 \n"
                 " 888 Y888P 888 888  888 888 Y888P 888     \"Y88b.   d88P  888 \n"
                 " 888  Y8P  888 888  888 888  Y8P  888       \"888  d88P   888 \n"
                 " 888   \"   888 Y88b 888 888   \"   888 Y88b  d88P d8888888888 \n"
                 " 888       888  \"Y88888 888       888  \"Y8888P\" d88P     888 \n"
                 "                    888                                      \n"
                 "               Y8b d88P                                      \n"
                 "                \"Y88P\"                                       \n" << std::endl;

    std::cout << " Welcome to MyMSA, the multiple sequence alignment algorithm that tries its best\n"
                 " to compete against its professional counterparts (with varying levels of success)!" << std::endl << std::endl;

    // explain expected inputs and provide additional information
    std::cout << " MyMSA is implemented using Smith-Waterman, a dynamic programming local alignment\n"
                 " algorithm that is similar to Needleman-Wunsch. As such, there are four required\n"
                 " parameters:\n"
                 "    1. a match score\n"
                 "    2. a mismatch penalty\n"
                 "    3. a single gap penalty\n"
                 "    4. an extended gap penalty\n"
                 "\n"
                 " The match score is recommended to be a positive value. The mismatch penalty is\n"
                 " recommended to be either zero or a negative value. However, the single gap penalty\n"
                 " must be STRICTLY less than the mismatch penalty. The extended gap penalty must also\n"
                 " be STRICTLY less than the gap penalty. The algorithm may crash if these are not followed.\n"
                 "\n"
                 " If you are unsure what inputs to use, (1, -1, -2, -4) is a good first option to try!" << std::endl;

    // variables to hold inputs
    int match;
    int mismatch;
    int singleGap;
    int extendedGap;

    // take in user inputs
    std::cout << "   Enter match score: ";
    std::cin >> match;
    std::cout << "   Enter mismatch penalty: ";
    std::cin >> mismatch;
    std::cout << "   Enter gap-opening penalty: ";
    std::cin >> singleGap;
    std::cout << "   Enter gap-extension penalty: ";
    std::cin >> extendedGap;

    // message to let users know that the inputs were successfully captured, as well as the number of sequences being aligned
    std::cout << std::endl << " Aligning " << numSequences << " sequences...";

    // create empty variables that will be used throughout the rest of this function
    std::string tempSequence1;
    std::string tempSequence2;

    // make AlignmentPairs of every two Sequence*'s in sequenceList (results in every unique pair of Sequence*)
    for (int i = 0; i < sequenceList.size() - 1; i++)
    {
        for (int j = i + 1; j < sequenceList.size(); j++)
        {
            auto* ap = new AlignmentPair(sequenceList.at(i), sequenceList.at(j));
            pairs.push_back(ap);

            tempSequence1 = sequenceList.at(i)->getSequence();
            tempSequence2 = sequenceList.at(j)->getSequence();

            // align the two Sequence*'s that are now in an AlignmentPair using the Smith-Waterman (local) alignment algorithm
            smithWatermanAlignment(tempSequence1, tempSequence2, match, mismatch, singleGap, extendedGap, ap);
        }
    }

    // addition to the previous cout that acts as a sort of progression bar
    std::cout << "...";

    // create variables in preparation for the scope of the while loop
    int highestScore;
    int index = 0;
    std::string consensusSequence;

    // until there are no more sequences left in sequenceList
    while (sequenceList.size() > 1)
    {
        // set the current highest score to 0
        highestScore = 0;

        // iterate over the AlignmentPair vector for the highest score
        for (int i = 0; i < pairs.size(); i++)
        {
            if (pairs.at(i)->getScore() > highestScore)
            {
                highestScore = pairs.at(i)->getScore();
                index = i; // keep track of the index of the highest score
            }
        }

        // get the two sequences from the AlignmentPair that has the highest score
        std::string refSequence1 = pairs.at(index)->getS1();
        std::string refSequence2 = pairs.at(index)->getS2();

        // get the two aligned sequences from the AlignmentPair that has the highest score;
        std::string alignedSequence1 = pairs.at(index)->getAlignedSequence1();
        std::string alignedSequence2 = pairs.at(index)->getAlignedSequence2();

        // the two alignedSequence strings above will be used to construct their respective consensus sequence
        consensusSequence = "";
        srand(time(nullptr));
        double rand;

        // for the entire length of alignedSequence1 (as well as alignedSequence2 since they are the same length)
        for (int i = 0; i < alignedSequence1.length(); i++)
        {
            // if the nucleotides match
            if (alignedSequence1.at(i) == alignedSequence2.at(i))
            {
                consensusSequence += alignedSequence1.at(i); // consensus adds that nucleotide
            }
            else if (alignedSequence1.at(i) != alignedSequence2.at(i)) // if the nucleotides do not match
            {
                // take the character that is not a '-'
                if (alignedSequence1.at(i) == '-')
                {
                    consensusSequence += alignedSequence2.at(i);
                }
                else if (alignedSequence2.at(i) == '-')
                {
                    consensusSequence += alignedSequence1.at(i);
                }
                else // in the case that either are gaps but are still not matching
                {
                    // use a random value to decide which nucleotide to take (since it doesn't matter which is taken in this case)
                    rand = std::rand();

                    if(rand <= 0.5)
                        consensusSequence += alignedSequence1.at(i);
                    else if(rand > 0.5)
                        consensusSequence += alignedSequence2.at(i);
                }
            }
        }

        // once the consensus sequence is complete, remove the original sequences that were just used from sequenceList

        // int to keep track of the position of the sequence to remove
        int toRemove = -1;

        // search for the first sequence
        for (int i = 0; i < sequenceList.size(); i++)
        {
            if (sequenceList.at(i)->getSequence() == refSequence1)
            {
                toRemove = i;
                break;
            }
        }

        // erase it if found
        if (toRemove != -1)
            sequenceList.erase(sequenceList.begin() + toRemove);

        // reset tracker to -1
        toRemove = -1;
        // search for the second sequence
        for (int i = 0; i < sequenceList.size(); i++)
        {
            if (sequenceList.at(i)->getSequence() == refSequence2)
            {
                toRemove = i;
                break;
            }
        }

        // erase it if found as well
        if (toRemove != -1)
            sequenceList.erase(sequenceList.begin() + toRemove);

        // remove the pairs with the highest scores
        for (int i = 0; i < pairs.size(); i++)
        {
            if (pairs.at(i)->getS1() == refSequence1 || pairs.at(i)->getS1() == refSequence2 ||
                pairs.at(i)->getS2() == refSequence1 || pairs.at(i)->getS2() == refSequence2)
            {
                pairs.erase(pairs.begin() + i);
            }
        }

        // using the consensus sequence previously generated, complete alignments against the remaining sequences in sequenceList
        for (auto & sequence : sequenceList)
        {
            auto* ap = new AlignmentPair(sequence, consensusSequence);
            pairs.push_back(ap);

            tempSequence1 = sequence->getSequence();

            smithWatermanAlignment(tempSequence1, consensusSequence, match, mismatch, singleGap, extendedGap, ap);
        }

        // add the consensus sequence to the sequenceList to be aligned against in further iterations
        sequenceList.push_back(new Sequence("Consensus", consensusSequence));
    }

    // additional progress bar once outside of the while loop
    std::cout << "...";

    // clear out the pairs vector in preparation for the final alignments
    pairs.clear();

    // lastly, align all the initial sequences in the untouched sequencesStored vector against the consensus sequence
    for (auto & i : sequencesStored)
    {
        auto* ap = new AlignmentPair(i, consensusSequence);
        pairs.push_back(ap);

        tempSequence1 = i->getSequence();

        smithWatermanAlignment(tempSequence1, consensusSequence, match, mismatch, singleGap, extendedGap, ap);
    }

    // last cout progression bar
    std::cout << "..." << std::endl << std::endl;

    // create file name for the output
    std::string DATA_FILE_NAME_OUT = "msa-output";

    // run the function, with all the relevant parameters, to write to a file
    std::string success = writeToOutput(DATA_FILE_NAME_OUT, pairs, consensusSequence, numSequences, match, mismatch, singleGap, extendedGap);

    // if the return string is not empty, it is successful and an informational message is displayed
    if (!success.empty())
        std::cout << " Success! Refer to \"" << success << "\" for results (located in the same directory as MyMSA.exe)." << std::endl;
    else // else print an error message
        std::cout << " Error occurred. Output file could not be created/opened." << std::endl;
}

// performs local alignment through the Smith-Waterman dynamic programming algorithm, similar to Needleman-Wunsch
void Alignment::smithWatermanAlignment(std::string sequence1, std::string sequence2, int match, int mismatch, int singleGap, int extendedGap, AlignmentPair* ap)
{
    // variables for the length of each sequence for ease of access
    int sequenceLength1 = sequence1.length();
    int sequenceLength2 = sequence2.length();

    // allocate two matrices:
        // scoring - for calculating an alignment score based on the input parameters
        // traceback - for keeping track of the alignment path
    int** scoringMatrix = new int*[sequenceLength2 + 1];
    for (int i = 0; i < sequenceLength2 + 1; ++i)
        scoringMatrix[i] = new int[sequenceLength1 + 1];

    int** tracebackMatrix = new int*[sequenceLength2 + 1];
    for (int i = 0; i < sequenceLength2 + 1; ++i)
        tracebackMatrix[i] = new int[sequenceLength1 + 1];

    // parameters used in the tracebackMatrix to keep track of which direction the alignment is going
    int DIAGONAL = 1;
    int VERTICAL = 2;
    int HORIZONTAL = 3;
    int VERTICAL_DIAGONAL_TIE = 4;
    int HORIZONTAL_DIAGONAL_TIE = 5;
    int VERTICAL_HORIZONTAL_TIE = 6;
    int ALL_TIE = 7;

    // fill the matrices with values to avoid garbage values
    for (int i = 0; i < sequenceLength2 + 1; i++)
    {
        for (int j = 0; j < sequenceLength1 + 1; j++)
        {
            scoringMatrix[i][j] = 0;

            // specific values based on the algorithm
            if (i == 0 && j >= 1)
                tracebackMatrix[i][j] = 3;
            else if (i >= 1 && j == 0)
                tracebackMatrix[i][j] = 2;
        }
    }

    /* calculating the scoring matrix */

    // variables for keeping track of current scores and choosing the maximum out of them, again similar to Needleman-Wunsch
    int diagonalScore;
    int horizontalScore;
    int verticalScore;

    // calculate the scoring matrix starting at (1, 1)
    for (int i = 1; i < sequenceLength2 + 1; i++)
    {
        for (int j = 1; j < sequenceLength1 + 1; j++)
        {
            // if there is a match, set all the scores to the match value
            if (sequence1.at(j - 1) == sequence2.at(i - 1))
            {
                diagonalScore = match;
                horizontalScore = match;
                verticalScore = match;
            } else { // else set all the scores to mismatch
                diagonalScore = mismatch;
                horizontalScore = mismatch;
                verticalScore = mismatch;
            }

            // the value of diagonalScore = (match or mismatch value) + the score from the previous diagonal (i - 1, j - 1)
            diagonalScore += (scoringMatrix[i - 1][j - 1]);

            // calculate the highest score from the left position (i, j - 1) to the current position (i, j)
            // the value of horizontalScore = (match or mismatch value) + the gap penalty for moving horizontally
            if (i == sequenceLength2)
                horizontalScore += scoringMatrix[i][j - 1];
            else if (tracebackMatrix[i][j-1] == 3 || tracebackMatrix[i][j-1] == 5 || tracebackMatrix[i][j-1] == 6 || tracebackMatrix[i][j-1] == 7)
                horizontalScore += extendedGap + scoringMatrix[i][j - 1];
            else
                horizontalScore += (singleGap + scoringMatrix[i][j - 1]);

            // calculate the highest score from the top position (i - 1, j) to the current position (i, j)
            // the value of verticalScore = (match or mismatch value) + the gap penalty for moving vertically
            if (j == sequenceLength1)
                verticalScore += scoringMatrix[i - 1][j];
            else if (tracebackMatrix[i-1][j] == 2 || tracebackMatrix[i-1][j] == 4 || tracebackMatrix[i-1][j] == 6 || tracebackMatrix[i-1][j] == 7)
                verticalScore += (extendedGap + scoringMatrix[i - 1][j]);
            else
                verticalScore += (singleGap + scoringMatrix[i - 1][j]);

            // set highestScore to the maximum of diagonalScore, horizontalScore, and verticalScore
            int highestScore = std::max(diagonalScore, std::max(horizontalScore, verticalScore));
            scoringMatrix[i][j] = highestScore; // set the current position of scoringMatrix to highestValue

            /* calculating the traceback matrix */

            // assigns values to the traceback matrix depending on if there are ties in scores or not
            if (diagonalScore > horizontalScore && diagonalScore > verticalScore)
                tracebackMatrix[i][j] = DIAGONAL;
            else if (verticalScore > diagonalScore && verticalScore > horizontalScore)
                tracebackMatrix[i][j] = VERTICAL;
            else if (horizontalScore > diagonalScore && horizontalScore > verticalScore)
                tracebackMatrix[i][j] = HORIZONTAL;
            else if (diagonalScore == verticalScore && diagonalScore > horizontalScore)
                tracebackMatrix[i][j] = VERTICAL_DIAGONAL_TIE;
            else if (diagonalScore == horizontalScore && diagonalScore > verticalScore)
                tracebackMatrix[i][j] = HORIZONTAL_DIAGONAL_TIE;
            else if (horizontalScore == verticalScore && horizontalScore > diagonalScore)
                tracebackMatrix[i][j] = VERTICAL_HORIZONTAL_TIE;
            else if (diagonalScore == verticalScore && diagonalScore == horizontalScore)
                tracebackMatrix[i][j] = ALL_TIE;
        }
    }

    // set the alignment score as the bottom-rightmost element in the scoring matrix
    int finalScore = scoringMatrix[sequenceLength2][sequenceLength1];

    /* using the traceback matrix to calculate the aligned sequences */

    // initial variable setup for iterating and containing alignment information
    int j = sequenceLength1;
    int i = sequenceLength2;
    std::string alignedSequence1;
    std::string alignedSequence2;
    srand(time(nullptr));
    double rand;

    // goes from the bottom-right element to either the top or left edge of the matrix
    while (i != 0 || j != 0)
    {
        // coming from diagonal, so prepend both nucleotides (since starting from bottom-right element) and decrement both counters
        if (tracebackMatrix[i][j] == DIAGONAL)
        {
            alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
            alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
            i--;
            j--;
        }
        // coming from the top, so add a gap onto alignedSequence1 and decrement i
        else if (tracebackMatrix[i][j] == VERTICAL)
        {
            alignedSequence1 = "-" + alignedSequence1;
            alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
            i--;
        }
        // coming from the left, so add a gap to alignedSequence2 and decrement j
        else if (tracebackMatrix[i][j] == HORIZONTAL)
        {
            alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
            alignedSequence2 = "-" + alignedSequence2;
            j--;
        }
        // tied in score for coming from top and diagonal
        else if (tracebackMatrix[i][j] == VERTICAL_DIAGONAL_TIE)
        {
            // resolve the tie using a random value (since it doesn't matter which nucleotide is taken in these tie cases)
            rand = std::rand();

            if (rand >= 0.5) // top case
            {
                alignedSequence1 = "-" + alignedSequence1;
                alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                i--;
            } else { // diagonal case
                alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                i--;
                j--;
            }
        }
        // tied in score for coming from left and diagonal
        else if(tracebackMatrix[i][j] == HORIZONTAL_DIAGONAL_TIE)
        {
            rand = std::rand();

            if (rand >= 0.5) // left case
            {
                alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                alignedSequence2 = "-" + alignedSequence2;
                j--;
            } else { // diagonal case
                alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                i--;
                j--;
            }
        }
        // tied in score for coming from top and left
        else if (tracebackMatrix[i][j] == VERTICAL_HORIZONTAL_TIE)
        {
            rand = std::rand();

            if (rand >= 0.5) // top case
            {
                alignedSequence1 = "-" + alignedSequence1;
                alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                i--;
            } else { // left case
                alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                alignedSequence2 = "-" + alignedSequence2;
                j--;
            }
        }
        // tied in score for coming from diagonal, left, and top
        else if (tracebackMatrix[i][j] == ALL_TIE)
        {
            // starting element case takes the diagonal nucleotide
            if (i == 1 && j == 1)
            {
                alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                i--;
                j--;
            }
            else // otherwise if it is not the starting element, decide by a random value again
            {
                rand = std::rand();

                if (rand <= 0.33) // diagonal case
                {
                    alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                    alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                    i--;
                    j--;
                }
                else if (rand > 0.33 && rand < 0.66) // left case
                {
                    alignedSequence1 = sequence1.at(j - 1) + alignedSequence1;
                    alignedSequence2 = "-" + alignedSequence2;
                    j--;
                }
                else if (rand >= 0.66) // top case
                {
                    alignedSequence1 = "-" + alignedSequence1;
                    alignedSequence2 = sequence2.at(i - 1) + alignedSequence2;
                    i--;
                }
            }
        }
    }

    // store final results and data in the AlignmentPair parameter
    ap->setScore(finalScore);
    ap->setAlignedSequence1(alignedSequence1);
    ap->setAlignedSequence2(alignedSequence2);
}
