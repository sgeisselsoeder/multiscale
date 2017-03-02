#include "candidate.hpp"
#include "signalFirstInternal.h"
#include <vector>
#include "config.h"
#include <algorithm>
#include "sgcu/sgcu.hpp"

int main(int argc, char** argv)
{
	sgcu::argumentParser<char> myP;
	// Declare the supported options.
	myP.desc_.add_options()
		("help", "produce help message")
		("input", boost::program_options::value<std::string>(), "set input file")
		("seed", boost::program_options::value<double>(), "seed")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));
	unsigned int seed = myP.parseVariable("seed", time(NULL));
	srand(seed);

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);

	std::random_shuffle(allCandidates.begin(), allCandidates.end());

	dumpEvents(allCandidates, inputFile+"Shuffled.txt");

	return 0;
}
