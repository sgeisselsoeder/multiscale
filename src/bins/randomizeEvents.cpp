#include "sgcu/sgcu.hpp"
#include "candidate.hpp"
#include "signalFirstInternal.h"
#include <vector>
#include "config.h"

int main(int argc, char** argv)
{
	sgcu::argumentParser<char> myP;
	// Declare the supported options.
	myP.desc_.add_options()
		("help", "produce help message")
		("input", boost::program_options::value<std::string>(), "set input file")
		("seed", boost::program_options::value<unsigned int>(), "seed")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));
	unsigned int seed = myP.parseVariable("seed", static_cast<unsigned int>( time(NULL) ) );
	srand(seed);

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);

	// use randomized events
	std::vector<SG::candidate_t> randomizedCandidates;
	getRandomEventsFromEvents(allCandidates, randomizedCandidates);

	dumpEvents(randomizedCandidates, inputFile+"Randomized.txt");

	return 0;
}
