#include "candidate.hpp"
#include "signalFirstInternal.h"
#include <vector>
#include "config.h"
#include "sgcu/sgcu.hpp"
#include <utility>

int main(int argc, char** argv)
{
	sgcu::argumentParser<char> myP;
	// Declare the supported options.
	myP.desc_.add_options()
		("help", "produce help message")
		("input", boost::program_options::value<std::string>(), "set input file")
		("splitRatio", boost::program_options::value<double>(), "set split ratio")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));
	double splitRatio = myP.parseVariable("splitRatio", double(0.5) );

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);
	unsigned int splitIndex = static_cast<unsigned int>( static_cast<double>(allCandidates.size())*splitRatio + 0.5 );
	assert(splitIndex < allCandidates.size());
	std::vector<SG::candidate_t> split1( allCandidates.begin(), allCandidates.begin()+splitIndex );
	std::vector<SG::candidate_t> split2( allCandidates.begin()+splitIndex, allCandidates.end() );

	dumpEvents(split1, inputFile+"SplitPart1.txt");
	dumpEvents(split2, inputFile+"SplitPart2.txt");

	return 0;
}

