#include "signalFirstInternal.h"
#include <vector>
#include "config.h"
#include "sgcu/sgcu.hpp"

int main(int argc, char** argv)
{
	sgcu::argumentParser<char> myP;
	// Declare the supported options.
	myP.desc_.add_options()
		("help", "produce help message")
		("input", boost::program_options::value<std::string>(), "set input file")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));

	std::string outputPath = "intermediate/step1/normalizationDataStep1/";
	if (system(std::string("mkdir -p " + outputPath).c_str()) != 0) return 1;

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);

	// set up an object to handle all computations with a simple interface to this module
	SG::signalFirstInternal intern;
	intern.init(SG::fine);

	// Evaluate the clustering/density of the found candidates here
	intern.computeNormalizationsFromDataStep1(allCandidates, SG::fine);

	return 0;
}
