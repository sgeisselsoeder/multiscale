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
		("seed", boost::program_options::value<double>(), "seed")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));
	unsigned int seed = myP.parseVariable("seed", time(NULL));
	srand(seed);

	std::string outputPath = "output/scenarioEvaluation/";
	if (system(std::string("mkdir -p " + outputPath).c_str()) != 0) return 1;

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);

	// set up an object to handle all computations with a simple interface to this module
	SG::signalFirstInternal intern;
	intern.init(SG::fine);

	intern.setOutputFolder(outputPath);
	// Evaluate the clustering/density of the found candidates here
	intern.evaluateEventsToSkymap(allCandidates, SG::fine);

	return 0;
}
