#include "candidate.hpp"
#include "signalFirstInternal.h"
#include <vector>
#include "config.h"

int main(int argc, char** argv)
{
	sgcu::argumentParser<char> myP;
	myP.desc_.add_options()
		("help", "produce help message")
		("input", boost::program_options::value<std::string>(), "set input file")
		("dec", boost::program_options::value<double>(), "set declination of the center of the source in degrees")
		("ra", boost::program_options::value<double>(), "set right ascension of the center of the source in degrees")
		("decsize", boost::program_options::value<double>(), "set half of declination size of the source in degrees")
		("radius", boost::program_options::value<double>(), "set radius or half of right ascension size of the source in degrees")
		("number", boost::program_options::value<unsigned int>(), "set the number of events within the source")
		("seed", boost::program_options::value<double>(), "seed")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));
	double dec = myP.parseVariable("dec", 15.0);
	double ra = myP.parseVariable("ra", 260.0);
	double radius = myP.parseVariable("radius", 8.0);
	double decSize = myP.parseVariable("decsize", 0.0);
	unsigned int numberSourceEvents = myP.parseVariable("number", static_cast<unsigned int>(80) );
	unsigned int seed = myP.parseVariable("seed", time(NULL));
	srand(seed);

	std::string outputPath = "output/scenarioRandomWSource"+boost::lexical_cast<std::string>(seed)+"/";
	if (system(std::string("mkdir -p " + outputPath).c_str()) != 0) return 1;

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);

	// use randomized events
	std::vector<SG::candidate_t> randomizedCandidates;
	getRandomEventsFromEvents(allCandidates, randomizedCandidates);

	// add an artificial source
	const double oneDeg = M_PI/180.0;
	if (decSize == 0.0)
	{
		addSourceGaussian(randomizedCandidates, ra*oneDeg, dec*oneDeg, radius*oneDeg, numberSourceEvents);
	}
	else
	{
		addSource(randomizedCandidates, ra*oneDeg, dec*oneDeg, radius*oneDeg, decSize*oneDeg, numberSourceEvents);
	}
	
	// set up an object to handle all computations with a simple interface to this module
	SG::signalFirstInternal intern;
	intern.init(SG::fine);

	intern.setOutputFolder(outputPath);

	// Evaluate the clustering/density of the found candidates here
	intern.evaluateEventsToSkymap(randomizedCandidates, SG::fine);

	return 0;
}
