#include "candidate.hpp"
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
		("ra", boost::program_options::value<double>(), "set right ascension of the center of the source in degrees")
		("dec", boost::program_options::value<double>(), "set declination of the center of the source in degrees")
		("radius", boost::program_options::value<double>(), "set radius or half of right ascension size of the source in degrees")
		("decsize", boost::program_options::value<double>(), "set half of declination size of the source in degrees")
		("number", boost::program_options::value<unsigned int>(), "set the number of events within the source")
		("templatefile", boost::program_options::value<std::string>(), "use a template from a file (image) to define the source")
		;
	myP.defaultParse(argc, argv);
	std::string inputFile = myP.parseVariable("input", std::string("."));
	double dec = myP.parseVariable("dec", 15.0);
	double ra = myP.parseVariable("ra", 260.0);
	double radius = myP.parseVariable("radius", 8.0);
	double decSize = myP.parseVariable("decsize", 0.0);
	unsigned int numberSourceEvents = myP.parseVariable("number", static_cast<unsigned int>(80) );
	std::string templateFile = myP.parseVariable("templatefile", std::string(""));

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile(inputFile, allCandidates);

	// add an artificial source
	const double oneDeg = M_PI/180.0;
	if (templateFile == "")
	{
		if (decSize == 0.0)
		{
			addSourceGaussian(allCandidates, ra*oneDeg, dec*oneDeg, radius*oneDeg, numberSourceEvents);
		}
		else
		{
			addSource(allCandidates, ra*oneDeg, dec*oneDeg, radius*oneDeg, decSize*oneDeg, numberSourceEvents);
		}
	}
	else
	{
		SG::sphericalGrid<double> sourceShape;
		// TODO: read image and set source shape accordingly
		addSourceTemplate(allCandidates, sourceShape, numberSourceEvents);
	}
	dumpEvents(allCandidates, inputFile+"With"+boost::lexical_cast<std::string>(numberSourceEvents)+"SourceEvents.txt");
	dumpEvents(allCandidates, inputFile+"WithSourceEvents.txt");

	return 0;
}
