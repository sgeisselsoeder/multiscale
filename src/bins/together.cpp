#include "signalFirstInternal.h"
#include <vector>
#include "config.h"

int main()
{
	srand(123);
//	srand(time(NULL));

	std::vector<SG::candidate_t> allCandidates;
	readCandidatesFromFile("./input/events.txt", allCandidates);

	// set up an object to handle all computations with a simple interface to this module
	SG::signalFirstInternal intern;
	intern.init(SG::fine);

// step1
	// Evaluate the clustering/density of the found candidates here
	intern.computeNormalizationsFromDataStep1(allCandidates, SG::fine);

	// it causes problems if the same scenario is used over and over again for renormalization
	srand(time(NULL));

// step2
	const unsigned int numRepetitions = 1;
	// Evaluate the clustering/density of the found candidates here
	intern.computeNormalizationsFromDataStep2(allCandidates, SG::fine, numRepetitions);

	srand(123);
//	srand(time(NULL));

// step4
	intern.renormalizeAllPseudoexperimentsAllClustersAllMetrics();

// step5
	intern.performEvaluate(allCandidates, SG::fine);

	return 0;
}
