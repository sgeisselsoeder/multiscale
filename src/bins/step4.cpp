#include "signalFirstInternal.h"
#include "config.h"

int main()
{
	// set up an object to handle all computations with a simple interface to this module
	SG::signalFirstInternal intern;
	intern.init(SG::fine);

	std::string outputPath = "intermediate/step4/normalizationDataStep3/";
	if (system(std::string("mkdir -p " + outputPath).c_str()) != 0) return 1;

	// Evaluate the clustering/density of the found candidates here
	intern.renormalizeAllPseudoexperimentsAllClustersAllMetrics();

	return 0;
}
