/**
 * sgeisselsoeder
 *
 * June 2016
 *
 */
	
#ifndef SGALLSCALES_H
#define SGALLSCALES_H

#include <vector>
#include <string>

#include "candidate.hpp"
#include "sphericalGrid.hpp"

namespace SG{

const double maximumDistanceToEvaluateInDegree = 90.0;

struct allScales{
	unsigned int numBins_;
	double stepDegree_;
	double stepSize_;
	
	// the current evaluations for each searchpoint and each distance bin
	std::vector<std::vector<double> > aroundRatioAll;
	
	unsigned int size() const;

	allScales();

	void computePoissonProbabilityForAllRatios(unsigned int eventSize);
	void computePoissonProbabilityForAllRatios(SG::allScales const & normalization, unsigned int eventSize);
	int loadNormalization2(unsigned int gridsize, unsigned int numBins, unsigned int numEvents, std::string folder);
	int loadFileToData(std::string file, unsigned int gridsize, unsigned int numBins, std::string searchStringEnding, std::vector<std::vector<double> >& data);
	int dumpRatiosAsNormalization(std::string folder, unsigned int numEvents);
	
	long double computePValue(long double mean, long double found);
	
	void clear();
	void init(double stepDeg, unsigned int numGridPoints = 0);
	void set(double stepDeg, unsigned int numGridPoints = 0);
	
	void computeAllRatios(const SG::sphericalGrid<double>& grid, const std::vector<SG::candidate_t>& events);
	int dumpNormalizationToDisk(std::string folder = ".");
	void evaluateUnmergedNormalizations();
	
	int dumpNormalizationResultsToDisk(std::string folder = ".");
	int getNormalizationResults(std::string file, unsigned int gridsize, unsigned int numBins);
	
	void compute(const SG::sphericalGrid<double>& grid, unsigned int currentGridIndex, const std::vector<SG::candidate_t>& events);
	int dumpToDisk();
	
	int dumpAllRatiosToDisk();
	int loadAllRatios(std::string file, unsigned int gridsize, unsigned int numBins);
	
	void segment(unsigned int numRepetitions = 1);
	
	void findClusters(std::vector<std::vector<unsigned int> >& clusters) const;
	
//	void selfCheck();

	void addRatios(struct allScales const & b);
	void divideRatios(double value);

	void assertNoNansAndInfs() const;
};

void segmentAllRatiosScalewise(allScales & scales);

void setToZero(allScales & scales);

unsigned int countNumberOfNonZeros(struct allScales const & allScales);

} // end of namespace SG

#endif 
