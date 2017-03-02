/**
 * class: signalFirstInternal.h
 * Separating the internals from the Seatray module
 *
 * sgeisselsoeder
 *
 * Nov 2013
 *
 * (c) 2013 Antares Collaboration
 */
	
#ifndef SignalFirstInternal_H
#define SignalFirstInternal_H

#include <vector>
#include <string>
#include <set>
#include <utility>

#include "candidate.hpp"
#include "sphericalGrid.hpp"
#include "allScales.h"
#include "vtkWriter.h"
#include "sgcu/histogram.hpp"
#include "clusters.h"
#include "rememberBoostHistogram.h"

namespace SG{

void addSource(std::vector<SG::candidate_t>& artificial, double az, double zen, double delta, unsigned int num);
void addSource(std::vector<SG::candidate_t>& artificial, double az, double zen, double deltaAz, double deltaZen, unsigned int num);
void addSourceGaussian(std::vector<SG::candidate_t>& artificial, double az, double zen, double sigma, unsigned int num);
void addSourceTemplate(std::vector<SG::candidate_t>& artificial, SG::sphericalGrid<double> const & sourceShape, unsigned int num);
void dumpEvents(std::vector<SG::candidate_t> const & events, std::string outputFile);

void getRandomEventsFromEvents(const std::vector<SG::candidate_t>& events, std::vector<SG::candidate_t>& eventsRandomized);
void addRandomEventsFromEvents(const std::vector<SG::candidate_t>& events, std::vector<SG::candidate_t>& eventsRandomized);
void sumScalesToGrid(SG::allScales const & scales, SG::sphericalGrid<double> & grid);
void invertRatiosAndTakeLog10(allScales & scales);

void writeRatiosToGrid(allScales const & scales, unsigned int distanceIndex, SG::sphericalGrid<double> & grid);
void washOutWithinScales(allScales & scales, SG::sphericalGrid<double> const & grid);
void washOutBetweenScales(allScales & scales);


/**
 * @brief 
 */
class signalFirstInternal
{
public:
	/**
	 * Builds an instance of this class
	 */
	signalFirstInternal() : stepDegree_(0.5)
	{
	}
  
	/**
	 * Destroys an instance of this class
	 */
	~signalFirstInternal(){
		// NOTHING HERE
	}
	
	void init(unsigned int fine = 20);

	void computeNormalizationsFromDataStep1(const std::vector<SG::candidate_t>& events, unsigned int fine);
	void computeNormalizationsFromDataStep2(const std::vector<SG::candidate_t>& events, unsigned int fine, unsigned int numberOfRepetitions);
	void renormalizeAllPseudoexperimentsAllClustersAllMetrics();
	void performEvaluate(std::vector<SG::candidate_t>& events, unsigned int fine);
	void readStencil(unsigned int fine, std::string fileName = "stencilGrid.grid", double significanceThreshold = 0.0);
	void setOutputFolder(std::string newOutputFolder);
	void evaluateEventsToSkymap(const std::vector<SG::candidate_t>& events, unsigned int fine);


private:
	SG::sphericalGrid<double> grid_;
	SG::sphericalGrid<double> stencilGrid_;
	std::string outputFileName_;
	double stepDegree_;
	allScales ratios_;
	allScales normalization_;
	vtkWriter myVtkWriter_;
	std::string outputFolder_;
	
	std::vector<std::vector<cluster> > finalClustersMulti_;
	std::vector<double> segmentationPoints_;
	
	// default, assignment, and copy constructor declared private
	signalFirstInternal(const signalFirstInternal&);
	signalFirstInternal& operator=(const signalFirstInternal&);
	
	void computeOneNormalization(const std::vector<SG::candidate_t>& events, unsigned int fine);
	void writeNormalizationGrid(unsigned int distanceBinToPlot, std::string const & filenameWithoutNumberOrEnding);
	void writeNormalizationGrids(std::string folderPrefix);

	void evaluateDensity(const std::vector<SG::candidate_t>& events);
	void debugWriteRatios();
	
	// this function doesn't need to set the segmentationPoints, this can be done in init
	void evaluateCore(const std::vector<SG::candidate_t>& events, unsigned int fine, std::string outputFolder = "");
	
	void evaluateEventsToRatios(const std::vector<SG::candidate_t>& events, unsigned int fine);
	
	void remapAllRatiosNew();

	void extendNormalizationFromDataStep2();

	std::vector<double> computeRemapGridDistanceFractions();
};


} // end of namespace SG

#endif 
