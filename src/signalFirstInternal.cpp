/*
 * class: signalFirstInternal.cpp
 *
 * sgeisselsoeder
 *
 * Date Okt 2013
 *
 * (c) 2013 Antares Collaboration
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>
#include <boost/lexical_cast.hpp>

#include "signalFirstInternal.h"
#include "allScales.h"
#include "candidate.hpp"
#include "clusters.h"
#include "rememberBoostHistogram.h"
#include "utils.hpp"
#include "vtkWriter.h"

#include "sgcu/sgcu.hpp"

namespace SG{


void signalFirstInternal::setOutputFolder(std::string newOutputFolder)
{
	outputFolder_ = newOutputFolder;
}


void addSourceGaussian(std::vector<SG::candidate_t>& artificial, double az, double zen, double sigma, unsigned int num){
	unsigned int counter = 0;
	double sigMax = 3.0*sigma;
	candidate_t center(az, zen);
	while (counter < num)
	{
		double offAz = sgcu::goodRand(2.0*sigMax)-sigMax;
		double offZen = sgcu::goodRand(2.0*sigMax)-sigMax;
		candidate_t temp(az+offAz, zen+offZen);

		double distance = distanceRad(center, temp);
		double sig = distance/sigma;
		double prob = std::exp(-0.5*sig*sig);

		double keepIt = sgcu::goodRand(1.0);
		if (keepIt <= prob)
		{
			artificial.push_back(temp);
			++counter;
		}
	}
}

void addSourceTemplate(std::vector<SG::candidate_t>& artificial, SG::sphericalGrid<double> const & sourceShape, unsigned int num)
{
	unsigned int counter = 0;
	unsigned int counter2 = 0;
	// TODO: check that the shape isn't empty

	// TODO: normalize template to 0:1

	while (counter < num)
	{
		++counter2;
		double offAz = sgcu::goodRand(2.0*M_PI)-M_PI;
		double offZen = sgcu::goodRand(2.0*M_PI)-M_PI;
		candidate_t temp(offAz, offZen);

		// TODO: keep the candidate according to the probability at its closest pixel in the source template
		double prob = sourceShape.findGridpoint(offAz, offZen)->proximity_;
		double keepIt = sgcu::goodRand(1.0);
		if (keepIt <= prob)
		{
			artificial.push_back(temp);
			++counter;
		}
		if (counter2 > 10000*num)
		{
			std::cerr << "DEBUG: addSourceTemplate needs insanely many tries to fill the source template. Something seems to be wrong? Aborting the computation here." << std::endl;
			throw ("party");
		}
	}
}

void addSource(std::vector<candidate_t>& artificial, double az, double zen, double deltaAz, double deltaZen, unsigned int num){
	for (unsigned int i = 0; i < num; ++i){
		double d = sgcu::goodRand(2.0*deltaAz)-deltaAz;
		double d2 = sgcu::goodRand(2.0*deltaZen)-deltaZen;
		artificial.push_back(candidate_t(az+d, zen+d2));
	}
}


void addSource(std::vector<candidate_t>& artificial, double az, double zen, double delta, unsigned int num)
{
	addSource(artificial, az, zen, delta, delta, num);
}

void dumpEvents(std::vector<candidate_t> const & events, std::string outputFile){
	std::cout << "dumpEvents: dumping " << events.size() << " events to file " << outputFile << std::endl;
	std::ofstream raus(outputFile.c_str());
	for (unsigned int i = 0; i < events.size(); ++i)
	{
		raus << events[i] << std::endl;
	}
	raus.close();
}

void dumpDebugOutputOnZeroOrTwenty(unsigned int n, std::vector<double> data, std::string const & filenameWithoutNumberOrEnding)
{
	if ( (n==0) || (n==20) )
	{
//		sgcu::saveVectorToFile(data, std::string(filenameWithoutNumberOrEnding + std::to_string(n) + ".txt") );
		sgcu::saveVectorToFile(data, std::string(filenameWithoutNumberOrEnding + boost::lexical_cast<std::string>(n) + ".txt") );
	}
}

// TODO substitute this?
void signalFirstInternal::writeNormalizationGrid(unsigned int distanceBinToPlot, std::string const & filenameWithoutNumberOrEnding)
{
	for (unsigned int i = 0; i < ratios_.aroundRatioAll.size(); ++i)
		grid_[i].proximity_ = ratios_.aroundRatioAll[i][distanceBinToPlot];
//	myVtkWriter_.write(grid_.gridpoints_, std::string(filenameWithoutNumberOrEnding + std::to_string(distanceBinToPlot) + ".vtk"));
	myVtkWriter_.write(grid_.gridpoints_, std::string(filenameWithoutNumberOrEnding + boost::lexical_cast<std::string>(distanceBinToPlot) + ".vtk"));
}

void signalFirstInternal::writeNormalizationGrids(std::string folderPrefix)
{
	writeNormalizationGrid(0, folderPrefix+"gridNormalizationNew");
	writeNormalizationGrid(20, folderPrefix+"gridNormalizationNew");
	writeNormalizationGrid(100, folderPrefix+"gridNormalizationNew");
}

void computeRadialSymmetricValuesOneScale(allScales& ratios, SG::sphericalGrid<double>& grid, unsigned int scale, std::vector<double>& radialSymmetricValues, double& allSum)
{
	//double allSum = 0.0;
	//std::vector<double> radialSymmetricValues;
	std::vector<unsigned int> counters;
	double currentDec = -99999999999.9;
	// Now store the normalization ratios on a grid
	for (unsigned int i = 0; i < ratios.aroundRatioAll.size(); ++i)
	{
		double value = ratios.aroundRatioAll[i][scale];
			grid[i].proximity_ = value;
			grid[i].multiPurposeValue_ = value;
			allSum += value;

			if (currentDec != grid[i].zenith_){
				currentDec = grid[i].zenith_;
				radialSymmetricValues.push_back(0.0);
				counters.push_back(0);
				//std::cerr << "newly filling in " << radialSymmetricValues.size()-1 << " " << i << std::endl;
			}
			unsigned int index = radialSymmetricValues.size()-1;
			counters[index] += 1;
			radialSymmetricValues[index] += value;
	}

	if (radialSymmetricValues.size() != counters.size()){
		std::cerr << "sizes do not match" << radialSymmetricValues.size() << " " << counters.size() << std::endl;
		throw -3;
	}

	for (unsigned int i = 0; i < radialSymmetricValues.size(); ++i)
	{
		radialSymmetricValues[i] /= static_cast<double>(counters[i]);
	}
}

void smoothRadialSymmetricValues(std::vector<double>& radialSymmetricValues, unsigned int n = 0, std::string outputFolder = "")
{
	if (outputFolder != "") dumpDebugOutputOnZeroOrTwenty(n, radialSymmetricValues, std::string(outputFolder+"radialSymmetricNormalization") );

	// smoothing away the last bins on both ends (don't trust the poles with their small statistic!)
	unsigned int numSmoothAways = 3;
	for (unsigned int i = 0; i < numSmoothAways; ++i)
	{
		radialSymmetricValues[i] = radialSymmetricValues[numSmoothAways];
		radialSymmetricValues[radialSymmetricValues.size()-1-i] = radialSymmetricValues[radialSymmetricValues.size()-1-numSmoothAways];
	}

	if (outputFolder != "") dumpDebugOutputOnZeroOrTwenty(n, radialSymmetricValues, std::string(outputFolder + "radialSymmetricNormalization1Smoothed") );

	sgcu::smoothMedian(radialSymmetricValues);	// this used to be median5 !
	// sgcu::smoothMedian5(radialSymmetricValues);

	if (outputFolder != "")	dumpDebugOutputOnZeroOrTwenty(n, radialSymmetricValues, std::string(outputFolder + "radialSymmetricNormalization2Smoothed") );

	unsigned int numSmoothNoCloseToZeros = 5;	// for ANTARES 2 and IC40 2
//	numSmoothNoCloseToZeros = 3;	// for IceCube, large dataset

	for (unsigned int i = 0; i < numSmoothNoCloseToZeros; ++i)
	{
		smoothNoCloseToZerosKillingEnds(radialSymmetricValues);
	}

	if (outputFolder != "") dumpDebugOutputOnZeroOrTwenty(n, radialSymmetricValues, std::string(outputFolder+"radialSymmetricNormalization3Smoothed") );


	// #CODETODETERMINETHENUMBEROFLOWPASS
	unsigned int numSmoothNoZeros = 60;	// for ANTARES
	numSmoothNoZeros = 15;	// for IceCube, IC40
	numSmoothNoZeros = 2;	// for ANTARES 2 and IC40 2
//	numSmoothNoZeros = 1;	// for IceCube, large dataset
	// std::cout << "debug: smoothRadialSymmetricValues: number of smoothing operations : " << numSmoothNoZeros << std::endl;
	for (unsigned int i = 0; i < numSmoothNoZeros; ++i)
	{
		// sgcu::gaussianSmooth(radialSymmetricValues);	// used to be smoothNoZeros (due to bug without killing ends)
		smoothNoZerosKillingEnds(radialSymmetricValues);
	}

	//sgcu::smooth(radialSymmetricValues);

	if (outputFolder != "") dumpDebugOutputOnZeroOrTwenty(n, radialSymmetricValues, std::string(outputFolder + "radialSymmetricNormalization5Smoothed") );
}

void distributeSymmetricValuesOnSphere(std::vector<double> const & radialSymmetricValues, SG::sphericalGrid<double>& grid)
{
	unsigned int index = 0;
	double currentDec = grid[0].zenith_;
	for (unsigned int i = 0; i < grid.size(); ++i)
	{
		if (grid[i].zenith_ != currentDec)
		{
			currentDec = grid[i].zenith_;
			++index;
		}
		grid[i].proximity_ = radialSymmetricValues[index];
	}
}

void signalFirstInternal::computeOneNormalization(const std::vector<SG::candidate_t>& events, unsigned int fine)
{
	std::string folderBkp = outputFolder_;
	outputFolder_ = "intermediate/step1/";

	timeval startt, endt;
	gettimeofday(&startt, 0);

	// Set up the grid to search and evaluate on
	grid_.setupGrid(fine);
	ratios_.init(stepDegree_);

	// get randomized events
	std::vector<SG::candidate_t> eventsRandomized;
	getRandomEventsFromEvents(events, eventsRandomized);
	myVtkWriter_.write(eventsRandomized, std::string(outputFolder_+"eventsNorm3D") );

	gettimeofday(&endt, 0); std::cout << "Runtime of setup: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	// evaluate on every grid point the distance to all events
	ratios_.computeAllRatios(grid_, eventsRandomized);

	gettimeofday(&endt, 0); std::cout << "Runtime of computeAllRatios: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	unsigned int numDistanceBins = ratios_.aroundRatioAll[0].size();
	for (unsigned int n = 0; n < numDistanceBins; ++n)
	{
		// Check if the ratios and the grid match in size
		if (ratios_.aroundRatioAll.size() != grid_.size()){
			std::cerr << "Size mismatch ratios: " << ratios_.aroundRatioAll.size() << " vs grid:" << grid_.size() << std::endl;
			return;
		}

		double allSum = 0.0;
		std::vector<double> radialSymmetricValues;
		computeRadialSymmetricValuesOneScale(ratios_, grid_, n, radialSymmetricValues, allSum);

		smoothRadialSymmetricValues(radialSymmetricValues, n, outputFolder_);

		distributeSymmetricValuesOnSphere(radialSymmetricValues, grid_);

		renormalizeGrid(grid_, allSum);

		// store the result for this scale in the full set of ratios
		for (unsigned int i = 0; i < grid_.size(); ++i)
		{
			ratios_.aroundRatioAll[i][n] = grid_[i].proximity_;
		}
	}

	gettimeofday(&endt, 0); std::cout << "Runtime of the computing expectations from ratios: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	outputFolder_ = folderBkp;
}



void signalFirstInternal::computeNormalizationsFromDataStep1(const std::vector<SG::candidate_t>& events, unsigned int fine)
{
	std::cout << "Normalization Step1 with " << events.size() << " events" << std::endl;

	std::string outputPath = "intermediate/step1/";
	setOutputFolder(outputPath);
	computeOneNormalization(events, fine);

//	writeNormalizationGrids(folderPrefix+"gridNormalizationNew");
	ratios_.dumpRatiosAsNormalization(outputPath+"normalizationDataStep1", events.size());
}


void signalFirstInternal::init(unsigned int fine)
{
	grid_.setupGrid(fine);
	segmentationPoints_.clear();

	for (unsigned int i = 1; i < 5; ++i)
	{
		segmentationPoints_.push_back( std::pow(0.5, i) );
	}
}







void signalFirstInternal::performEvaluate(std::vector<SG::candidate_t>& events, unsigned int fine) {
	std::string outputBkp = outputFolder_;
	outputFolder_ = "intermediate/step5/";

	std::cout << "Evaluating " << events.size() << " events" << std::endl;

	myVtkWriter_.write(events, outputFolder_+"events3D");
	writeEquatorialProjection(events, outputFolder_+"hammerProjEventsEquatorial.txt");
	writeGalacticProjection(events, outputFolder_+"hammerProjEventsGalactic.txt");

	evaluateCore(events, fine, "");
	
	printMaxRelevanceSegmentationwise(finalClustersMulti_, grid_, segmentationPoints_);

	std::string folderPrefixStep2 = "intermediate/step2/";
	renormalizeStep2MaximaOnly(finalClustersMulti_, segmentationPoints_, folderPrefixStep2+"normalizationDataStep2", "clusters_segment");

 	std::cout << "Maximal relevance, all segmentations, all clusters, excl. trial factor: " << getMaxRelevance(finalClustersMulti_) << std::endl;
 	printMaxRelevanceSegmentationwise(finalClustersMulti_, grid_, segmentationPoints_);

//	for (unsigned int i = 0; i < finalClustersMulti_.size(); ++i)
//	{
//		drawClusters(finalClustersMulti_[i], grid_);
//		myVtkWriter_.write(grid_.gridpoints_, outputFolder_+"grid3DIntermed" + boost::lexical_cast<std::string>(++globalOutputCounter));
//	}
	// output the stencil of the large cluster
//	drawClusters(finalClustersMulti_[0], grid_);
	//grid_.save("stencilGrid.grid");

//	for (unsigned int i = 0; i < finalClustersMulti_.size(); ++i){
//		std::cerr << "DEBUG drawing debug output" << i << std::endl;
//		drawClusters(finalClustersMulti_[i], grid_);
//		std::cerr << "DEBUG writing debug output" << i << std::endl;
//		myVtkWriter_.write(grid_.gridpoints_, outputFolder_+"grid3DIntermed" + boost::lexical_cast<std::string>(++globalOutputCounter));
//	}
	
//	std::cerr << "DEBUG renorm3" << std::endl;
	renormalizeClusterRelevances(finalClustersMulti_);

	//std::cerr << "DEBUG finished renorm3" << std::endl;
 	std::cout << "Maximal relevance, all segmentations, all clusters, incl. trial factor, with metrics mask: " << getMaxRelevance(finalClustersMulti_) << std::endl;
 	printMaxRelevanceSegmentationwise(finalClustersMulti_, grid_, segmentationPoints_);
	//std::cerr << "printing clusters in ASCII files" << std::endl;
	printClustersToASCII(finalClustersMulti_, grid_);
	
	printEventInfoPerCluster(events, grid_, finalClustersMulti_, outputFolder_);

	std::cerr << "Drawing finalized clusters" << std::endl;
	for (unsigned int i = 0; i < finalClustersMulti_.size(); ++i)
	{
//		std::cerr << "Drawing finalized cluster " << i << " / " << finalClustersMulti_.size() << std::endl;
		drawClusters(finalClustersMulti_[i], grid_);
//		myVtkWriter_.write(grid_.gridpoints_, outputFolder_+"grid3DIntermed" + boost::lexical_cast<std::string>(++globalOutputCounter));
////	grid_.writeAsGalacticHammerProj(outputFolder_);
		grid_.writeAsEquatorialHammerProj(outputFolder_);
	}

	outputFolder_ = outputBkp;
}

void setFirstScaleToZero(allScales & scales)
{
	for (unsigned int i = 0; i < scales.aroundRatioAll.size(); ++i)
	{
		scales.aroundRatioAll[i][0] = 0.0;
	}
}

void evaluateRatiosToGrid(SG::allScales const & scales, SG::sphericalGrid<double> & grid, std::string outputFolder)
{
	//TODO: check if grid is actually set up at this stage
	sumScalesToGrid(scales, grid);
	grid.bkpProximity();

	std::cout << "OF: " << outputFolder << std::endl;
//	myVtkWriter_.write(grid.gridpoints_, outputFolder+"grid3DIntermed" + boost::lexical_cast<std::string>(++globalOutputCounter));
	grid.writeAsEquatorialHammerProj(outputFolder);
//	writeGridAsGalacticHammerProj(grid);

	// dump grid without hammerProjection
	save(grid, outputFolder+"/gridComputationCoordsDetailed.grid");
}


void checkConsistentSize(std::vector<std::vector<std::vector<SG::cluster> > > const & allEntries, std::vector<double> const & segmentationPoints)
{
	unsigned int numberOfPseudoExperimentsPerSegmentation = allEntries[0].size();
	for (unsigned int n = 1; n < allEntries.size(); ++n)
	{
		if (allEntries[n].size() != numberOfPseudoExperimentsPerSegmentation)
		{
			std::cerr << "Not all Segmentation have the same number of pseudoexperiments!" <<
					"I expected " << numberOfPseudoExperimentsPerSegmentation << " but for segmentation " << n  <<
					" with threshold " << segmentationPoints[n] << " I got " << allEntries[n].size() <<
					"! Cannot proceed to compute normalization3 from normalization 2 values." << std::endl;
			throw -1;
		}
	}
	if (numberOfPseudoExperimentsPerSegmentation == 0)
	{
		std::cerr << "Failed to load any clusters from the pseudo experiments." << std::endl;
		throw -1;
	}
}


// renormalize all clusters of this one pseudo experiment with all their metrics, overwrite the raw metrics values by the renormalized results
// and return the highest renormalized relevance of any cluster and any metric of this pseudo experiment
double renormalizeOnePseudoExperimentAllClustersAllMetrics(std::vector<std::vector<struct rememberBoostHistogram > > const & histsAllSegmentationsAllRelevances,
		std::vector<double> const & segmentationPoints, unsigned int pseudoexperimentIndex,
		std::vector<std::vector<std::vector<SG::cluster> > > & allClustersFromPseudoExperiments)
{
	std::vector<double> relevances;
	// check all segmentations for one pseudo experiment
	for (unsigned int n = 0; n < segmentationPoints.size(); ++n)
	{
		// std::cerr << "in seg " << n << " / " << segmentationPoints_.size() << std::endl;

		std::vector<SG::cluster> comparisonClusters;
		// Now consider all clusters that have been produced in one pseudo experiment with one segmentation threshold
		comparisonClusters = allClustersFromPseudoExperiments[n][pseudoexperimentIndex];
		if (comparisonClusters.size() == 0)
		{
			std::cerr << "No cluster loaded from segmentation " << n << " at " << segmentationPoints[n] << std::endl;
			continue;
		}

		// now translate the relevance values (size, highest value, etc.) to
		// a p-value (how likely is it to observe such a value or greater) and therefore
		// to a significance (sometimes also called renormalized relevance or pre-trial significance)
		renormalizeCore(comparisonClusters, histsAllSegmentationsAllRelevances[n]);	// renormalize all metrics of all clusters of this pseudo experiment of this segmentation

		// use the renormalized clusters from now on
		allClustersFromPseudoExperiments[n][pseudoexperimentIndex] = comparisonClusters;

		for (unsigned int j = 0; j < comparisonClusters.size(); ++j)
		{	// just to make the relevances vector shorter, inprinciple the maximum is not needed here
			relevances.push_back(*std::max_element(comparisonClusters[j].relevances.begin(), comparisonClusters[j].relevances.end()));
		}
	}
	// store the maximal renormalized relevance value of any cluster and any metric for this pseudo experiment
	return ( *max_element(relevances.begin(), relevances.end()) );
}


//void signalFirstInternal::extendNormalizationFromStep2ToStep3()
void signalFirstInternal::renormalizeAllPseudoexperimentsAllClustersAllMetrics()
{
	std::string inputFolder = "intermediate/step2/normalizationDataStep2";
	std::string searchString1 = "clusters_segment";
	
	std::cout << "Computing normalization3 distribution from normalization2 pseudo-experiments." << std::endl;
	std::cout << "Loading comparison values ..." << std::endl;

	// Load all comparison events
	std::vector<std::vector<struct rememberBoostHistogram > > histsAllSegmentationsAllRelevances;
	std::vector<std::vector<std::vector<SG::cluster> > > allEntries;
	getComparisonValuesAndHistograms(inputFolder, searchString1, segmentationPoints_, allEntries, histsAllSegmentationsAllRelevances);
	checkConsistentSize(allEntries, segmentationPoints_);

	std::vector<double> maximalValues;
	for (unsigned int i = 0; i < allEntries[0].size(); ++i){	// for all pseudo experiments of this segmentation
		// std::cerr << "Processing pseudo-experiment " << i << " / " << allEntries[0].size() << std::endl;
		maximalValues.push_back( renormalizeOnePseudoExperimentAllClustersAllMetrics(histsAllSegmentationsAllRelevances, segmentationPoints_, i, allEntries) );
	}
	
	const std::string normalizationFolderOutput = "intermediate/step4/normalizationDataStep3/";
	sgcu::saveVectorToFile(maximalValues, normalizationFolderOutput + "normalizationStep3_Maxima.txt");
}



void signalFirstInternal::evaluateEventsToRatios(const std::vector<SG::candidate_t>& events, unsigned int fine) {
	std::cout << "Evaluate events to ratios" << std::endl;

	timeval startt, endt;
	gettimeofday(&startt, 0);

	// Set up the grid to search and evaluate on
	grid_.setupGrid(fine);
	ratios_.init(stepDegree_);

	gettimeofday(&endt, 0); std::cout << "runtime of init: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	// The evaluateDensity sets the ratios, not the grid
	evaluateDensity(events);
	ratios_.assertNoNansAndInfs();

	gettimeofday(&endt, 0); std::cout << "runtime of evaluateDensity: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	debugWriteRatios();

	segmentAllRatiosScalewise(ratios_);
	ratios_.assertNoNansAndInfs();

	debugWriteRatios();

	gettimeofday(&endt, 0); std::cout << "runtime of segmentAllRatiosScalewise: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	remapAllRatiosNew();
	ratios_.assertNoNansAndInfs();

	gettimeofday(&endt, 0); std::cout << "runtime of remapAllRatiosNew: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);
}



std::vector<cluster> evaluateGridToClusters(SG::sphericalGrid<double> & grid, double keepFraction)
{
	std::cout << "grid non zero before keepTop  " << grid.countNumberOfNonZeroValues() << std::endl;
	grid.keepTop(keepFraction);
	std::cout << "grid non zero after keepTop  " << grid.countNumberOfNonZeroValues() << std::endl;
//	if (doDebugOutput) myVtkWriter_.write(grid.gridpoints_, outputFolder_+"grid3DIntermed" + boost::lexical_cast<std::string>(++globalOutputCounter)); //21

	grid.medianFilter();
//	if (doDebugOutput) myVtkWriter_.write(grid.gridpoints_, outputFolder_+"grid3DIntermed" + boost::lexical_cast<std::string>(++globalOutputCounter)); //21

	return clustersFromGrid(grid);
}




void evaluateGridToSegmentsToClusters(SG::sphericalGrid<double> const & grid, std::vector<double> const & segmentationPoints,
		std::string outputFolder, std::vector<std::vector<cluster> > & finalClustersMulti)
{
	std::cout << "Evaluate grid to segments to clusters" << std::endl;
	grid.assertNoNansAndInfs();

	finalClustersMulti.clear();
	SG::sphericalGrid<double> gridToWork = grid;

	unsigned int thisTimeIndex = time(NULL);
	std::cout << "segmentationPoints_.size() " << segmentationPoints.size() << std::endl;
	for (unsigned int i = 0; i < segmentationPoints.size(); ++i)
	{
		std::cout << i+1 << "/" << segmentationPoints.size() << "  keepFraction=" << segmentationPoints[i] << std::endl;

		gridToWork = grid;
		std::vector<cluster> finalClusters = evaluateGridToClusters(gridToWork, segmentationPoints[i]);

		if (outputFolder != "")
		{
			writeClusters(outputFolder + "/clusters_segment" + boost::lexical_cast<std::string>(segmentationPoints[i])
					+ "_" + boost::lexical_cast<std::string>(thisTimeIndex) + ".txt", finalClusters);
		}
		
		finalClustersMulti.push_back(finalClusters);
	}
}


void signalFirstInternal::evaluateEventsToSkymap(const std::vector<SG::candidate_t>& events, unsigned int fine){
//	myVtkWriter_.write(events, outputFolder_+"events3D");
//	writeGalacticProjection(events, outputFolder_+"hammerProjEventsGalactic.txt");
	writeEquatorialProjection(events, outputFolder_+"hammerProjEventsEquatorial.txt");
	dumpEvents(events, outputFolder_+"events.txt");

	timeval startt, endt;
	gettimeofday(&startt, 0);

	evaluateEventsToRatios(events, fine);

	gettimeofday(&endt, 0); std::cout << "runtime of evalEventsToRatios: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);

	evaluateRatiosToGrid(ratios_, grid_, outputFolder_);

	writeEquatorialProjection(grid_.gridpoints_, outputFolder_+"hammerProjGridDebugEquatorial.txt");

	gettimeofday(&endt, 0); std::cout << "runtime of evaluateRatiosToGrid: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);
}


void signalFirstInternal::evaluateCore(const std::vector<SG::candidate_t>& events, unsigned int fine, std::string outputFolder)
{
	evaluateEventsToSkymap(events, fine);
	evaluateGridToSegmentsToClusters(grid_, segmentationPoints_, outputFolder, finalClustersMulti_);
}


void signalFirstInternal::evaluateDensity(const std::vector<SG::candidate_t>& events)
{
	setToZero(ratios_);
	unsigned int numBins = static_cast<unsigned int>(maximumDistanceToEvaluateInDegree/stepDegree_)+1;
	if (normalization_.size() == 0)
	{
		std::cout << "Automatically trying to load normalization from step1." << std::endl;
		int ret = normalization_.loadNormalization2(grid_.size(), numBins, events.size(), "intermediate/step1/normalizationDataStep1");
		if (ret != 0)
		{
			std::cerr << "Failed to load normalization for step1 from intermediate/step1/normalizationDataStep1 with error code " << ret << std::endl;
			throw -1;
		}
	}
	if (normalization_.size() == 0)
	{
		std::cerr << "evaluateDensity is only possible if a normalization has been computed and loaded before!" << std::endl;
		throw -2;
	}
	
	ratios_.computeAllRatios(grid_, events);

	std::cout << countNumberOfNonZeros(ratios_) << " non-zero ratio values" << std::endl;

//	debugWriteRatios();

	ratios_.computePoissonProbabilityForAllRatios(normalization_, events.size());

	invertRatiosAndTakeLog10(ratios_);

	washOutWithinScales(ratios_, grid_);
	washOutBetweenScales(ratios_);
}



void addRandomEventsFromEvents(const std::vector<SG::candidate_t>& events, std::vector<SG::candidate_t>& eventsRandomized)
{
	// use events to generate a lot of pseudo random events here!
	for (unsigned int i = 0; i < events.size(); ++i)
	{
		candidate_t c = events[i];
		c.azimuth_ = sgcu::goodRand(2.0*M_PI);
		eventsRandomized.push_back(c);
	}
}

void getRandomEventsFromEvents(const std::vector<SG::candidate_t>& events, std::vector<SG::candidate_t>& eventsRandomized)
{
	eventsRandomized.clear();
	addRandomEventsFromEvents(events, eventsRandomized);
}




void signalFirstInternal::computeNormalizationsFromDataStep2(const std::vector<SG::candidate_t>& eventsRaw, unsigned int fine, unsigned int numberOfRepetitions){
	std::cout << "Normalization Step2 with " << eventsRaw.size() << " events and " << numberOfRepetitions << " repetitions" << std::endl;

	std::string folderPrefix = "intermediate/step2/";
	setOutputFolder(folderPrefix);
	std::string normalizationFolder = folderPrefix+"normalizationDataStep2";

	for (unsigned int n = 0; n < numberOfRepetitions; ++n){
		std::cout << "Normalizations " << n+1 << "/" << numberOfRepetitions << std::endl;
	
		std::vector<SG::candidate_t> eventsRandomized;
		getRandomEventsFromEvents(eventsRaw, eventsRandomized);	// only give a new random azimuth

		evaluateCore(eventsRandomized, fine, normalizationFolder);
	}
}

//// This outputs a series of grids , each with the ratios found for one distance bin
//void signalFirstInternal::dumpRatiosForDebug(){
//	std::vector<std::vector<float> > temp;
//	temp.resize(ratios_.aroundRatioAll.size());
//	unsigned int numFirst = 6;
//	unsigned int numSecond = 14;
//	for (unsigned int i = 0; i < temp.size(); ++i){
//		temp[i].resize(numFirst+numSecond+1);
//	}
//	for (unsigned int i = 0; i < temp.size(); ++i){
//		for (unsigned int j = 0; j < numFirst; ++j){
//			temp[i][j] = ratios_.aroundRatioAll[i][j];
//		}
//	}
//	for (unsigned int i = 0; i < temp.size(); ++i){
//		for (unsigned int j = numFirst; j < temp[i].size(); ++j){
//			temp[i][j] = ratios_.aroundRatioAll[i][(j-numFirst)*10];
//		}
//	}
//	for (unsigned int i = 0; i < temp.size(); ++i){
//		unsigned int j = numFirst+numSecond;
//		temp[i][j] = ratios_.aroundRatioAll[i][ratios_.aroundRatioAll[i].size()-1];
//	}
//	myVtkWriter_.write(grid_.gridpoints_, temp, "grid3D");
//}

void washOutWithinScales(allScales & scales, SG::sphericalGrid<double> const & grid)
{
	const unsigned int numberOfGeneralSmearings = 1;
	const unsigned int numSmearCentered = 1;
	
	SG::sphericalGrid<double> tempGrid = grid;
	for (unsigned int k = 0; k < numberOfGeneralSmearings; ++k)
	{
		// apply spatial filters here for every distance bin
		for (unsigned int i = 0; i < scales.aroundRatioAll[0].size(); ++i)
		{
			
			// copy ratios to the grid	// TODO check if this can't use a function that can be reused
			for (unsigned int j = 0; j < tempGrid.size(); ++j)
			{
				tempGrid[j].proximity_ = scales.aroundRatioAll[j][i];
			}
			
			// perform some filters
			for (unsigned int rep = 0; rep < numSmearCentered; ++rep)
				tempGrid.smearCentered();
			
			// copy values back from the grid to the ratios
			for (unsigned int j = 0; j < tempGrid.size(); ++j)
			{
				scales.aroundRatioAll[j][i] = tempGrid[j].proximity_;
			}
		}
	}
}

void washOutBetweenScales(allScales & scales){	// TODO: exchange tempRatios and scales (saves a copy of the whole thing)
	const unsigned int numDistanceGaussians = 1;

	allScales tempRatios = scales;
	// smooth the results over neighbouring distance bins
	unsigned int numBins = scales.aroundRatioAll[0].size();
	for (unsigned int rep = 0; rep < numDistanceGaussians; ++rep){
		for (unsigned int i = 0; i < scales.aroundRatioAll.size(); ++i){
			// narrow centered smear
			for (unsigned int j = 1; j < numBins-1; ++j){
				tempRatios.aroundRatioAll[i][j] = (scales.aroundRatioAll[i][j-1] + 4.0*scales.aroundRatioAll[i][j] + scales.aroundRatioAll[i][j+1]) / 6.0;
			}
			// approximately correct the first and last
			tempRatios.aroundRatioAll[i][0] = (2.0*scales.aroundRatioAll[i][0] + scales.aroundRatioAll[i][1]) / 3.0;
			tempRatios.aroundRatioAll[i][numBins-1] = (scales.aroundRatioAll[i][numBins-2] + 2.0*scales.aroundRatioAll[i][numBins-1])/3.0;
		}
	}
	scales = tempRatios;
}

std::vector<double> signalFirstInternal::computeRemapGridDistanceFractions()
{
	// std::cout << "Computing the remapping grid numbers how many points are at a certain distance" << std::endl;
	// This is the mapping for every distance how many other gridpoints are at this distance
	std::vector<double> mapping;
	unsigned int numBins = 180;
	mapping.resize(numBins);
	
	// we compute a mapping for every distance
	for (unsigned int i = 0; i < numBins; ++i)
	{
		// how many gridpoints should be at that distance?
		double tempAngle = static_cast<double>( M_PI*i / (18*grid_.fine()) );
		double densityFactor = sin(tempAngle);
		unsigned int howManyGridpoints = static_cast<unsigned int>(grid_.fine()*36*densityFactor + 0.5);
		if (howManyGridpoints == 0) howManyGridpoints = 1;
		mapping[i] = 1.0/static_cast<double>(howManyGridpoints);
	}
	return mapping;
}


void signalFirstInternal::remapAllRatiosNew(){
//	cout << "why aren't you using remapAllRatiosNewAgain() ???" << std::endl;

	unsigned int numGridpoints = grid_.size();
	if (numGridpoints == 0){
		std::cerr << "Remapping ratios is only possible if there is a grid!" << std::endl;
		return;
	}
	unsigned int numBins = ratios_.aroundRatioAll[0].size();
	std::cout << "Remapping ratios" << std::endl;

	double stepRad = 0.5/180*M_PI;
	allScales tempRatios = ratios_;
	for (unsigned int i = 0; i < numGridpoints; ++i){
		for (unsigned int j = 0; j < numBins; ++j){
			tempRatios.aroundRatioAll[i][j] = 0.0;
		}
	}

	static std::vector<double> distanceFractions = computeRemapGridDistanceFractions();

	// remap proximities on every gridpoint
//	omp_set_num_threads(16);
	#pragma omp parallel for
	for (int i = 0; i < (int)numGridpoints; ++i){
		// check with every other gridpoint if it has exactly the distance the found cluster
		for (unsigned int j = 0; j < numGridpoints; ++j){
			double dist = distanceRad(grid_[i], grid_[j]);

			// what is this distance in gridpoint bins?
			unsigned int binIndex = static_cast<unsigned int>(dist/stepRad+0.5);
			if (binIndex < numBins){
				double proximity = ratios_.aroundRatioAll[i][binIndex];
				// do not consider gridpoints which have NO proximity to redistribute!
				if (proximity <= 0.0){
					continue;
				}
				// else: if they have a proximity (meaning e.g. they are part of a found cluster) then process this gridpoint

				// always set the minimal multipurpose value that helped create the resulting remapped gridpoint
//				double multiPurposeValue = binIndex;

				double proximitySplit = proximity*distanceFractions[binIndex];

				#pragma omp critical
				{
					// then add some of the proximity of this cluster point to the found point in the right distance
					tempRatios.aroundRatioAll[j][binIndex] += proximitySplit;
				}
			}
		}
	}
	// finally store our work
	ratios_ = tempRatios;
}


// TODO speed up the remapping, but nothing really worked so far
//void signalFirstInternal::remapAllRatiosNewAgain(){
//	unsigned int numGridpoints = grid_.size();
//	if (numGridpoints == 0){
//		std::cerr << "Remapping ratios is only possible if there is a grid!" << std::endl;
//		return;
//	}
//	unsigned int numBins = ratios_.aroundRatioAll[0].size();
//	std::cout << "Remapping ratios" << std::endl;
//
//	timeval startt, endt;
//	gettimeofday(&startt, 0);
//
//	double stepRad = ratios_.stepSize_;
//	std::cout << "debug: stepRad = " << stepRad << " , used to be " << 0.5/180.0*M_PI << std::endl;
//	allScales tempRatios = ratios_;
//	for (unsigned int i = 0; i < numGridpoints; ++i){
//		for (unsigned int j = 0; j < numBins; ++j){
//			tempRatios.aroundRatioAll[i][j] = 0.0;
//		}
//	}
//
//	gettimeofday(&endt, 0); std::cout << "runtime of zeroing: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);
//
//	static std::vector<double> distanceFractions = computeRemapGridDistanceFractions();
//
//	gettimeofday(&endt, 0); std::cout << "runtime of computeRemapGridDistanceFractions: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);
//
//	// remap proximities on every gridpoint
////	omp_set_num_threads(16);
//	#pragma omp parallel for
//	for (int i = 0; i < 1000; ++i){
////	for (int i = 0; i < (int)numGridpoints; ++i){
//		// check with every other gridpoint if it has exactly the distance the found cluster
//		for (unsigned int j = 0; j < numGridpoints; ++j){
//			double dist = distanceRad(grid_[i], grid_[j]);
//
//			// what is this distance in gridpoint bins?
//			unsigned int binIndex = static_cast<unsigned int>(dist/stepRad+0.5);
//			if (binIndex < numBins){
//				double proximity = ratios_.aroundRatioAll[i][binIndex];
//				// do not consider gridpoints which have NO proximity to redistribute!
//				if (proximity <= 0.0){
//					continue;
//				}
//				// else: if they have a proximity (meaning e.g. they are part of a found cluster) then process this gridpoint
//
//				// always set the minimal multipurpose value that helped create the resulting remapped gridpoint
////				double multiPurposeValue = binIndex;
//
//				double proximitySplit = proximity*distanceFractions[binIndex];
//
//				#pragma omp critical
//				{
//					// then add some of the proximity of this cluster point to the found point in the right distance
//					tempRatios.aroundRatioAll[j][binIndex] += proximitySplit;
//				}
//			}
//		}
//	}
//
//	gettimeofday(&endt, 0); std::cout << "runtime of remapping core: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);
//
//	// finally store our work
//	ratios_ = tempRatios;
//
//	gettimeofday(&endt, 0); std::cout << "runtime of copying ratios: " << endt.tv_sec - startt.tv_sec << std::endl; gettimeofday(&startt, 0);
//}

// TODO check
void sumScalesToGrid(allScales const & scales, SG::sphericalGrid<double> & grid)
{
	for (unsigned int i = 0; i < scales.aroundRatioAll.size(); ++i)
	{
		grid[i].proximity_ = sgcu::accumulate(scales.aroundRatioAll[i]);
	}
}

void invertRatiosAndTakeLog10(allScales & scales)
{
	for (unsigned int i = 0; i < scales.aroundRatioAll.size(); ++i)
	{
		for (unsigned int j = 0; j < scales.aroundRatioAll[i].size(); ++j)
		{
			scales.aroundRatioAll[i][j] = log10( 1.0/scales.aroundRatioAll[i][j] );
		}
	}
}


void writeRatiosToGrid(allScales const & scales, unsigned int distanceIndex, SG::sphericalGrid<double> & grid)
{
	static unsigned int outputCounter = 0;
	++outputCounter;
	for (unsigned int i = 0; i < grid.size(); ++i)
	{
		grid[i].proximity_ = scales.aroundRatioAll[i][distanceIndex];
	}
	writeEquatorialProjection(grid.gridpoints_, "hammerProjEquatorial_distance" +
			boost::lexical_cast<std::string>(distanceIndex) + "_nr"+boost::lexical_cast<std::string>(outputCounter)+".txt");
//	myVtkWriter_.write(grid_.gridpoints_, outputFolder_ + "grid3DIntermed" + boost::lexical_cast<std::string>(outputCounter));
}


void signalFirstInternal::debugWriteRatios()
{
	//std::cerr << " We don't write debugRatios here anymore!" << std::endl;
	return;

	writeRatiosToGrid(ratios_, 0, grid_);
	writeRatiosToGrid(ratios_, 6, grid_);
	writeRatiosToGrid(ratios_, 20, grid_);
	writeRatiosToGrid(ratios_, 100, grid_);
}


//void signalFirstInternal::discardClustersOutsideOfTemplate(double minFraction){
//	//std::cerr << "DEBUG Starting discard cluster based on stencil" << std::endl;
//	// if there is a stencil: discard every cluster that is not in the stencil (at least to a minimum fraction)
//	if (stencilGrid_.gridpoints_.size() != 0){
//		std::cerr << "Beware! Stencil found! Cutting away everything that is not in this stencil!" << std::endl;
//		//for (unsigned int i = 0; i < finalClustersMulti_.size(); ++i){
//		// we operate on finalClusters_ now, no outer loop neccessary
//			//for (unsigned int j = 0; j < finalClustersMulti_[i].size(); ++j){
//			for (unsigned int j = 0; j < finalClusters_.size(); ++j){
//				// check this one cluster
//				//std::vector<unsigned int> clusterpoints = finalClustersMulti_[i][j].gridpoints;
//				std::vector<unsigned int> clusterpoints = finalClusters_[j].gridpoints;
//				//std::cerr << "Testing cluster no " << j << " of size " << clusterpoints.size() << std::endl;
//
//				unsigned int insideCounter = 0;
//				// count how many pixels are inside
//				for (unsigned int k = 0; k < clusterpoints.size(); ++k){
//					unsigned int index = clusterpoints[k];
//					if (stencilGrid_.gridpoints_[index].proximity_ > 0.0){
//						++insideCounter;
//					}
//				}
//				// compute the overlaps-with-stencil-fraction:
//				double fraction = static_cast<double>(insideCounter)/static_cast<double>(clusterpoints.size());
//				//std::cerr << "Cluster has a fraction of " << fraction << std::endl;
//				// discard all gridpoints for this cluster if it doesn't match well enough
//				if (fraction < minFraction){
//					//std::cerr << "Discarding cluster! ###########################" << std::endl;
//					//finalClustersMulti_[i][j].gridpoints.clear();
//					//finalClusters_[j].gridpoints.clear();
//					//for (unsigned int k = 0; k < finalClustersMulti_[i][j].relevance.size(); ++k){
//					//std::cerr << "relevances before: " ;
//					for (unsigned int k = 0; k < finalClusters_[j].relevances.size(); ++k){
//						//std::cerr << finalClusters_[j].relevance[k] << " " << std::endl;
//						//finalClustersMulti_[i][j].relevance[k] = 0.0;
//						finalClusters_[j].relevances[k] = 0.0;
//					}
//					//std::cerr << std::endl;
//				}
//			}
//		//}
//	}
//	//std::cerr << "DEBUG finished discarding cluster" << std::endl;
//}

//void signalFirstInternal::readTemplate(unsigned int fine, std::string fileName, double significanceThreshold) {
//	stencilGrid_.setupGrid(fine);	// shouldn't be neccessary, but I think it wouldn't work without
//	stencilGrid_.load(fileName);
//	//stencilGrid_.correctEquatorialRightAscensionRange();
//
//	for (unsigned int i = 0; i < stencilGrid_.gridpoints_.size(); ++i){
//		if (stencilGrid_.gridpoints_[i].proximity_ < significanceThreshold){
//			stencilGrid_.gridpoints_[i].proximity_ = 0.0;
//		}
//	}
//	stencilGrid_.writeHammerProjection("hammerProjStencilEquatorial.txt");
//}






// put points closely together on the coordinates grid lines
void writeGridLines(){
	std::cerr << "writing grid Lines" << std::endl;
	const double oneDeg = M_PI/180.0;
	const double delta = 0.5*oneDeg;
	const int decStart = -90;
	const int decEnd = 90;
	const int raStart = -180;
	const int raEnd = 180;

	std::vector<SG::candidate<double> > gridLines;
	SG::candidate<double> tempPoint;
	for (int i = decStart; i <= decEnd; i += 10){	// declination from -90 to 90
		double dec = i*oneDeg;
		tempPoint.zenith_ = dec;
		double ra = raStart*oneDeg;
		while (ra <= raEnd*oneDeg){
			tempPoint.azimuth_ = ra;
			tempPoint.proximity_ = ra;
			gridLines.push_back(tempPoint);
			ra += delta;
		}
	}
	std::cerr << "dec finished" << std::endl;
	for (int i = raStart; i <= raEnd; i += 10){ // right ascension from 0 to 360
		double ra = i*oneDeg;
		tempPoint.azimuth_ = ra;
		tempPoint.proximity_ = ra;
		double dec = decStart*oneDeg;
		while (dec <= decEnd*oneDeg){
			tempPoint.zenith_ = dec;
			gridLines.push_back(tempPoint);
			dec += delta;
		}

	}
	//findRange(gridLines);
	std::cerr << "ra finished" << std::endl;
	writeHammerProjection(gridLines, "gridLines.grid");
	std::cerr << "file written" << std::endl;
}

void writeGridLabelsGalactic(){
	std::cerr << "writing grid Labels" << std::endl;
	const double oneDeg = M_PI/180.0;
//	const double delta = 1.0*oneDeg;
	const int decStart = -90;
	const int decEnd = 90;
	const int raStart = -180;
	const int raEnd = 180;
	const int stepSize = 30;

	std::vector<SG::candidate<double> > gridLines;
	SG::candidate<double> tempPoint;
	for (int i = decStart; i <= decEnd; i += 10){	// declination from -90 to 90
		double dec = i*oneDeg;
		double ra = raStart*oneDeg;
		tempPoint.zenith_ = dec;
		tempPoint.azimuth_ = ra;
		tempPoint.proximity_ = i;
		gridLines.push_back(tempPoint);
	}
	std::cerr << "dec finished" << std::endl;
	for (int i = raStart+stepSize; i <= raEnd-stepSize; i += stepSize){ // right ascension from 0 to 360
		double ra = i*oneDeg;
		double dec = 0.5*(decStart+decEnd)*oneDeg;
		tempPoint.azimuth_ = ra;
		tempPoint.zenith_ = dec;
		tempPoint.proximity_ = i;
		gridLines.push_back(tempPoint);
	}
	//findRange(gridLines);
	std::cerr << "ra finished" << std::endl;
//	writeHammerProjection("gridLabelsGal.grid", gridLines);
	writeHammerProjection(gridLines, "gridLabelsGal.grid");
	std::cerr << "file written" << std::endl;
}

void writeGridLabelsEquatorial(){
	std::cerr << "writing grid Labels" << std::endl;
	const double oneDeg = M_PI/180.0;
//	const double delta = 1.0*oneDeg;
	const int decStart = -90;
	const int decEnd = 90;
	const int raStart = -180;
	const int raEnd = 180;
	const int stepSize = 30;

	std::vector<SG::candidate<double> > gridLines;
	SG::candidate<double> tempPoint;
	for (int i = decStart; i <= decEnd; i += 10){	// declination from -90 to 90
		double dec = i*oneDeg;
		double ra = raStart*oneDeg;
		tempPoint.zenith_ = dec;
		tempPoint.azimuth_ = ra;
		tempPoint.proximity_ = i;
		gridLines.push_back(tempPoint);
	}
	std::cerr << "dec finished" << std::endl;
	for (int i = raStart+stepSize; i <= raEnd-stepSize; i += stepSize){ // right ascension from 0 to 360
		double ra = i*oneDeg;
		double dec = 0.5*(decStart+decEnd)*oneDeg;
		tempPoint.azimuth_ = ra;
		tempPoint.proximity_ = i+180;
		tempPoint.zenith_ = dec;
		gridLines.push_back(tempPoint);
	}
	//findRange(gridLines);
	std::cerr << "ra finished" << std::endl;
//	writeHammerProjection("gridLabelsEquat.grid", gridLines);
	writeHammerProjection(gridLines, "gridLabelsEquat.grid");
	std::cerr << "file written" << std::endl;
}





} // end of namespace SG
