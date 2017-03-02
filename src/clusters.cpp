#include <vector>
#include <iostream>
#include <cmath>	// for M_PI
#include <boost/lexical_cast.hpp>
//#include "candidate.hpp"
#include "clusters.h"
#include "searchpoint.h"
#include <dirent.h>	// for file listing
#include <algorithm>    // std::sort
#include <boost/math/special_functions/erf.hpp>

using namespace std;

namespace SG
{


void fuseAlreadyOverlappingClusters(std::vector<std::vector<std::pair<unsigned int, unsigned int> > >& clusters)
{
	// go through every cluster
	for (unsigned int i = 0; i < clusters.size(); ++i){
		// go through every gridpoint of this cluster
		for (unsigned int j = 0; j < clusters[i].size(); ++j){
			// search for this gridpoint in every other cluster
			for (unsigned int k = 0; k < clusters.size(); ++k){
				if (k != i){
					bool theyHaveOverlap = false;
					// check if the other cluster also contains this gridpoint:
					if (find(clusters[k].begin(), clusters[k].end(), clusters[i][j]) != clusters[k].end()){
						// if they overlap: add the second cluster to the first: look at each point of the second cluster:
						theyHaveOverlap = true;
						for (unsigned int l = 0; l < clusters[k].size(); ++l){
							// check if this point is already contained in the first cluster
							if (find(clusters[i].begin(), clusters[i].end(), clusters[k][l]) == clusters[i].end()){
								// if it is not contained in the first cluster: add it
								clusters[i].push_back(clusters[k][l]);
							}
						}
					}
					// if the two checked clusters overlapped and are now fused to the first:
					if (theyHaveOverlap){
						// delete the second
						clusters.erase(clusters.begin()+k);
						// since we have erased the old cluster at position k,
						// we need to check the new cluster at position k in the next iteration
						--k;
					}
				}
			}
		}
	}
}

std::ostream& operator<<(std::ostream& lhs, SG::cluster const& rhs)
{
	lhs << rhs.relevances.size() << " ";
	for (unsigned int i = 0; i < rhs.relevances.size(); ++i)
	{
		lhs << rhs.relevances[i] << " ";
	}
	lhs << rhs.gridpoints.size() << " ";
	for (unsigned int i = 0; i < rhs.gridpoints.size(); ++i)
	{
		lhs << rhs.gridpoints[i] << " ";
	}
	return lhs;
}

void writeClusters(std::string fileName, std::vector<SG::cluster> clusters)
{
	std::cout << "Writing to " << fileName << std::endl;
	std::ofstream raus;
	raus.open(fileName.c_str());
	if ( (! raus.is_open()) || ( ! raus.good()) )
	{
		cerr << "Failed to open file " << fileName << " for loading of clusters." << std::endl;
		throw -1;
	}
	for (unsigned int i = 0; i < clusters.size(); ++i)
	{
		raus << clusters[i] << std::endl;
		if ( ! raus.good() )
		{
			cerr << "Failed to write cluster " << i << "/" << clusters.size() << " to file " << fileName << std::endl;
			throw -1;
		}
	}
	raus.close();
}


std::istream& operator>>(std::istream& lhs, SG::cluster& rhs)
{
	rhs.relevances.clear();
	rhs.gridpoints.clear();

	std::string temp;
	lhs >> temp;
	if ((! lhs.good()) || (temp == ""))
	{
		return lhs;
	}
	unsigned int number = boost::lexical_cast<unsigned int>(temp);
	for (unsigned int i = 0; i < number; ++i)
	{
		lhs >> temp;
		double value = boost::lexical_cast<double>(temp);
		rhs.relevances.push_back(value);
	}
	lhs >> temp;
	number = boost::lexical_cast<unsigned int>(temp);
	for (unsigned int i = 0; i < number; ++i)
	{
		lhs >> temp;
		unsigned int value = boost::lexical_cast<unsigned int>(temp);
		rhs.gridpoints.push_back(value);
	}
	return lhs;
}

void loadClusters(std::string fileName, std::vector<SG::cluster>& vec)
{
// 	cout << "Loading clusters from file " << fileName << endl;
	ifstream rein(fileName.c_str());
	if ( (! rein.is_open()) || ( ! rein.good()) )
	{
		cerr << "Failed to open file " << fileName << " for loading of clusters." << endl;
		throw -1;
	}
	vec.clear();
	while (rein.good())
	{
		SG::cluster temp;
		rein >> temp;
		if (rein.good())
		{
			vec.push_back(temp);
		}
	}
	rein.close();
// 	cout << "Loaded " << vec.size() << " clusters." << endl;
}


double sumTopN(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster, unsigned int n)
{
	if (n == 0){
		return 0.0;
	}
	if (n > finalCluster.gridpoints.size()){
		n = finalCluster.gridpoints.size();
	}
	double relevanceVal = 0.0;
	std::vector<double> allRelevances;
	for (std::vector<unsigned int>::const_iterator it = finalCluster.gridpoints.begin(); it != finalCluster.gridpoints.end(); ++it){
		allRelevances.push_back(grid[*it].proximity_);
	}
	std::sort(allRelevances.begin(), allRelevances.end());
	for (unsigned int i = 0; i < n; ++i){
		relevanceVal += allRelevances[allRelevances.size()-1-i];
	}
	return relevanceVal;
}

// Computes the relevance of a cluster in L1 metric
double relevanceSum(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster)
{
	return sumTopN(grid, finalCluster, finalCluster.gridpoints.size());
}

// Computes the relevance of a cluster in L inf metric (maximum)
double relevanceMax(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster)
{
	return sumTopN(grid, finalCluster, 1);
}

// Computes the relevance of a cluster as average of the top N contained gridpoints
double relevanceAverageOfTopN(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster, unsigned int n)
{
	return sumTopN(grid, finalCluster, n)/n;
}

// Computes the relevance of a cluster as average metric
double relevanceSize(SG::cluster& finalCluster)
{
	return finalCluster.gridpoints.size();
}

// Computes the relevance of a cluster as average metric
double relevanceAverage(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster)
{
	return relevanceSum(grid, finalCluster)/finalCluster.gridpoints.size();
}

void computeRelevance(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster)
{
	finalCluster.relevances.push_back( relevanceSum(grid, finalCluster) );
	finalCluster.relevances.push_back( relevanceAverage(grid, finalCluster) );
	finalCluster.relevances.push_back( relevanceMax(grid, finalCluster) );
	finalCluster.relevances.push_back( finalCluster.gridpoints.size() );

	for (unsigned int i = 0; i < finalCluster.relevances.size(); ++i)
	{
		if (std::isnan(finalCluster.relevances[i]) || std::isinf(finalCluster.relevances[i]))
		{
			cerr << "Cluster relevance no " << i << " is nan or inf! Setting it to 0.0" << endl;
			finalCluster.relevances[i] = 0.0;
		}
	}
}


void computeRelevance(SG::sphericalGrid<double> const & grid, std::vector<SG::cluster>& finalClusters)
{
	for (unsigned int i = 0; i < finalClusters.size(); ++i)
	{
		computeRelevance(grid, finalClusters[i]);
	}
}




void printEventInfoPerCluster(std::vector<SG::candidate_t> const & events, SG::sphericalGrid<double> const & grid,
		std::vector<std::vector<cluster> > const & finalClustersMulti, std::string outputFolder)
{
//	cerr << "printing events per clusters" << endl;
	for (unsigned int i = 0; i < finalClustersMulti.size(); ++i)
	{
		for (unsigned int j = 0; j < finalClustersMulti[i].size(); ++j)
		{
			const double minRel = 0.05;
			double rel = sgcu::maxValue(finalClustersMulti[i][j].relevances);
			// cout << "maxRel= " << rel << endl;
			if (rel > minRel)
			{
				std::vector<SG::candidate_t> eventsInCluster;
				getEventsInCluster(events, grid, finalClustersMulti[i][j], eventsInCluster);
				ofstream raus(string(outputFolder+"info/EventFromSegmentation" +
						boost::lexical_cast<string>(i)+"Cluster"+boost::lexical_cast<string>(j) + ".dat").c_str());
				raus << eventsInCluster.size() << endl;
				for (unsigned int k = 0; k < eventsInCluster.size(); ++k)
				{
					raus << eventsInCluster[k] << endl;
				}
				raus.close();
				ofstream raus2(string(outputFolder+"info/EventtimesFromSegmentation" +
						boost::lexical_cast<string>(i)+"Cluster"+boost::lexical_cast<string>(j) + ".dat").c_str());
				raus2.precision(20);
				for (unsigned int k = 0; k < eventsInCluster.size(); ++k)
				{
					raus2 << eventsInCluster[k].time_ << endl;
				}
				raus2.precision(6);
				raus2.close();
			}
		}
	}
}

void getEventsInCluster(const std::vector<SG::candidate_t>& events, SG::sphericalGrid<double> const & grid, const SG::cluster& cluster,
		std::vector<SG::candidate_t>& eventsInCluster)
{
	const double togetherDistance = 0.25*M_PI/180.0;
	eventsInCluster.clear();
	for (unsigned int i = 0; i < cluster.gridpoints.size(); ++i)
	{
		SG::candidate<double> point = grid[cluster.gridpoints[i]];
		for (unsigned int j = 0; j < events.size(); ++j)
		{
			if (distanceRad(point, events[j]) <= togetherDistance)
			{
				eventsInCluster.push_back(events[j]);
			}
		}
	}
}



void renormalizeClusterRelevances(std::vector<std::vector<SG::cluster> > & finalClustersMulti)
{
//	std::cerr << "DEBUG renormalizeStep3 started" << endl;
	if (finalClustersMulti.size() == 0)
	{
		cerr << "renormalizeStep3 is not possible on empty data" << endl;
		return;
	}

	unsigned int numberSegmentations = finalClustersMulti.size();
	string normalizationFile = "intermediate/step4/normalizationDataStep3/normalizationStep3_Maxima.txt";
	std::vector<double> comparisonVals;
	sgcu::loadVectorFromFile(normalizationFile, comparisonVals);
	if (comparisonVals.size() == 0)
	{
		cerr << "renormalizeStep3 loaded empty file " << normalizationFile << ". Unable to renormalize." << endl;
		return;
	}

	struct rememberBoostHistogram bestFit;
	//getFitHistogram(comparisonVals, bestFit);
	getFitHistogramOwn(comparisonVals, bestFit);

	for (unsigned int j = 0; j < numberSegmentations; ++j)
	{
		for (unsigned int i = 0; i < finalClustersMulti[j].size(); ++i)
		{
			for (unsigned int k = 0; k < finalClustersMulti[j][i].relevances.size(); ++k)
			{
				double oldValue = finalClustersMulti[j][i].relevances[k];
				double finalRelevance = getRelevanceFromFit(oldValue, bestFit);
				finalClustersMulti[j][i].relevances[k] = finalRelevance;
			}
		}
	}
}



// mark the gridpoints as belonging to the clusters
void drawClusters(std::vector<cluster>& finalClusters, SG::sphericalGrid<double> grid)
{
	// set everything to "no cluster"
	grid.reset();

	// go through all clusters
	for (unsigned int i = 0; i < finalClusters.size(); ++i)
	{
		// for each cluster mark every point in the cluster
		for (unsigned int j = 0; j < finalClusters[i].gridpoints.size(); ++j)
		{
			unsigned int index = finalClusters[i].gridpoints[j];
			if (index >= grid.size())
			{
				std::cerr << "Drawing gridpoint " << index << " (of max. " << grid.size() << ") for finalCluster " << i << " gridpoint " << j << endl;
			}
			if (grid[index].proximity_ != 0.0)
			{
				cerr << "gridpoint " << index << " is marked multiple times as gridpoint!" << endl;
			}
			// set every point to the relevance of this cluster
			// store the value of the most unlikely metric as proximity
			grid[index].proximity_ = *std::max_element(finalClusters[i].relevances.begin(), finalClusters[i].relevances.end());
			// store which metric has been the most unlikely as multipurpose
			grid[index].multiPurposeValue_ = std::max_element( finalClusters[i].relevances.begin(), finalClusters[i].relevances.end() ) - finalClusters[i].relevances.begin();
		}
	}
}

void findClustersFromGrid(SG::sphericalGrid<double> const & grid, std::vector<std::vector<std::pair<unsigned int, unsigned int> > >& clusters)
{
	if (grid.size() == 0){
		cerr << "There are no clusters to find in an empty grid." << endl;
		return;
	}
	clusters.clear();

	std::vector<searchpoint> searchpoints;
	for (unsigned int i = 0; i < grid.size(); ++i)
	{
		if (grid[i].proximity_ != 0.0)
		{
			searchpoints.push_back( SG::searchpoint(grid[i].proximity_, i, 0) );
		}
	}

	std::sort(searchpoints.begin(), searchpoints.end(), wayToSortSearchpoints);

// 	cout << "Found " << searchpoints.size() << " potential cluster seeds" << endl;

	unsigned int maxNum = 15000;
	if (searchpoints.size() > maxNum){
		searchpoints.resize(maxNum);
// 		cout << "Only using top " << searchpoints.size() << " potential cluster seeds" << endl;
	}

	for (unsigned int i = 0; i < searchpoints.size(); ++i)
	{
		unsigned int maxIndexGrid = searchpoints[i].indexGrid;
		unsigned int maxIndexDist = 0;

		// process one maximum (=cluster) and the points around it
		vector<pair<unsigned int, unsigned int> > cluster;
		set<pair<unsigned int, unsigned int> > convexHull;

		// check if this gridpoint is already part of any previous cluster
		pair<unsigned int, unsigned int> center = pair<unsigned int, unsigned int>(maxIndexDist, maxIndexGrid);
		bool alreadyIntegrated = false;
		for (unsigned int j = 0; j < clusters.size(); ++j){
			if (find(clusters[j].begin(), clusters[j].end(), center) != clusters[j].end()){
				alreadyIntegrated = true;
				break;
			}
		}

		// if it is not part of a previous cluster, we start a new one
		if (! alreadyIntegrated){
			// start the cluster
			cluster.push_back(center);

			// the neighbours need to be found to check if they should become part of the cluster
			addNeighboursToConvexHullGrid(grid, center, convexHull);

			// check the convex hull for searchpoints which are also part of the foreground/signal/have a proximity != 0.0
			checkClusterForSuperSimpleConvexHullExpansionGrid(grid, cluster, convexHull);

			// now this cluster is done
			clusters.push_back(cluster);
		}
	}
}

void checkClusterForSuperSimpleConvexHullExpansionGrid(SG::sphericalGrid<double> const & grid,
		std::vector<std::pair<unsigned int, unsigned int> >& cluster, std::set<std::pair<unsigned int, unsigned int> >& convexHull)
{
	// check every gridpoint in the convex hull
	for (set<pair<unsigned int, unsigned int> >::iterator currentGridpoint = convexHull.begin(); currentGridpoint != convexHull.end(); ++currentGridpoint)
	{
		// see if this gridpoint is part of the segmented foreground (= its proximity is not 0)
		if (grid[currentGridpoint->second].proximity_ != 0.0)
		{
			// if we haven't added this point to the cluster already
			if (find(cluster.begin(), cluster.end(), *currentGridpoint) == cluster.end())
			{
				// then add the gridpoint to the cluster
				cluster.push_back(*currentGridpoint);
				// and expand the convex hull for the neighbours of the new gridpoint
				addNeighboursToConvexHullGrid(grid, *currentGridpoint, convexHull);
			}
		}
	}
}

void addNeighboursToConvexHullGrid(SG::sphericalGrid<double> const & grid,
		const std::pair<unsigned int, unsigned int>& center, std::set<std::pair<unsigned int, unsigned int> >& convexHull)
{
	// the grid is just a helping construct here to find the neighbours since this function is already present
	static SG::sphericalGrid<double> tempGrid = grid;
	std::vector<SG::candidate_t>::iterator gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
	gC = tempGrid.begin()+center.second;
	tempGrid.getNeighbours(gC, gNW, gN, gNE, gW, gE, gSW, gS, gSE);

	// add the neighbours on this searchBin
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gNW-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gN-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gNE-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gW-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gE-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gSW-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gS-tempGrid.begin()));
	convexHull.insert(pair<unsigned int, unsigned int>(center.first, gSE-tempGrid.begin()));
}

void printClustersToASCII(std::vector<std::vector<cluster> > const & finalClustersMulti, SG::sphericalGrid<double> const & grid)
{
	for (unsigned int n = 0; n < finalClustersMulti.size(); ++n)
	{
		for (unsigned int i = 0; i < finalClustersMulti[n].size(); ++i)
		{
			double maxRel = sgcu::maxValue(finalClustersMulti[n][i].relevances);
			if (maxRel > 0.0)
			{
				std::string fileName = "info2/segmentation" + boost::lexical_cast<std::string>(n) + "cluster" +  boost::lexical_cast<std::string>(i) + ".dat";
				ofstream raus(fileName.c_str());
				raus << "cluster significance: " << maxRel << endl;
				raus << finalClustersMulti[n][i].gridpoints.size() << " gridpoints, format: declination rightAscension" << endl;
				for (unsigned int j = 0; j < finalClustersMulti[n][i].gridpoints.size(); ++j)
				{
					raus << grid[finalClustersMulti[n][i].gridpoints[j]].zenith_ << " " << grid[finalClustersMulti[n][i].gridpoints[j]].azimuth_  << endl;
				}
				raus.close();
			}
		}
	}
}

double getMaxRelevance(std::vector<std::vector<SG::cluster> > const & finalClustersMulti, int segmentation)
{
	double theTotalFinalMax = 0.0;
	for (unsigned int n = 0; n < finalClustersMulti.size(); ++n)
	{	// all segmentations
		if ((segmentation >= 0) && (segmentation != (int)n)) continue;
		for (unsigned int i = 0; i < finalClustersMulti[n].size(); ++i)
		{	// all clusters
			for (unsigned int j = 0; j < finalClustersMulti[n][i].relevances.size(); ++j)
			{	// all relevances
				if (finalClustersMulti[n][i].relevances[j] > theTotalFinalMax)
					theTotalFinalMax = finalClustersMulti[n][i].relevances[j];
			}
		}
	}
	return theTotalFinalMax;
}


void printMaxRelevanceSegmentationwise(std::vector<std::vector<SG::cluster> > const & finalClustersMulti,
		SG::sphericalGrid<double> const & grid, std::vector<double> const & segmentationPoints)
{
//	cerr << "DEBUG signalFirstInternal::printMaxRelevanceSegmentationwise : finalclustersmulti.size() = " << finalClustersMulti_.size() << endl;
	for (unsigned int n = 0; n < finalClustersMulti.size(); ++n)
	{
		double theSegmentationFinalMax = 0.0;
		unsigned int where = 0;
		double declination = 0.0;
		double rightAscension = 0.0;
		for (unsigned int i = 0; i < finalClustersMulti[n].size(); ++i)
		{
			for (unsigned int j = 0; j < finalClustersMulti[n][i].relevances.size(); ++j)
			{
				if (finalClustersMulti[n][i].relevances[j] > theSegmentationFinalMax)
				{
					theSegmentationFinalMax = finalClustersMulti[n][i].relevances[j];
					where = j;
					declination = grid[finalClustersMulti[n][i].gridpoints[0]].zenith_;
					rightAscension = grid[finalClustersMulti[n][i].gridpoints[0]].azimuth_;
				}
			}
		}
		double rightAscensionCorrected = rightAscension;
		if (rightAscensionCorrected < 0.0){
			rightAscensionCorrected += 2.0*M_PI;
		}
		cout << "Segmentation" << n << " ( " << segmentationPoints[n] << " ) maximal significance: " << theSegmentationFinalMax <<
				" cluster location (Dec=" << declination/M_PI*180.0 << ", Ra=" << rightAscensionCorrected/M_PI*180.0 << ", Org. Ra=" <<
				rightAscension/M_PI*180.0 << ")  @ metric " << where << endl;
	}
}



// Renormalizes a set of finalclustersMulti (all segmentations)
// against the (basic or extended, unnormalized) clusters from normalization2 files
void renormalizeStep2MaximaOnly(std::vector<std::vector<SG::cluster> >& finalClustersMulti, std::vector<double> const & segmentationPoints,
		std::string folderName, std::string searchString1)
{
	if (finalClustersMulti.size() == 0)
	{
		cerr << "renormalizeStep2MaximaOnly is not possible on empty data" << endl;
		return;
	}
	if (segmentationPoints.size() != finalClustersMulti.size())
	{
		cerr << "renormalizeStep2MaximaOnly does not have the same number of segmentations for final clusters ( " << finalClustersMulti.size() <<
				" ) than for normalization sets ( " << segmentationPoints.size() << " ). This doesn't match, not renormalizing anything." << endl;
		return;
	}

	cout << "Renormalize found clusters against distributiuons from pseudo experiments." << endl;

	// Load all comparison events
	std::vector<std::vector<struct rememberBoostHistogram > > histsAllSegmentationsAllRelevances;
	std::vector<std::vector<std::vector<SG::cluster> > > allEntries;
	getComparisonValuesAndHistograms(folderName, searchString1, segmentationPoints, allEntries, histsAllSegmentationsAllRelevances);

	for (unsigned int n = 0; n < segmentationPoints.size(); ++n)
	{
		cout << "Renormalizing " << finalClustersMulti[n].size() << " clusters with searchstrings " << searchString1 <<
				" and " << boost::lexical_cast<std::string>(segmentationPoints[n]) << endl;
		renormalizeCore(finalClustersMulti[n], histsAllSegmentationsAllRelevances[n]);
	}
}

// Renormalizes a set of finalclusters (one segmentation)
// against the (basic or extended, unnormalized) clusters comparisonclusters
void renormalizeCore(std::vector<SG::cluster>& finalClusters, std::vector<struct rememberBoostHistogram > const & comparisonHists)
{
	if (finalClusters.size() == 0){
		cerr << "No cluster to renormalize! Not doing anything here!" << endl;
		return;
	}
	if (comparisonHists.size() == 0){
		cerr << "No cluster to renormalize against! Not doing anything here!" << endl;
		return;
	}
	if (finalClusters[0].relevances.size() == 0){
		cerr << "No metric to renormalize! Not doing anything here!" << endl;
		return;
	}
	if (comparisonHists.size() != finalClusters[0].relevances.size()){
		cerr << "The number of metrics of the found clusters (" << finalClusters[0].relevances.size() << ") does not match the number of metrics to compare against (" << comparisonHists.size() << ")! Not doing anything here!" << endl;
		return;
	}

	// evaluate one finalCluster at a time
	for (unsigned int i = 0; i < finalClusters.size(); ++i)
	{
		// cerr << "processing cluster "<< i << " of " << finalClusters.size() << endl;

		// Look at the how this cluster scored at the different metrics for relevance
		std::vector<double> relevances = finalClusters[i].relevances;
		for (unsigned int l = 0; l < relevances.size(); ++l)
		{
			// Test the statistics for the found relevance in this metric against that relevances found for clusters similar to the current cluster in all other metrics
			// Store the corrected relevance for this cluster and this metric
			finalClusters[i].relevances[l] = getRelevanceFromFit(relevances[l], comparisonHists[l]);
		}	// end of loop over all relevance metrics
	}	// end of loop over all found finalClusters
}

double getPValueFromFit(double val, struct rememberBoostHistogram const & fit)
{
	double probabilitySum = 0.0;
	for (unsigned int i = 0; i < fit.size(); ++i)
	{
		if (val > fit.indices[i])
		{
			probabilitySum += fit.values[i];
		}
	}
	double pValue = 1.0-probabilitySum;
	const double minPValue = 1.0/390682215445.0; // sevenSigmaPValue
	if (pValue < minPValue)
	{
		pValue = minPValue;
	}
	return pValue;
}

double getRelevanceFromFit(double val, struct rememberBoostHistogram const & fit)
{
	double pValue = getPValueFromFit(val, fit);
	double probabilitySum = 1.0-pValue;
	return boost::math::erf_inv(probabilitySum) * sqrt(2.0);
}


void condenseMaxima(std::vector<SG::cluster> const & comparisonClustersOneRun, SG::cluster& clusterusMaximus)
{
	if (comparisonClustersOneRun.size() == 0)
	{
		cerr << "Cannot condense maxima of zero clusters." << endl;
		return;
	}
	clusterusMaximus.relevances.resize(comparisonClustersOneRun[0].relevances.size());
	for (unsigned int k = 0; k < clusterusMaximus.relevances.size(); ++k)
	{
		clusterusMaximus.relevances[k] = -9999999.99;
	}
	for (unsigned int j = 0; j < comparisonClustersOneRun.size(); ++j)
	{
		for (unsigned int k = 0; k < clusterusMaximus.relevances.size(); ++k)
		{
			if (clusterusMaximus.relevances[k] < comparisonClustersOneRun[j].relevances[k])
			{
				clusterusMaximus.relevances[k] = comparisonClustersOneRun[j].relevances[k];
			}
		}
	}
}

std::vector<cluster> clustersFromGrid(SG::sphericalGrid<double> const & grid)
{
	std::vector<std::vector<std::pair<unsigned int, unsigned int> > > clusters;

	// Find the clusters
	findClustersFromGrid(grid, clusters);
// 	cout << "Found " << clusters_.size() << " clusters." << endl;

	// Fuse those which were created twice because the neighbour function for convex hulls has aliasing
	fuseAlreadyOverlappingClusters(clusters);
	cout << "Finally found " << clusters.size() << " clusters." << endl;

	return clustersToClusterFormat(grid, clusters);
}

std::vector<cluster> clustersToClusterFormat(SG::sphericalGrid<double> const & grid, std::vector<std::vector<std::pair<unsigned int, unsigned int> > >& clusters)
{
	std::vector<cluster> clustersWithFormat(clusters.size());
	for (unsigned int i = 0; i < clusters.size(); ++i)
	{
		for (unsigned int j = 0; j < clusters[i].size(); ++j)
		{
			clustersWithFormat[i].gridpoints.push_back(clusters[i][j].second);
		}
	}
	computeRelevance(grid, clustersWithFormat);
	return clustersWithFormat;
}





void listDirectory(std::string folderName, std::string searchString1, std::string searchString2, std::vector<std::string>& fileNames){
	fileNames.clear();
	DIR *folder;
	struct dirent *datei;
	folder = opendir(folderName.c_str());
	if (folder == NULL)
	{
		cerr << "Failed to open folder " << folderName << endl;
		throw -8;
	}
	while ( (datei = readdir(folder)) ) {
		std::string entry = std::string(datei->d_name);
		bool useThis = true;
		if (searchString1 != ""){
			if (entry.find(searchString1) == string::npos){
				useThis = false;
			}
		}
		if (searchString2 != ""){
			if (entry.find(searchString2) == string::npos){
				useThis = false;
			}
		}
		if (useThis == true){
			fileNames.push_back(folderName + "/" + std::string(datei->d_name));
		}
	}
}


void addClustersForOneSegmentationAndOnePseudoExperiment(std::vector<std::vector<std::string> > const & fileNames,
		unsigned int currentThresholdIndex, unsigned int currentFilenameIndex,
		std::vector<std::vector<std::vector<SG::cluster> > > & allEntries, std::vector<std::vector<double> > & relevanceComparisonVals)
{
	// load the clusters from one pseudo experiment and one segmentation threshold (contained in one file)
	std::vector<SG::cluster> comparisonClustersOneFile;
	try
	{
		loadClusters(fileNames[currentThresholdIndex][currentFilenameIndex], comparisonClustersOneFile);
	}
	catch (...)
	{
		cerr << "Failed to load clusters for segmentation " << currentThresholdIndex << ". Skipping this segmentation." << endl;
		return;
	}
	if (comparisonClustersOneFile.size() == 0)
	{
		cerr << "No cluster loaded from " << fileNames[currentThresholdIndex][currentFilenameIndex] << endl;
		return;
	}

	// TODO: check if allEntries is not needed! maybe the collection of all clusterusMaximus would do the same trick!

	// add these clusters to the total number of clusters found for this segmentation
	allEntries[currentThresholdIndex].push_back(comparisonClustersOneFile);

	// get one cluster that, for each metric, contains the maximal relevance values of this metric of any considered cluster
	SG::cluster clusterusMaximus;
	condenseMaxima(comparisonClustersOneFile, clusterusMaximus);

	// if relevanceComparisonVals doesn't have the right size (the number of relevance metrics) then it is adopted here
	relevanceComparisonVals.resize( clusterusMaximus.relevances.size() );

	// for each metric: store the maximal value of any cluster
	for (unsigned int j = 0; j < clusterusMaximus.relevances.size(); ++j)
	{
		relevanceComparisonVals[j].push_back(clusterusMaximus.relevances[j]);
	}
}

void getComparisonValuesAndHistograms(std::string folderName, std::string searchString1, std::vector<double> segmentationPoints,
		std::vector<std::vector<std::vector<SG::cluster> > >& allEntries,
		std::vector<std::vector<struct rememberBoostHistogram > >& histsAllSegmentationsAllRelevances)
{
	cout << "Loading comparison values to histograms for folder " << folderName << " and searchString " << searchString1 << endl;

	std::vector<std::vector<std::string> > fileNames;
	fileNames.resize(segmentationPoints.size());

	allEntries.clear();
	allEntries.resize(segmentationPoints.size());

	histsAllSegmentationsAllRelevances.clear();
	histsAllSegmentationsAllRelevances.resize(segmentationPoints.size());

	// TO DO: allow readin of already computed histograms

	unsigned int lastSize = 0;
	for (unsigned int n = 0; n < segmentationPoints.size(); ++n)
	{
		cout << "Loading files for segmentation " << n+1 << "/" << segmentationPoints.size() << " (" << segmentationPoints[n] << ")" << endl;
		std::vector<std::vector<double> > relevanceComparisonVals;
		string searchStringSegmentationThreshold = boost::lexical_cast<std::string>(segmentationPoints[n]);

		// Get a list of all files in the folder that match the current configuration (threshold,
		listDirectory(folderName, searchString1, searchStringSegmentationThreshold, fileNames[n]);
		cout << "Found " << fileNames[n].size() << " matching files." << endl;

		// Now load the clusters from all found pseudo experiments
		bool somethingFailed = false;
		for (unsigned int i = 0; i < fileNames[n].size(); ++i)
		{
			try
			{
				// try to load clusters for one pseudo experiment
				// cerr << i << ": file " << fileNames[n][i] << endl;
				addClustersForOneSegmentationAndOnePseudoExperiment(fileNames, n, i, allEntries, relevanceComparisonVals);
			}
			catch (std::exception& e)
			{
				// Do not crash instantly if a file cannot be read, as multiple might be damaged and that way they are found all at once
				cerr << "Caught exception " << e.what() << " while loading from file " << fileNames[n][i] << endl;
				somethingFailed = true;
			}
		}
		if (somethingFailed)
		{
			throw(-1);
		}

		// All clusters matching the desired segmentation threshold from all pseudo experiments have been read here
		if (n > 0)
		{
			// check if the clusters had the same number of metrics as last time
			if (relevanceComparisonVals.size() != lastSize)	// all thresholds must have the same number of metrics
			{
				cerr << "Found clusters with " << relevanceComparisonVals.size() << " metrics for segmentation " << n << " but expected " << lastSize << endl;
				throw(-1);
			}
		}
		lastSize = relevanceComparisonVals.size();

		// Now make histograms from the loaded clusters
		computeAndOutputHistogramsFromClusters(folderName, n, relevanceComparisonVals, searchStringSegmentationThreshold, histsAllSegmentationsAllRelevances);
	}
}













} // end of namespace SG
