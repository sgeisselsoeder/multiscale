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
	
#ifndef SGCLUSTERING_H
#define SGCLUSTERING_H

#include <vector>
#include <string>

#include "candidate.hpp"
#include "sphericalGrid.hpp"
#include "allScales.h"
#include "vtkWriter.h"
#include "searchpoint.h"
#include "rememberBoostHistogram.h"

namespace SG{

struct cluster
{
	std::vector<double> relevances;
	std::vector<unsigned int> gridpoints;
};

void computeRelevance(SG::sphericalGrid<double> const & grid, std::vector<SG::cluster>& finalClusters);

std::ostream& operator<<(std::ostream& lhs, SG::cluster const& rhs);
std::istream& operator>>(std::istream& lhs, SG::cluster& rhs);
void writeClusters(std::string fileName, std::vector<SG::cluster> clusters);
void loadClusters(std::string fileName, std::vector<SG::cluster>& vec);
void printEventInfoPerCluster(std::vector<SG::candidate_t> const & events, SG::sphericalGrid<double> const & grid,
		std::vector<std::vector<cluster> > const & finalClustersMulti, std::string outputFolder);
void printClustersToASCII(std::vector<std::vector<cluster> > const & finalClustersMulti, SG::sphericalGrid<double> const & grid);
void condenseMaxima(std::vector<SG::cluster> const & comparisonClustersOneRun, SG::cluster& clusterusMaximus);
void renormalizeStep2MaximaOnly(std::vector<std::vector<SG::cluster> >& finalClustersMulti, std::vector<double> const & segmentationPoints,
		std::string folderName, std::string searchString1 = "clusters_segment");
void printMaxRelevanceSegmentationwise(std::vector<std::vector<SG::cluster> > const & finalClustersMulti, SG::sphericalGrid<double> const & grid, std::vector<double> const & segmentationPoints);
void getEventsInCluster(const std::vector<SG::candidate_t>& events, SG::sphericalGrid<double> const & grid,
		const SG::cluster& cluster, std::vector<SG::candidate_t>& eventsInCluster);
void addNeighboursToConvexHullGrid(SG::sphericalGrid<double> const & grid,
		const std::pair<unsigned int, unsigned int>& center, std::set<std::pair<unsigned int, unsigned int> >& convexHull);
void checkClusterForSuperSimpleConvexHullExpansionGrid(SG::sphericalGrid<double> const & grid,
		std::vector<std::pair<unsigned int, unsigned int> >& cluster, std::set<std::pair<unsigned int, unsigned int> >& convexHull);
void renormalizeCore(std::vector<SG::cluster>& finalClusters, std::vector<struct rememberBoostHistogram > const & comparisonHists);
std::vector<cluster> clustersToClusterFormat(SG::sphericalGrid<double> const & grid,
		std::vector<std::vector<std::pair<unsigned int, unsigned int> > >& clusters);
double getRelevanceFromFit(double val, struct rememberBoostHistogram const & fit);

void getComparisonValuesAndHistograms(std::string folderName, std::string searchString1, std::vector<double> segmentationPoints,
		std::vector<std::vector<std::vector<SG::cluster> > >& allEntries,
		std::vector<std::vector<struct rememberBoostHistogram > >& histsAllSegmentationsAllRelevances);




void drawClusters(std::vector<cluster>& finalClusters, SG::sphericalGrid<double> grid);
	
void fuseAlreadyOverlappingClusters(std::vector<std::vector<std::pair<unsigned int, unsigned int> > >& clusters);
	
void getComparisonValuesAndHistograms(std::string folderName, std::string searchString1, std::vector<std::vector<std::vector<SG::cluster> > >& allEntries, std::vector<std::vector<struct rememberBoostHistogram > >& histsAllSegmentationsAllRelevances);

void getFitHistogramOwn(const std::vector<double>& data, struct rememberBoostHistogram& fit);
void getFitHistogramOwn(const struct rememberBoostHistogram& data, struct rememberBoostHistogram& fit);
double getMaxRelevance(std::vector<std::vector<SG::cluster> > const & finalClustersMulti, int segmentation = -1);

void printMaxRelevanceSegmentationwise(std::vector<std::vector<SG::cluster> > const & finalClustersMulti,
		SG::sphericalGrid<double> const & grid, std::vector<double> const & segmentationPoints);

void renormalizeClusterRelevances(std::vector<std::vector<SG::cluster> > & finalClustersMulti);

// Computes the relevance of a cluster in L1 metric
double relevanceSum(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster);
// Computes the relevance of a cluster in L inf metric (maximum)
double relevanceMax(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster);
// Computes the relevance of a cluster as average of the top N contained gridpoints
double relevanceAverageOfTopN(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster, unsigned int n);
// Computes the relevance of a cluster as average metric
double relevanceSize(SG::cluster& finalCluster);
// Computes the relevance of a cluster as average metric
double relevanceAverage(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster);
// Computes the relevance of a cluster in all available metrics
void computeRelevance(SG::sphericalGrid<double> const & grid, SG::cluster& finalCluster);

std::vector<cluster> clustersFromGrid(SG::sphericalGrid<double> const & grid);
	
} // end of namespace SG

#endif 
