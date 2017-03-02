/*
 * sgeisselsoeder
 *
 * Date June 2016
 *
 */
 
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <limits>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <boost/lexical_cast.hpp>
#include "candidate.hpp"
#include "vtkWriter.h"
#include "allScales.h"
#include <omp.h>
#include <boost/foreach.hpp>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cmath>
#include "sgcu/sgcu.hpp"
#include "rememberBoostHistogram.h"

#include <cmath>	// for isnan

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::stringstream;

// implementation of struct aroundRatios

namespace SG{



//double computeSegmentationThreshold(std::vector<double>& valuesToSegment, unsigned int numberOfSmearings = 50, double useOnlyUpperFraction = 0.5);

// NOTE this function can be kept for now to segment the individual scales as it worked well,
// but it could be replaced by a more deterministic segmentation as will be the case for the cluster search in the end

// This function looks at the distribution and tries to separate the majority of random non clustering gridpoints
// from the parts where high clustering might occur
// it fits the falling flank of the distribution linearly, computes the intersection with y=0
// and treats everything below this as background
double computeSegmentationThreshold(std::vector<double>& valuesToSegment, unsigned int numberOfSmearings = 50, double upperFractionToUse = 0.5)
{
	static unsigned int counter = 0;
//	if (counter == 0)
//	{
//		// build a histogram for the densities of all gridpoints
//		// create an accumulator
//		// 100 bins, approx. the same range as the actual distribution
//		acc myAccumulator( tag::density::num_bins = 100, tag::density::cache_size = valuesToSegment.size()-2);
//		// fill accumulator
//		for (unsigned int j = 0; j < valuesToSegment.size(); ++j){
//			myAccumulator(valuesToSegment[j]);
//		}
//		// do the boost histogram
//		histogram_type hist = density(myAccumulator);
//		ofstream raus(std::string(outputFolder_+"segmentationHistBeforeFull" + boost::lexical_cast<std::string>(counter) + ".txt").c_str());
//		for( unsigned int i = 0; i < (unsigned int) hist.size(); i++ ) {
//			raus << hist[i].first << " " << hist[i].second << endl;
//		}
//		raus.close();
//	}

	std::vector<double> valuesCopy = valuesToSegment;
	std::sort(valuesCopy.begin(), valuesCopy.end());
	valuesToSegment.clear();
	int numToUse = upperFractionToUse*valuesCopy.size();
	valuesToSegment.insert(valuesToSegment.begin(), valuesCopy.end()-numToUse, valuesCopy.end());

	// the threshold for background against signal in the end
	double xThres = 0.0;

	double mean = 0.0;
	double stdDev = 0.0;
	for (unsigned int j = 0; j < valuesToSegment.size(); ++j)
		mean += valuesToSegment[j];
	mean /= valuesToSegment.size();

	// build a histogram for the densities of all gridpoints
	// create an accumulator
	// 100 bins, approx. the same range as the actual distribution
	acc myAccumulator( boost::accumulators::tag::density::num_bins = 100, boost::accumulators::tag::density::cache_size = valuesToSegment.size()-2);
	// fill accumulator
	for (unsigned int j = 0; j < valuesToSegment.size(); ++j){
		myAccumulator(valuesToSegment[j]);
		double temp = valuesToSegment[j] - mean;
		stdDev += temp*temp;
	}
	stdDev /= valuesToSegment.size();
	stdDev = sqrt(stdDev);

// 	cout << "mean +- stddev: " << mean << " " << stdDev << endl;
// 	cout << "1,2,3,4,5 sigma borders: " << mean+stdDev << " , " << mean+2*stdDev << " , " << mean+3*stdDev << " , " << mean+4*stdDev << " , " << mean+5*stdDev << endl;
// 	xThres = mean+5*stdDev;
// 	return xThres;

	// do the boost histogram
	histogram_type hist = boost::accumulators::density(myAccumulator);
//	histogram_type histOrg = hist;
//	if (counter == 0)
//	{
//		ofstream raus(std::string(outputFolder_+"segmentationHistBefore" + boost::lexical_cast<std::string>(counter) + ".txt").c_str());
//		for( unsigned int i = 0; i < (unsigned int) hist.size(); i++ ) {
//			raus << hist[i].first << " " << hist[i].second << endl;
//		}
//		raus.close();
//	}

	for (unsigned int i = 0; i < numberOfSmearings; ++i){
		vector<double> vals(hist.size());
		for (unsigned int j = 0; j < (unsigned int) hist.size(); ++j){
			vals[j] = hist[j].second;
		}
		for (unsigned int j = 1; j < (unsigned int) hist.size()-1; ++j){
			hist[j].second = (vals[j-1] + vals[j] + vals[j+1]) / 3.0;
		}
		hist[0].second = (vals[0] + vals[1]) / 2.0;
		hist[hist.size()-1].second = (vals[hist.size()-2] + vals[hist.size()-1]) / 2.0;
	}

//	if (counter == 0)
//	{
//	 	ofstream raus(std::string(outputFolder_+"segmentationHist" + boost::lexical_cast<std::string>(counter) + ".txt").c_str());
//		for( unsigned int i = 0; i < (unsigned int) hist.size(); i++ ) {
//			raus << hist[i].first << " " << hist[i].second << " " << histOrg[i].second << endl;
//		}
//		raus.close();
//	}

	// find the maximum of the distribution and the point where only half of the maximum is
	std::vector<double> vals(hist.size());
	double max = 0.0;
	unsigned int maxIndex = 0;
//	double halfMax = 0.0;
	unsigned int halfMaxIndex = 0;
	for( unsigned int i = 0; i < (unsigned int) hist.size(); i++ ) {
		vals[i] = hist[i].second;
		if (hist[i].second > max){
			max = hist[i].second;
			maxIndex = i;
		}
		if (hist[i].second > 0.5*max){
//			halfMax = hist[i].second;
			halfMaxIndex = i;
		}
	}

	// now average the next broadness bins from the maximum on
	// to get a more robust point for (approx.) the maximum and the half maximum
	const unsigned int broadness = 3;
	double topX = 0.0;
	double topY = 0.0;
	for (unsigned int i = 0; ((i < broadness) && (i < (unsigned int) hist.size())); ++i){
		topX += hist[maxIndex+i].first;
		topY += hist[maxIndex+i].second;
	}
	topX /= broadness;
	topY /= broadness;

	double middleX = 0.0;
	double middleY = 0.0;
	for (unsigned int i = 0; i < (i < broadness) && (i < (unsigned int) hist.size()); ++i){
		middleX += hist[halfMaxIndex+i].first;
		middleY += hist[halfMaxIndex+i].second;
	}
	middleX /= broadness;
	middleY /= broadness;

	// compute the slope between the two points
	double diffX = middleX - topX;
	double diffY = middleY - topY;
	double slope = diffY/diffX;

	// compute the intersection of a line from the maximum to through the half maximum with the x axis
	double dx = fabs(topY/slope);
	xThres = topX+dx;

	++counter;
	return xThres;
}

void allScales::assertNoNansAndInfs() const
{
	unsigned int counter = 0;
	// cross check if there is a nan or inf value somewhere in the grid
	for (unsigned int i = 0; i < aroundRatioAll.size(); ++i)
	{
		for (unsigned int j = 0; j < aroundRatioAll[i].size(); ++j)
		{
			if (std::isnan(aroundRatioAll[i][j]) || std::isinf(aroundRatioAll[i][j]))
			{
				++counter;
			}
		}
	}
	if (counter > 0)
	{
		std::cerr << counter << " ratios have a nan or inf value! This is not how it works!" << std::endl;
		throw -1;
	}
}

unsigned int countNumberOfNonZeros(struct allScales const & allScales)
{
	unsigned int counter  = 0;
	for (unsigned int i = 0; i < allScales.aroundRatioAll.size(); ++i)
	{
		for (unsigned int j = 0; j < allScales.aroundRatioAll[i].size(); ++j)
		{
			if (allScales.aroundRatioAll[i][j] != 0.0)
				++counter;
		}
	}
	return counter;
}


void segmentAllRatiosScalewise(allScales & scales)
{
	unsigned int numPoints = scales.aroundRatioAll.size();
	unsigned int numDistanceBins = scales.aroundRatioAll[0].size();
	if (numPoints == 0)
	{
		cerr << "There is nothing to segment in an empty set of ratios." << endl;
		return;
	}
	if (numDistanceBins == 0)
	{
		cerr << "There is nothing to segment in a set of ratios without distance bins." << endl;
		return;
	}

	for (unsigned int j = 0; j < numDistanceBins; ++j)
	{
		std::vector<double> allRatioValues(numPoints);
		//cout << "######### scale " << j << endl;
		for (unsigned int i = 0; i < numPoints; ++i)
		{
			allRatioValues[i] = scales.aroundRatioAll[i][j];
		}
		double threshold = computeSegmentationThreshold(allRatioValues);
		if (threshold < sgcu::minValue(allRatioValues))
			throw -1;
		if (threshold > sgcu::maxValue(allRatioValues))
			throw -2;
			
		for (unsigned int i = 0; i < numPoints; ++i)
		{
			if (scales.aroundRatioAll[i][j] < threshold)
			{
				scales.aroundRatioAll[i][j] = 0.0;
			}
		}
	}
}





unsigned int allScales::size() const
{
	return aroundRatioAll.size();
}

void allScales::divideRatios(double value)
{
	if (value == 1.0)
		return;
	if (value == 0.0)
	{
		std::cerr << "Division by zero destroys all ratios and the universe!" << std::endl;
		throw (-1);
	}
	double invValue = 1.0/value;
	for (unsigned int i = 0; i < aroundRatioAll.size(); ++i)
	{
		for (unsigned int n = 0; n < aroundRatioAll[0].size(); ++n)
		{
			aroundRatioAll[i][n] *= invValue;
		}
	}
}

allScales::allScales(){
	clear();
}

void allScales::clear(){
	numBins_ = 0;
	stepDegree_ = 0.0;
	stepSize_ = 0.0;

	aroundRatioAll.clear();
}

void allScales::set(double stepDeg, unsigned int numGridpoints)
{
	stepDegree_ = stepDeg;
	numBins_ = static_cast<unsigned int>(90.0/stepDegree_)+1;
	const double oneDegree = M_PI/180.0;
	stepSize_ = stepDegree_*oneDegree;

	if (numGridpoints != 0)
	{
		if (aroundRatioAll.size() != numGridpoints){
			aroundRatioAll.resize(numGridpoints);
			for (unsigned int i = 0; i < aroundRatioAll.size(); ++i){
				aroundRatioAll[i].resize(numBins_);
			}
		}
	}
}

void allScales::init(double stepDeg, unsigned int numGridpoints)
{
 	clear();
 	set(stepDeg, numGridpoints);
}

void setToZero(allScales & scales)
{
	for (unsigned int i = 0; i < scales.aroundRatioAll.size(); ++i)
	{
		for (unsigned int j = 0; j < scales.aroundRatioAll[i].size(); ++j)
		{
			scales.aroundRatioAll[i][j] = 0.0;
		}
	}
}

long double allScales::computePValue(long double mean, long double found){
	long double p = 0.0;
	const long double changeThreshold = 120.0;
	if ((found <= changeThreshold) && (mean <= changeThreshold)){
		// plain poisson way:
		long double exponent = expl(-mean);
		long double power = powl(mean, found);
		long double fac = sgcu::faculty(found);
		int expExp = log10l(exponent);
		int expPow = log10l(power);
		int expFac = log10l(fac);
		int finalExp = expExp+expPow-expFac;
		exponent /= powl(10, expExp);
		power /= powl(10, expPow);
		fac /= powl(10, expFac);
		p = exponent*power/fac;
		p *= powl(10, finalExp);
	} else {
		if (found == 0.0) {
			p = 0.0;
		} else {
			// poissonian with stirling formula:
			long double frac = mean/found;
			long double logarithm = logl(frac);
			long double exponent2 = found*(1.0+logarithm)-mean;
			long double efunc = expl(exponent2);
			long double denominator = sqrtl(2.0*M_PI*(found+1.0/6.0));
			p = efunc / denominator;
		}
	}
	return p;
}


void allScales::computeAllRatios(const SG::sphericalGrid<double>& grid, const std::vector<SG::candidate_t>& events)
{
	timeval starttTotal, endtTotal;
	timeval startt, endt;
	gettimeofday(&starttTotal, 0);
	gettimeofday(&startt, 0);
 	cout << "DEBUG allScales::computeAllRatios started" << endl;

	if (aroundRatioAll.size() != grid.size())
	{
		aroundRatioAll.resize(grid.size());
		for (unsigned int i = 0; i < aroundRatioAll.size(); ++i)
		{
			aroundRatioAll[i].resize(numBins_);
		}
	}
	
	unsigned int eventsSize = events.size();
	int gridSize = static_cast<int>(grid.size());
	double evtSize = static_cast<double>(events.size());
	
	gettimeofday(&endt, 0); cout << "Runtime of ratio resize: " << endt.tv_sec - startt.tv_sec << endl; gettimeofday(&startt, 0);

	double distanceLimit = static_cast<double>(numBins_)*stepSize_;

//	double addValue = 1.0/static_cast<double>(events.size());
//	omp_set_num_threads(4);
//	#pragma omp parallel for // schedule(dynamic,1)	// parallelization even slows this step down if many events are used - prob. because each thread requires a distances vector and they don't all fit into a faster cache level (sorting becomes the extreme bottle neck), anyways, it doesn't get faster even without distance vector, so maybe the cache misses when writing
	for (int j = 0; j < gridSize; ++j)
	{
// 		if ((j % 99) == 0)
// 		{
// 			gettimeofday(&endt, 0); cout << j << " ... Runtime of last 100 loops : " << endt.tv_sec - startt.tv_sec << endl; gettimeofday(&startt, 0);
// 		}
// 		cout << j << " / " << grid.size() << endl;

		// evaluate all distances
		vector<double> distances;
		for (unsigned int k = 0; k < eventsSize; ++k)
		{
// 			if (quickDist(grid[j], events[k]) < distanceLimit)	// TODO: test if this actually speeds up the computation!
// 			{
		// TODO: speed up by continue if the conditions are not me// TODO: speed up by continue if the conditions are not mett
				double d = distanceRad(grid[j], events[k]);
				if (d < distanceLimit)
				{
					distances.push_back(d);
	//				aroundRatioAll[j][static_cast<unsigned int>(d/stepSize_)] += addValue;
				}
// 			}
		}
//		for (unsigned int i = 0; i < distances.size(); ++i)
//		{
//			++aroundRatioAll[j][distances[i]/stepSize_];
//		}
//
//		for (unsigned int i = 0; i < numBins_; ++i)
//		{
//			aroundRatioAll[j][i] /= evtSize;
//		}

		// sort them ascending
		std::sort(distances.begin(), distances.end());

		// fill the 180 ratios for this gridpoint accordingly
		vector<double>::const_iterator it = distances.begin();
		unsigned int counter = 0;
		for (unsigned int i = 0; i < numBins_; ++i){
			counter = 0; // resetting the counter for each bin means no cumulative search
			while ( ((*it) < (i+1)*stepSize_) && (it != distances.end()) ){
				++counter;
				++it;
			}
			aroundRatioAll[j][i] = static_cast<double>(counter)/static_cast<double>(evtSize);
		}
// 		cout << "distances size = " << distances.size() << endl;
// 		if ( (j == 10) || (j == 10000) )
// 		{
// 			for (unsigned int i = 0; i < numBins_; ++i)
// 				std::cout << "DEBUG " << i << " " << j << " : " << aroundRatioAll[j][i] << std::endl;
// 		}
	}

	gettimeofday(&endtTotal, 0); cout << "Total runtime of computeAllRatios: " << endtTotal.tv_sec - starttTotal.tv_sec << endl;
// 	cout << "DEBUG allScales::computeAllRatios finished" << endl;
}



void allScales::computePoissonProbabilityForAllRatios(SG::allScales const & normalization, unsigned int eventSize)
{
	timeval startt, endt;
	gettimeofday(&startt, 0);

	double evtSize = static_cast<double>(eventSize);

	//omp_set_num_threads(16);
	#pragma omp parallel for
	for (int j = 0; j < (int)aroundRatioAll.size(); ++j){
		for (unsigned int i = 0; i < numBins_; ++i){
			long double mean = normalization.aroundRatioAll[j][i]*evtSize;
			double found = aroundRatioAll[j][i]*evtSize;

			// Sum the p values for all values lower than the found number of events
			long double pSum = 0.0;
			for (unsigned int x = 0; x < static_cast<unsigned int>(found+0.5); ++x){
				long double thisFound = static_cast<long double>(x);
				pSum += computePValue(mean, thisFound);
			}

			const long double maxPSum = 1.0 - 1.0/100000000000;
			if (pSum > (maxPSum)){
				pSum = maxPSum;
// 				cerr << "No Probability " << pSum << " greater one allowed!" << endl;
			}
			aroundRatioAll[j][i] = 1.0-pSum;
		}
	}

	gettimeofday(&endt, 0); cout << "Runtime of computePoissonProbabilityForAllRatios: " << endt.tv_sec - startt.tv_sec << endl;
}


// load every .norm file in folder folder (that matches size and stepsize) and try to treat it as a valid normalization
int allScales::loadNormalization2(unsigned int gridsize, unsigned int numBins, unsigned int numEvents, std::string folder)
{
	if (folder == "")
	{
		cerr << "No or empty directory specified for data reading.";
		return -1;
	}
	string searchedEnding = ".nres2";
	string searchSize = "size" + boost::lexical_cast<string>(gridsize) + "_";
	string searchNumBins = "bins" + boost::lexical_cast<string>(numBins) + "_";
	string searchNumEvents = "events" + boost::lexical_cast<string>(numEvents) + "_";

	vector<string> inputFiles;
	int len = folder.length();
	if (folder.at(len-1) != '/')
		folder += "/";

	DIR* myDir;
	char* dirName = (char*) folder.c_str();
	struct dirent* entry;
	if (!(myDir = opendir(dirName))) {
		cerr << "Failed to open directory " << dirName << endl;
		return -1;
	}
	
	while ((entry = readdir(myDir)))
	{
		string temp = string(entry->d_name);
		string lowerFileName = temp;	// Assign lowerFileName to temp to have the right size reserved (transform does not allocate memory)
		std::transform(temp.begin(), temp.end(), lowerFileName.begin(), ::tolower);
		
		std::string::size_type found = lowerFileName.rfind(searchedEnding);
		if (found == (lowerFileName.length()-searchedEnding.length())){	// find it at the end
			if (lowerFileName.find(searchSize) != string::npos){	// only with the right gridSize
				if (lowerFileName.find(searchNumBins) != string::npos){	// only with the right stepSize
					if (lowerFileName.find(searchNumEvents) != string::npos){	// only for the right number of events
						inputFiles.push_back( temp );
					}
				}
			}
		}
	}
	closedir(myDir);
	if (inputFiles.empty()){
		cerr << "No " << searchedEnding << " file found in directory " << folder << endl;
		return -4;
	}
	
	int ret = 0;
	ret = getNormalizationResults(folder + *inputFiles.begin(), gridsize, numBins);
	if (ret != 0){
		cout << "Failed loading normalization file " << *inputFiles.begin() << " from folder " << folder << endl;
		throw -5;
	}

	double stepDegree = maximumDistanceToEvaluateInDegree/(numBins-1);
	set(stepDegree);
	return 0;
}

int allScales::loadFileToData(std::string file, unsigned int gridsize, unsigned int numBins, std::string searchStringEnding, std::vector<std::vector<double> >& data){
	// cerr << "DEBUG: loadFileToData starts" << endl;
	if (file == ""){
		cerr << "No file specified for data reading.";
		return -1;
	}
	
	ifstream rein;
	rein.open(file.c_str());
	if (! rein.is_open() ){
		cerr << file << " could not be opened!" << endl;
		return -2;
	}
	
	string lowerFileName = file;	// Assign lowerFileName to temp to have the right size reserved (transform does not allocate memory)
	std::transform(file.begin(), file.end(), lowerFileName.begin(), ::tolower);
	std::string::size_type found2 = lowerFileName.rfind(searchStringEnding);
	if (found2 == std::string::npos){
		cerr << "Dealing with a " << searchStringEnding << " file (" << file << ") that does not contain \"" << searchStringEnding << "\"! This should not have happend!" << endl;
		return -3;
	}
	
	data.clear();

	const int anz = 64000;
	char temp2[anz];
	rein.getline(temp2, anz);
	while (rein.good()){
		string thisLineAsString = string(temp2);
		size_t wo = thisLineAsString.find("#");
		string currentPart = thisLineAsString.substr(0, wo);	// get the line from start to the #

		if (currentPart.length() == 0){
			rein.getline(temp2, anz);
			continue;
		}
		
		stringstream sstr2;
		sstr2 << currentPart;
		if (!sstr2.good()){
			cerr << "Stringstream failed to hold the following data: " << currentPart << endl;
			return -4;
		}
		
		vector<double> thisLine;
		while (sstr2.good()){
			double entry;
			
			string temp = "";
			sstr2 >> temp;
			
			if (temp == "" || temp == " "){
				continue;
			}
			if (temp == "nan" || temp == "-nan"){
				cerr << "WARNING: Read " << temp << " in file " << file << " ! Substituting this mean value with 0.0." << endl;
				entry = 0.0;
			} else {
				entry = boost::lexical_cast<double>(temp);
			}
			
			thisLine.push_back(entry);
		}

		// Event was accepted
		data.push_back(thisLine);
		
		rein.getline(temp2, anz);
	}

	// Here we have read the entire file and parsed it to data

	// Now check for consistency with what was desired
	if (data.size() != gridsize){
		cerr << "The read number of lines (" << data.size() << ") in file " << file << " does not match the desired size of the normalization grid (" << gridsize << "). Dropping this file." << endl;
		return -6;
	}
	for (unsigned int i = 0; i < gridsize; ++i){
		if (data[i].size() != numBins){
			cerr << "addNormalizationResults: The read number of entries (" << data[i].size() << ") in line " << i << " in file " << file << " does not match the desired number of bins " << numBins << " (and therefore angular resolution) of the normalization grid. Dropping this file." << endl;
			return -7;
		}
	}
	
	return 0;
}


int allScales::getNormalizationResults(std::string file, unsigned int gridsize, unsigned int numBins)
{
  	cout << "Normalization adding file " << file << endl;
	
//	std::vector<std::vector<double> > data;
	
	string searchStringEnding = ".nres2";
	int ret = loadFileToData(file, gridsize, numBins, searchStringEnding, aroundRatioAll);
	if (ret){
		cerr << "loadFileToData failed with error code " << ret << endl;
		return ret;
	}
//	aroundRatioAll = data;
	
	return 0;
}


int allScales::dumpRatiosAsNormalization(std::string folder, unsigned int numEvents){
	folder += "/";
	string ofFileName;
	ofstream raus;
	
	static unsigned int numOfNormalizationsCounter = 0;
	++numOfNormalizationsCounter;
	
	unsigned int numBins = 0;
	if ( aroundRatioAll.size() ){
		numBins = aroundRatioAll[0].size();
	}
	
	// Output the results to file
	ofFileName = folder + "normalizationResults_events" + boost::lexical_cast<std::string>(numEvents) + "_Size" + boost::lexical_cast<string>(aroundRatioAll.size()) + "_bins" + boost::lexical_cast<string>(numBins) +  "_nr" + boost::lexical_cast<string>(numOfNormalizationsCounter) + ".nres2";
//	ofFileName = folder + "normalizationResults_Size" + boost::lexical_cast<string>(aroundRatioAll.size()) + "_bins" + boost::lexical_cast<string>(numBins) +  "_nr" + boost::lexical_cast<string>(numOfNormalizationsCounter) + ".nres2";
	cout << "Dumping normalization results to file " << ofFileName << endl;
	raus.open(ofFileName.c_str());
	if (raus.is_open() ){
		for (unsigned int i = 0; i < aroundRatioAll.size(); ++i){
			for (unsigned int j = 0; j < aroundRatioAll[i].size(); ++j)
				raus << aroundRatioAll[i][j] << " ";
			raus << endl;
		}
	} else {
		cerr << "Failed to open file " << ofFileName << endl;
		return -1;
	}
	raus.close();
	return 0;
}

void allScales::addRatios(struct allScales const & b)
{
       if (aroundRatioAll.size() != b.aroundRatioAll.size())
       {
               std::cerr << "The two allScales do not have the same size!" << std::endl;
               throw (-1);
       }
       for (unsigned int i = 0; i < aroundRatioAll.size(); ++i)
       {
               if (aroundRatioAll[i].size() != b.aroundRatioAll[i].size())
               {
                       std::cerr << "The allScales at index " << i << " do not have the same size!" << std::endl;
                       throw (-2);
               }
       }
       for (unsigned int i = 0; i < aroundRatioAll.size(); ++i)
       {
               for (unsigned int n = 0; n < aroundRatioAll[0].size(); ++n)
               {
                       aroundRatioAll[i][n] += b.aroundRatioAll[i][n];
               }
       }
}



}	// end of namespace SG
