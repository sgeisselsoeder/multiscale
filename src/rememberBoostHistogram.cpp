#include "rememberBoostHistogram.h"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include "sgcu/sgcu.hpp"
#include "utils.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace SG{

//using namespace boost::accumulators;

// build a histogram for a data vector
template <typename T>
void getHistogram(std::vector<T> const & data, histogram_type & hist)
{
	// create a boost accumulator
	acc myAccumulator( boost::accumulators::tag::density::num_bins = data.size()/10.0, boost::accumulators::tag::density::cache_size = data.size()-2);
	// fill accumulator
	for (unsigned int j = 0; j < data.size(); ++j){
		myAccumulator(data[j]);
	}
	// do the boost histogram
	hist = boost::accumulators::density(myAccumulator);
}

unsigned int rememberBoostHistogram::size() const
{
	unsigned int thisSize = values.size();
	if (indices.size() < thisSize)
		thisSize = indices.size();
	return thisSize;
}

void rememberBoostHistogram::clear()
{
	values.clear();
	indices.clear();
}

void rememberBoostHistogram::rememberThisHistogram(histogram_type& hist)
{
	clear();
	for (unsigned int i = 0; i < (unsigned int) hist.size(); ++i)
	{
		indices.push_back(hist[i].first);
		values.push_back(hist[i].second);
	}
}

void rememberBoostHistogram::operator()(histogram_type& hist)
{
	rememberThisHistogram(hist);
}

rememberBoostHistogram::rememberBoostHistogram()
{
}

rememberBoostHistogram::rememberBoostHistogram(histogram_type& hist)
{
	rememberThisHistogram(hist);
}

double computeArea(struct rememberBoostHistogram const & data)
{
	double area = 0.0;
	for (unsigned int i = 0; i < data.size(); ++i)
	{
		area += data.values[i];
	}
	return area;
}

void getFitHistogramOwn(struct rememberBoostHistogram const & dataOrg, struct rememberBoostHistogram& fit)
{
	// NOTE: this implicitly assumes all histogram values to be positive. true for a histograms, not necessarily true for all inputs!
	if ( computeArea(dataOrg) == 0.0 )
	{
		fit = dataOrg;
		return;
	}

	const bool doDebugOutput = false;

	struct rememberBoostHistogram data = dataOrg;
	if (doDebugOutput) saveHistogramToFile(data, "comparisonValsNormalization3.txt");

	// make sure the histogram is normalized to 1.0
	normalizeHistogramAreaTo(data, 1.0);

	if (doDebugOutput) saveHistogramToFile(data, "comparisonValsNormalization3Renorm.txt");

	unsigned int lastIndex = computeLastIndex(data);

	//	cerr << "new: setting last entry of histogram to zero" << endl;
	data.values[data.size()-1] = 0.0;

	struct rememberBoostHistogram dataBkp = data;
	struct rememberBoostHistogram dataRenormalizedBkp = data;

	smoothHistogram(data);

	if (doDebugOutput) saveHistogramToFile(data, "comparisonValsNormalization3Smoothed.txt");

	unsigned int firstIndex = getFirstIndex(data, lastIndex);

	// do not consider values from before the tail to fit the tail!
	eraseEverythingButTail(data, firstIndex);
	eraseEverythingButTail(dataBkp, firstIndex);

	double areaBefore = computeArea(dataBkp);
	normalizeHistogramAreaTo(data, areaBefore);

	if (doDebugOutput) saveHistogramToFile(data, "comparisonValsNormalization3SmoothedRenorm.txt");

	double scalingFactor = 0.0;
	double exponent = 0.0;
	obtainFitParameters(data, firstIndex, lastIndex, dataRenormalizedBkp, scalingFactor, exponent);

	generateActualFit(data, scalingFactor, exponent, fit);

	mergeOriginalAndFit(dataRenormalizedBkp, firstIndex, lastIndex, fit);

	if (doDebugOutput) saveHistogramToFile(fit, "comparisonValsNormalization3Fitted.txt");
}

void getFitHistogramOwn(const std::vector<double>& data, struct rememberBoostHistogram& fit)
{
	// do a histogram
	acc myAccumulator( boost::accumulators::tag::density::num_bins = 100, boost::accumulators::tag::density::cache_size = data.size());
	for (unsigned int i = 0; i < data.size(); ++i)
	{
		myAccumulator(data[i]);
	}
	histogram_type hist = boost::accumulators::density(myAccumulator);
	struct rememberBoostHistogram histNew(hist);

	getFitHistogramOwn(histNew, fit);
}


void outputHistogram(histogram_type const & hist, std::string fileName)
{
	std::ofstream raus;
	raus.open(fileName.c_str());
	if ( (! raus.is_open()) || (! raus.good()) )
	{
		std::cerr << "Failed to open file " << fileName << " for histogram output." << endl;
		throw -1;
	}

	for (unsigned int i = 0; i < (unsigned int) hist.size(); ++i)
	{
		raus << hist[i].first << " " << hist[i].second << endl;
	}

	if (! raus.good() )
	{
		std::cerr << "Failed to write histogram to file " << fileName << endl;
		throw -1;
	}
	raus.close();

}

void outputRememberHistogram(rememberBoostHistogram const & fit, std::string fileName)
{
	std::ofstream raus;
	raus.open(fileName.c_str());
	if ( (! raus.is_open()) || (! raus.good()) )
	{
		std::cerr << "Failed to open file " << fileName << " for rememberBoostHistogram output." << endl;
		throw -1;
	}

	for (unsigned int i = 0; i < fit.size(); ++i)
	{
		raus << fit.indices[i] << " " << fit.values[i] << endl;
	}
	raus.close();

	if (! raus.good() )
	{
		std::cerr << "Failed to write rememberBoostHistogram to file " << fileName << endl;
		throw -1;
	}
	raus.close();
}

// TODO: check if both outputs are actually used!
// for a fixed segmentation threshold: compute histograms for each metric from the clusters found in all pseudo experiments
void computeAndOutputHistogramsFromClusters(std::string folderName, unsigned int n, std::vector<std::vector<double> > & relevanceComparisonVals,
		std::string searchStringSegmentationThreshold,
		std::vector<std::vector<struct rememberBoostHistogram > >& histsAllSegmentationsAllRelevances)
{
		histsAllSegmentationsAllRelevances[n].resize(relevanceComparisonVals.size());
		for (unsigned int j = 0; j < relevanceComparisonVals.size(); ++j)	// for all metrics
		{
			// create a histogram from all clusters that match this segmentation and metric
			acc myAcc(boost::accumulators::tag::density::num_bins = 100, boost::accumulators::tag::density::cache_size = relevanceComparisonVals[j].size());
			for (unsigned int i = 0; i < relevanceComparisonVals[j].size(); ++i)
			{
				myAcc(relevanceComparisonVals[j][i]);
			}
			histogram_type hist = boost::accumulators::density(myAcc);
			outputHistogram(hist, folderName + "/hist_Size" + searchStringSegmentationThreshold + "_Metric" + boost::lexical_cast<std::string>(j) + ".txt");

			struct rememberBoostHistogram data(hist);
			struct rememberBoostHistogram fit;
			// compute a fit that approximates the tail of the observed data
			getFitHistogramOwn(data, fit);
			outputRememberHistogram(fit, folderName + "/hist_Size" + searchStringSegmentationThreshold + "_Metric" + boost::lexical_cast<std::string>(j) + "fit.txt");

			histsAllSegmentationsAllRelevances[n][j] = fit;
		}
}















































void saveHistogramToFile(struct rememberBoostHistogram const & data, std::string const & fileName)
{
	std::ofstream raus(fileName.c_str());
	if (! raus.good())
	{
		std::cerr << "Failed to open file " << fileName << " for output of histogram data" << std::endl;
		throw (-1);
	}
	for (unsigned int i = 0; i < data.size(); ++i)
	{
		raus << data.indices[i] << " " << data.values[i] << endl;
		if (! raus.good())
		{
			std::cerr << "Failed to write histogram data entry " << i << " of " << data.size() << " to file " << fileName << std::endl;
			throw (-2);
		}
	}
	if (! raus.good())
	{
		std::cerr << "Failed to write last of " << data.size() << " histogram data entries to file " << fileName << std::endl;
		throw (-2);
	}
	raus.close();
}



void normalizeHistogramAreaTo(struct rememberBoostHistogram & data, double targetValue)
{
	double area = computeArea(data);
	double rescaleFactor = 0.0;
	if (area != 0.0)
	{
		rescaleFactor = targetValue/area;
		for (unsigned int i = 0; i < data.size(); ++i)
		{
			data.values[i] *= rescaleFactor;
		}
	}
}

unsigned int computeLastIndex(struct rememberBoostHistogram const & data)
{
	unsigned int lastIndex = 0;
	unsigned int notZeroIndex = 0;
	unsigned int notZeroIndexBefore = 0;
	for (unsigned int i = 0; i < data.size(); ++i)
	{
		if (data.values[i] > 0.0)
		{
			lastIndex = notZeroIndexBefore;
			notZeroIndexBefore = notZeroIndex;
			notZeroIndex = i;
		}
	}
	if (lastIndex == 0)
	{
		lastIndex = notZeroIndexBefore;
	}
	if (lastIndex == 0)
	{
		lastIndex = notZeroIndex;
	}
	if (lastIndex == 0)
	{
		lastIndex = data.size();
	}
	return lastIndex;
}

void smoothHistogram(struct rememberBoostHistogram & data)
{
	// NOTE: I just noted that the number of smoothings correlates with the size of the vector.
	// This has been intentional, but it hasn't been checked if this scales to less or more entries.
	for (unsigned int i = 0; i < 0.75*data.size(); ++i)
		smoothNoEnds(data.values);
	for (unsigned int i = 0; i < 0.25*data.size(); ++i)
		sgcu::smooth(data.values);
}



double evaluateFitL2(const std::vector<double>& data, const std::vector<double>& fit){
	if (data.size() != fit.size()){
		cerr << "evaluateFitL2 cannot compare data and fit of different size, data: " << data.size() << " and fit: " << fit.size() << endl;
		return std::numeric_limits<double>::max();
	}

	double diffSum = 0.0;
	for (unsigned int i = 0; i < data.size(); ++i){
		double diff = data[i] - fit[i];
		diffSum += (diff*diff);
	}
	return diffSum;
}

void obtainFitParameters(struct rememberBoostHistogram & data, unsigned int firstIndex, unsigned int lastIndex,
		struct rememberBoostHistogram const & dataRenormalizedBkp,
		double & scalingFactorFitFinal, double & exponent)
{
	// evaluate when each value of the smoothed right flanc is below half of its value - store the distance in x
	//cout << "maxIndex = " << maxIndex << endl;
	//cout << "firstIndex = " << firstIndex << endl;
	//cout << "lastIndex = " << lastIndex << endl;
	unsigned int num = 0;
	double howLongUntilBelowHalf = 0.0;
	for (unsigned int i = firstIndex; i < lastIndex; ++i){
		double currentY100 = data.values[i];
		if (currentY100 <= 0.0)
			break;
		double currentX = data.indices[i];
		for (unsigned int j = i+1; j < data.size(); ++j){
	//		for (unsigned int j = i+1; j < lastIndex; ++j){	// this was the default
			double currentY = data.values[j];
			if (currentY <= 0.5*currentY100){
				if (currentY > 0.0){
					double distance = data.indices[j] - currentX;
					//cout << i << " " << j << " distance: " << distance << "    histjfirst = " << data.indices[j] << " currentX = " << currentX << " currentY100 = " << currentY100 << " currentY = " << currentY << endl;
					howLongUntilBelowHalf += distance;
					++num;
				}
				break;
			}
		}
	}
	//cout << "nowLongUntil = " << howLongUntilBelowHalf << endl;
	//cout << "num = " << num << endl;
	if (num != 0){
		howLongUntilBelowHalf /= static_cast<double>(num);
	}

	double stepSize = data.indices[1] - data.indices[0];

	//cout << "Mean x distance until y value is below half of its starting value: " << howLongUntilBelowHalf << endl;
	howLongUntilBelowHalf  -= 0.5*stepSize;
	//cout << "x distance corrected for stepsize: " << howLongUntilBelowHalf << endl;
	//double
	exponent = -1.0/howLongUntilBelowHalf;
	//cout << "exponent of fit: " << exponent << endl;

	double minDistance = std::numeric_limits<double>::max();
	//double
	scalingFactorFitFinal = 0.0;

	std::vector<double> histValsFitComparison = dataRenormalizedBkp.values;
	for (unsigned int j = 0; j < firstIndex; ++j){
		histValsFitComparison[j] = 0.0;
	}
	for (unsigned int j = lastIndex; j < histValsFitComparison.size(); ++j){
		histValsFitComparison[j] = 0.0;
	}

	std::vector<double> fitVals = histValsFitComparison;
	for (unsigned int i = firstIndex; i < lastIndex; ++i){
		double startingXValue = data.indices[i];
		double expectedStartingValue = data.values[i];
		double plainStartingValue = pow(2.0, exponent*startingXValue);
		double scalingFactorFit = expectedStartingValue/plainStartingValue;

		for (unsigned int j = firstIndex; j < lastIndex; ++j){
			double currentX = data.indices[j];
			fitVals[j] = scalingFactorFit*pow(2.0, exponent*currentX);
		}

		double distance = evaluateFitL2(histValsFitComparison, fitVals);
		//cout << "index " << i << " distance " << distance << endl;
		if (distance < minDistance){
			minDistance = distance;
			scalingFactorFitFinal = scalingFactorFit;
		}
	}

	//cout << "scalingfactor: " << scalingFactorFitFinal << endl;
}



void generateActualFit(struct rememberBoostHistogram const & data, double scalingFactor, double exponent,
		struct rememberBoostHistogram& fit, double longerFactor, double finerFactor)
{
	double stepSize = data.indices[1] - data.indices[0];
	unsigned int targetSize = data.size()*longerFactor*finerFactor;
	fit.clear();
	fit.indices.resize(targetSize);
	fit.values.resize(targetSize);
	double startIndex = data.indices[0];
	for (unsigned int i = 0; i < targetSize; ++i){
		fit.indices[i] = startIndex+stepSize*(static_cast<double>(i)/finerFactor);
	}
	for (unsigned int i = 0; i < targetSize; ++i){
		fit.values[i] = scalingFactor*pow(2.0, exponent*fit.indices[i]);
	}
}


void mergeOriginalAndFit(struct rememberBoostHistogram const & dataRenormalizedBkp, unsigned int firstIndex,
		unsigned int lastIndex, struct rememberBoostHistogram& fit, double finerFactor)
{
	unsigned int firstIndexFiner = static_cast<unsigned int>(static_cast<double>(firstIndex)*finerFactor+0.5);
	unsigned int lastIndexFiner = static_cast<unsigned int>(static_cast<double>(lastIndex)*finerFactor+0.5);

	// only start using the fit when you don't trust the actual distribution anymore:
//	double shiftBack = 0.8;	// default
	double shiftBack = 0.6;

	unsigned int startReallyUsingIt = ((1.0-shiftBack)*firstIndex + shiftBack*lastIndex);
	unsigned int startReallyUsingItFiner = ((1.0-shiftBack)*firstIndexFiner + shiftBack*lastIndexFiner)-0.55*finerFactor;

	double area2 = 0.0;
	for (unsigned int i = startReallyUsingIt; i < dataRenormalizedBkp.size(); ++i){
		area2 += dataRenormalizedBkp.values[i];
	}

	double areaFit = 0.0;
	for (unsigned int i = startReallyUsingItFiner; i < fit.values.size(); ++i){
		areaFit += fit.values[i];
	}
	double rescalingFactor2 = 1.0/finerFactor;
	if (areaFit != 0.0)
		rescalingFactor2 = area2/areaFit;

	for (unsigned int i = 0; i < startReallyUsingItFiner; ++i){
		fit.values[i] = 0.0;
	}
	for (unsigned int i = startReallyUsingItFiner; i < fit.size(); ++i){
		fit.values[i] *= rescalingFactor2;
	}
	for (unsigned int i = 0; i < startReallyUsingIt; ++i){
		fit.values[static_cast<unsigned int>(static_cast<double>(i)*finerFactor+0.5)] = dataRenormalizedBkp.values[i];
	}

	// make sure the histogram is normalized to 1.0
	normalizeHistogramAreaTo(fit, 1.0);
}

unsigned int getFirstIndex(struct rememberBoostHistogram const & data, unsigned int lastIndex)
{
	unsigned int maxIndex = (std::max_element(data.values.begin(), data.values.end())-data.values.begin());
	return static_cast<unsigned int>((2.0*maxIndex+1.0*lastIndex)/3.0);
}

void eraseEverythingButTail(struct rememberBoostHistogram & data, unsigned int firstIndex)
{
	for (unsigned int i = 0; i < firstIndex; ++i)
	{
		data.values[i] = 0.0;
	}
}














} // end of namespace SG
