#ifndef SGRememberBoostHistogram_H
#define SGRememberBoostHistogram_H

#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>

namespace SG{

struct cluster;

typedef boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::density> > acc;
typedef boost::iterator_range<std::vector<std::pair<double, double> >::iterator > histogram_type;

struct rememberBoostHistogram
{
	std::vector<double> values;
	std::vector<double> indices;

	unsigned int size() const;
	void clear();
	void rememberThisHistogram(histogram_type& hist);
	void operator()(histogram_type& hist);
	rememberBoostHistogram();
	rememberBoostHistogram(histogram_type& hist);
};

void getFitHistogramOwn(struct rememberBoostHistogram const & dataOrg, struct rememberBoostHistogram& fit);

void getFitHistogramOwn(const std::vector<double>& data, struct rememberBoostHistogram& fit);


double computeArea(struct rememberBoostHistogram const & data);

void saveHistogramToFile(struct rememberBoostHistogram const & data, std::string const & fileName);

void normalizeHistogramAreaTo(struct rememberBoostHistogram & data, double targetValue);

unsigned int computeLastIndex(struct rememberBoostHistogram const & data);

void smoothHistogram(struct rememberBoostHistogram & data);

void obtainFitParameters(struct rememberBoostHistogram & data, unsigned int firstIndex, unsigned int lastIndex,
		struct rememberBoostHistogram const & dataRenormalizedBkp,
		double & scalingFactorFitFinal, double & exponent);

void generateActualFit(struct rememberBoostHistogram const & data, double scalingFactor, double exponent,
		struct rememberBoostHistogram& fit, double longerFactor = 5.0, double finerFactor = 10.0);

void mergeOriginalAndFit(struct rememberBoostHistogram const & dataRenormalizedBkp, unsigned int firstIndex,
		unsigned int lastIndex, struct rememberBoostHistogram& fit, double finerFactor = 10.0);

unsigned int getFirstIndex(struct rememberBoostHistogram const & data, unsigned int lastIndex);

void eraseEverythingButTail(struct rememberBoostHistogram & data, unsigned int firstIndex);

void computeAndOutputHistogramsFromClusters(std::string folderName,
		unsigned int n, std::vector<std::vector<double> > & relevanceComparisonVals, std::string searchString2,
		std::vector<std::vector<struct rememberBoostHistogram > >& histsAllSegmentationsAllRelevances);


} // end of namespace SG

#endif 
