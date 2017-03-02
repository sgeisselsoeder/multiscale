/*
 * sgcu.hpp
 * Contains a collection of sometimes useful common code sniplets
 *
 *  Created on: 2016-05-25
 *  Based on utils.hpp, 27.2.2012
 *  Author: sgeisselsoeder
 */

#ifndef SGCOMMONUTILS_VECTORS_H_
#define SGCOMMONUTILS_VECTORS_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>    // std::sort

namespace sgcu
{

template <typename T>
std::ostream& operator<<(std::ostream& raus, std::vector<T> const & vec)
{
	for (unsigned int i = 0; i < vec.size(); ++i)
	{
		raus << vec[i] << std::endl;
	}
	return raus;
}

template <typename T>
void saveVectorToFile(std::vector<T> const & data, std::string fileName)
{
	std::ofstream raus;
	raus.open(fileName.c_str());
	if ( (! raus.is_open()) || (! raus.good()) )
	{
		std::cerr << "Failed to open file " << fileName << " for vector output." << std::endl;
		throw (-1);
	}
	raus << data;
	if (! raus.good() )
	{
		std::cerr << "Failed to output vector to file " << fileName << std::endl;
		throw(-2);
	}
	raus.close();
}

template <typename T>
void loadVectorFromFile(std::string fileName, std::vector<T> & data)
{
	std::ifstream rein(fileName.c_str());
	if (! rein.good())
	{
		std::cerr << "Failed to open file " << fileName << std::endl;
		throw -1;
	}
	data.clear();
	while (rein.good())
	{
		T val;
		rein >> val;
		if (rein.good())
		{
			data.push_back(val);
		}
	}
	std::cout << "Loaded " << data.size() << " entries from file " << fileName << std::endl;
	rein.close();
}

template <typename T>
unsigned int maxIndex(std::vector<T> const & vec)
{
	return std::max_element(vec.begin(), vec.end())-vec.begin();
}

template <typename T>
unsigned int minIndex(std::vector<T> const & vec)
{
	return std::min_element(vec.begin(), vec.end())-vec.begin();
}

template <typename T>
T sum(std::vector<T> const & vec)
{
	return std::accumulate( vec.begin(), vec.end(), T(0.0) );
}

template <typename T>
T maxValue(std::vector<T> const & vec)
{
	return *std::max_element(vec.begin(), vec.end());
}

template <typename T>
T minValue(std::vector<T> const & vec)
{
	return *std::min_element(vec.begin(), vec.end());
}

template <typename T>
typename std::vector<T>::const_iterator max_element(std::vector<T> const & vec)
{
	return std::max_element(vec.begin(), vec.end());
}

template <typename T>
typename std::vector<T>::iterator max_element(std::vector<T> & vec)
{
	return std::max_element(vec.begin(), vec.end());
}

template <typename T>
typename std::vector<T>::const_iterator min_element(std::vector<T> const & vec)
{
	return std::min_element(vec.begin(), vec.end());
}

template <typename T>
typename std::vector<T>::iterator min_element(std::vector<T> & vec)
{
	return std::min_element(vec.begin(), vec.end());
}

template <typename T>
void derive(std::vector<T>& before, std::vector<T>& after){
// 	std::cerr << "derive start" << std::endl;
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins-1);
	for (unsigned int i = 0; i < numBins-1; ++i)
		after[i] = before[i+1]-before[i];
// 	std::cerr << "derive end" << std::endl;
}
template <typename T>
void derive(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	derive(vec2, vec);
}

template <typename T>
void normalize(std::vector<T>& before, std::vector<T>& after)
{
//	std::cerr << "normalize start" << std::endl;
	T max = *sgcu::max_element(before);
	T min = *sgcu::min_element(before);
	T diff = max - min;
	if (diff != 0.0)
	{	// normalize
		for (unsigned int i = 0; i < before.size(); ++i)
			after[i] = (before[i]-min)/diff;
	}
	else
	{	// if all values have been the same before, they all are 1.0 after
		for (unsigned int i = 0; i < before.size(); ++i)
					after[i] = 1.0;
	}
//	std::cerr << "normalize end" << std::endl;
}
template <typename T>
void normalize(std::vector<T>& vec)
{
	std::vector<T> vec2 = vec;
	normalize(vec2, vec);
}

template <typename T>
void gaussianSmooth(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = (before[0] + before[1])/2.0;
	after[numBins-1] = (before[numBins-2] + before[numBins-1])/2.0;
	for (unsigned int i = 1; i < numBins-1; ++i)
		after[i] = (before[i-1] + 2.0*before[i] + before[i+1])/4.0;
}

template <typename T>
void gaussianSmooth(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	gaussianSmooth(vec2, vec);
}

template <typename T>
void smoothNoEdge(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = before[0];
	after[numBins-1] = before[numBins-1];
	for (unsigned int i = 1; i < numBins-1; ++i)
		after[i] = (before[i-1] + before[i] + before[i+1])/3.0;
}

template <typename T>
void smoothNoEdge(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothNoEdge(vec2, vec);
}

template <typename T>
void smooth(std::vector<T>& before, std::vector<T>& after){
// 	std::cerr << "smooth2 start" << std::endl;
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = (before[0] + before[1])/2.0;
	after[numBins-1] = (before[numBins-2] + before[numBins-1])/2.0;
	for (unsigned int i = 1; i < numBins-1; ++i)
		after[i] = (before[i-1] + before[i] + before[i+1])/3.0;
// 	std::cerr << "smooth2 end" << std::endl;
}

template <typename T>
void smooth(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smooth(vec2, vec);
}

template <typename T>
T median(std::vector<T>& vals){
	if (vals.size() == 0) return 0.0;
	if (vals.size() == 1) return vals[0];

	std::sort(vals.begin(), vals.end());
	return vals[vals.size()*0.5];
}

template <typename T>
void smoothMedian5(std::vector<T>& before, std::vector<T>& after){
	if (before.size() <= 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = (before[0] + 2.0*before[1])/3.0;
	after[numBins-1] = (2.0*before[numBins-2] + 1.0*before[numBins-1])/3.0;

	for (unsigned int i = 2; i < numBins-2; ++i){
		std::vector<T> vals;
		vals.push_back(before[i]);
		vals.push_back(before[i-1]);
		vals.push_back(before[i+1]);
		vals.push_back(before[i-2]);
		vals.push_back(before[i+2]);
		after[i] = median(vals);
	}
	{
		unsigned int i = 1;
		std::vector<T> vals;
		vals.push_back(before[i]);
		vals.push_back(before[i-1]);
		vals.push_back(before[i+1]);
		after[i] = median(vals);
	}
	{
		unsigned int i = numBins-2;
		std::vector<T> vals;
		vals.push_back(before[i]);
		vals.push_back(before[i-1]);
		vals.push_back(before[i+1]);
		after[i] = median(vals);
	}
}

template <typename T>
void smoothMedian5(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothMedian5(vec2, vec);
}

template <typename T>
void smoothMedian3(std::vector<T>& before, std::vector<T>& after){
	if (before.size() <= 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);

	for (unsigned int i = 1; i < numBins-1; ++i){
		std::vector<T> vals;
		vals.push_back(before[i]);
		vals.push_back(before[i-1]);
		vals.push_back(before[i+1]);
		after[i] = median(vals);
	}
	after[0] = before[0];
	after[numBins-1] = before[numBins-1];
}

template <typename T>
void smoothMedian3(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothMedian3(vec2, vec);
}

template <typename T>
void smoothMedian(std::vector<T>& vec)
{
	smoothMedian3(vec);
}

template <typename T>
T accumulate(typename std::vector<T> const & vec)
{
	return std::accumulate( vec.begin(), vec.end(), T(0.0) );
}

//// sum all entries in the vector including the one at starting index, excluding the value at ending index
//template <typename T>
//T sum(typename std::vector<T>& vec, unsigned int starting = 0, unsigned int ending = 0){
//	if (starting > vec.size()){
//		if (vec.size() > 0){
//			std::cout << "Trying to sum entries after the end of the vector. Stopping after the last entry." << std::endl;
//		}
//		return 0.0;
//	}
//	if (ending == 0){
//		ending = vec.size();
//	}
//	if (ending > vec.size()){
//		std::cout << "Trying to sum entries after the end of the vector. Stopping after the last entry." << std::endl;
//		ending = vec.size();
//	}
//	T result = 0.0;
//	for (unsigned int i = starting; i < ending; ++i)
//	{
//		result += vec[i];
//	}
//	return result;
//}

template <typename T>
void meanAndStdDev(std::vector<T> const & data, double& mean, double& stdDev)
{
	mean = 0.0;
	stdDev = 0.0;
	if (data.size() == 0)
	{
		return;
	}
	double temp = 0.0;
	for (unsigned int i = 0; i < data.size(); ++i){
		double val = data[i]; 
		mean += val;
		temp += val*val;
	}
	mean /= data.size();
	temp /= data.size();
	stdDev = std::sqrt(std::fabs(temp - mean*mean));
}

} // end of namespace sg

#endif /* SGCOMMONUTILS_VECTORS_H_ */
