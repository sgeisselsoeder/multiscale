/*
 * sgcu.hpp
 * Contains a collection of sometimes useful common code sniplets
 *
 *  Created on: 2016-05-25
 *  Based on utils.hpp, 27.2.2012
 *  Author: sgeisselsoeder
 */

#ifndef SGCOMMONUTILS_H_
#define SGCOMMONUTILS_H_

#include "histogram.hpp"
#include "vectors.hpp"
#include "matrix.hpp"
#include "argumentParsing.hpp"
//#include "argumentParsingTemplates.h"

#include <iomanip>      // std::setprecision
#include <numeric>

namespace sgcu
{

template <typename T>
std::string convertNumberToString(T value, unsigned int numberDigits = 2)
{
	std::stringstream sstream;
	sstream << std::fixed << std::setprecision(numberDigits) << value;
	return sstream.str();
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
T faculty(T number){
	T result = 1;
	for (T i = 1; i <= number; ++i){
		result *= i;
	}
	return result;
}

// Spits out "good" random numbers from limitStart to limitEnd
template <typename T>
T goodRand(T limitStart, T limitEnd){
//	if (limitEnd < limitStart) std::swap(limitStart, limitEnd);
	T delta = fabs(limitEnd-limitStart);
	return (T) limitStart+goodRand(delta);
}

// Spits out "good" random numbers from 0 to limit-1		// Better than rand()%limit
template <typename T>
T goodRand(T limit){
	return (T) (limit * (((double)rand())/(RAND_MAX)) );
}

} // end of namespace sg

#endif /* SGCOMMONUTILS_H_ */
