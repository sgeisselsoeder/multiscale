/*
 * utils.hpp
 * Contains a collection of sometimes useful common code sniplets
 *
 *  Created on: 27.2.2012
 *  Author: sgeisselsoeder
 */

#ifndef SGUTILS_H_
#define SGUTILS_H_

template <typename T>
void gaussianSmoothNoEnds(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = before[0];
	after[numBins-1] = before[numBins-1];
	for (unsigned int i = 1; i < numBins-1; ++i)
		after[i] = (before[i-1] + 2.0*before[i] + before[i+1])/4.0;
}

template <typename T>
void gaussianSmoothNoEnds(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	gaussianSmoothNoEnds(vec2, vec);
}

template <typename T>
void smoothNoEnds(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = before[0];
	after[numBins-1] = before[numBins-1];
	for (unsigned int i = 1; i < numBins-1; ++i){
		after[i] = (before[i-1] + before[i] + before[i+1])/3.0;
	}
}
template <typename T>
void smoothNoEnds(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothNoEnds(vec2, vec);
}
template <typename T>
void smoothNoZeros(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = before[0];
	if ((before[0] != 0.0) && (before[1] != 0.0)) {
		after[0] = (before[0] + before[1])/2.0;
	}
	after[numBins-1] = before[numBins-1];
	if ((before[numBins-1] != 0.0) && (before[numBins-2] != 0.0)) {
		after[numBins-1] = (before[numBins-2] + before[numBins-1])/2.0;
	}
	for (unsigned int i = 1; i < numBins-1; ++i){
		if ((before[i-1] != 0.0) && (before[i] != 0.0) && (before[i+1] != 0.0)) {
			after[i] = (before[i-1] + before[i] + before[i+1])/3.0;
		} else {
			after[i] = before[i];
		}
	}
}
template <typename T>
void smoothNoZeros(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothNoZeros(vec2, vec);
}

template <typename T>
void smoothNoZerosKillingEnds(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = before[0];
	if ((before[0] != 0.0) && (before[1] != 0.0)) {
		after[0] = (before[0] + 1.0*before[1])/2.0;
	}
	after[numBins-1] = before[numBins-1];
	if ((before[numBins-1] != 0.0) && (before[numBins-2] != 0.0)) {
		after[numBins-1] = (1.0*before[numBins-2] + before[numBins-1])/2.0;
	}
	for (unsigned int i = 1; i < numBins-1; ++i){
		if ((before[i-1] != 0.0) && (before[i] != 0.0) && (before[i+1] != 0.0)) {
//			after[i] = (before[i-1] + 2.0*before[i] + before[i+1])/4.0;
			after[i] = (before[i-1] + 4.0*before[i] + before[i+1])/6.0;
		} else {
			after[i] = before[i];
		}
	}
}
template <typename T>
void smoothNoZerosKillingEnds(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothNoZerosKillingEnds(vec2, vec);
}

template <typename T>
void smoothNoCloseToZerosKillingEnds(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}

	const double whatsClose = 6.0/(190.0);
//	std::cout << "whatsClose " << whatsClose << std::endl;
	unsigned int numBins = before.size();
	unsigned int whatsCloseBins = static_cast<unsigned int>(whatsClose*static_cast<double>(numBins)+0.5);
//	std::cout << "whatsCloseBins " << whatsCloseBins << std::endl;

	after.resize(numBins);
	after[0] = before[0];
	if ((before[0] != 0.0) && (before[1] != 0.0)) {
		after[0] = (before[0] + 1.0*before[1])/2.0;
	}
	after[numBins-1] = before[numBins-1];
	if ((before[numBins-1] != 0.0) && (before[numBins-2] != 0.0)) {
		after[numBins-1] = (1.0*before[numBins-2] + before[numBins-1])/2.0;
	}
	for (unsigned int i = 1; i < numBins-1; ++i){
		// check if any value within the desired range is zero
		bool isCloseToZero = false;
		for (unsigned int j = 0; j < whatsCloseBins; ++j)
		{
			// do not check for zero values outside of the bounds
			int whereAmI = static_cast<int>(i)-static_cast<int>(j);
			if (whereAmI < 0) continue;

			if (before[whereAmI] == 0.0)
			{
				isCloseToZero = true;
				break;
			}
			if (i+j > numBins-1) continue;
			if (before[i+j] == 0.0)
			{
				isCloseToZero = true;
				break;
			}
		}

//		std::cout << i << " " << isCloseToZero << std::endl;

		if (isCloseToZero == false) {
			after[i] = (before[i-1] + before[i] + before[i+1])/3.0;
		} else {
			after[i] = before[i];
		}
	}
}
template <typename T>
void smoothNoCloseToZerosKillingEnds(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothNoCloseToZerosKillingEnds(vec2, vec);
}


template <typename T>
void smoothKillingEnds(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = (before[0] + 2.0*before[1])/3.0;
	after[numBins-1] = (2.0*before[numBins-2] + before[numBins-1])/3.0;
	for (unsigned int i = 1; i < numBins-1; ++i)
		after[i] = (before[i-1] + before[i] + before[i+1])/3.0;

// 	after[0] = after[1];
// 	after[numBins-1] = after[numBins-2];
}

template <typename T>
void smoothKillingEnds(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothKillingEnds(vec2, vec);
}

template <typename T>
void smoothKeepingEnds(std::vector<T>& before, std::vector<T>& after){
	if (before.size() < 2) {
		after = before;
		return;
	}
	unsigned int numBins = before.size();
	after.resize(numBins);
	after[0] = (2.0*before[0] + 1.0*before[1])/3.0;
	after[numBins-1] = (1.0*before[numBins-2] + 2.0*before[numBins-1])/3.0;
	for (unsigned int i = 1; i < numBins-1; ++i)
		after[i] = (before[i-1] + before[i] + before[i+1])/3.0;

// 	after[0] = after[1];
// 	after[numBins-1] = after[numBins-2];
}

template <typename T>
void smoothKeepingEnds(std::vector<T>& vec){
	std::vector<T> vec2 = vec;
	smoothKeepingEnds(vec2, vec);
}

#endif /* SGUTILS_H_ */
