/**
 * class: sphericalGrid.hpp
 * providing a consistent interface to the gridpoints within the regular spherical grid
 *
 * sgeisselsoeder
 *
 * 28.Nov 2013
 *
 * (c) 2013 Antares Collaboration
 */
	
#ifndef SGSPHERICALGRID_H
#define SGSPHERICALGRID_H

#include <vector>
#include <string>

#include <time.h>
#include <sys/time.h>

#include "candidate.hpp"
#include <iostream>
#include <fstream>
#include "sgcu/sgcu.hpp"

namespace SG{

/**
 * @brief 
 */
template <typename T>
class sphericalGrid
{
public:
	/**
	 * Builds an instance of this class
	 */
	sphericalGrid(){
		fine_ = 0;
	}
  
	/**
	 * Destroys an instance of this class
	 */
	~sphericalGrid(){
		// NOTHING HERE
	}

	void clear(){
		gridpoints_.clear();
		fine_ = 0;
	}
	
	void reset()
	{
		for (unsigned int i = 0; i < size(); ++i){
			gridpoints_[i].multiPurposeValue_ = -1.0;
			gridpoints_[i].proximity_ = 0.0;
		}
	}

	std::vector<SG::candidate<T> > gridpoints_;
	
	// set up the grid to search on
	// fine determines how fine to search ( in 10/fine degrees )
	void setupGrid(unsigned int fine) {
		fine_ = fine;
		gridpoints_.clear();
		for (unsigned int j = 0; j < fine*18+1; ++j){
			double tempZenith = static_cast<double>(M_PI*j/(18*fine)) - 0.5*M_PI;
			// the density factor is necessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			unsigned int numAz = static_cast<unsigned int>(fine*36*densityFactor + 0.5);
			if (numAz == 0) numAz = 1;
			for (unsigned int i = 0; i < numAz; ++i){
				gridpoints_.push_back(SG::candidate<T>(static_cast<double>(2.0*M_PI*i/numAz - M_PI), tempZenith));
			}
		}
	}
	
	double stepsizeInDegrees()
	{
		return 10.0/static_cast<double>(fine);
	}

	void assertNoNansAndInfs() const
	{
		unsigned int counter = 0;
		// cross check if there is a nan or inf value somewhere in the grid
		for (unsigned int i = 0; i < size(); ++i)
		{
			if (std::isnan(gridpoints_[i].proximity_) || std::isinf(gridpoints_[i].proximity_))
			{
				++counter;
			}
		}
		if (counter > 0)
		{
			std::cerr << counter << " gridpoints have a nan or inf value! This is not how it works!" << std::endl;
			throw -1;
		}
	}

	unsigned int countNumberOfNonZeroValues() const
	{
		unsigned int counter = 0;
		for (unsigned int i = 0; i < size(); ++i)
		{
			if (gridpoints_[i].proximity_ != 0.0)
				++counter;
		}
		return counter;
	}

	void correctEquatorialRightAscensionRange(){
                for (unsigned int i = 0; i < gridpoints_.size(); ++i){
                        double ra = gridpoints_[i].azimuth_;
                        if (ra < 0.0){
                                ra = 2.0*M_PI+ra;
                        }
                        gridpoints_[i].azimuth_ = ra;
                }
        }

	void invertEquatorialRightAscensionRange(){
                for (unsigned int i = 0; i < gridpoints_.size(); ++i){
                        gridpoints_[i].invertEquatorialRightAscensionRange();
                }
        }

	void correctEquatorialRightAscensionRangeFranzoesischInvers(std::vector<SG::candidate<T> >& points){
		for (unsigned int i = 0; i < points.size(); ++i){
			points[i].correctEquatorialRightAscensionRangeFranzoesischInvers();
		}
	}

	void correctEquatorialRightAscensionRangeFranzoesischInvers(){
		correctEquatorialRightAscensionRangeFranzoesischInvers(gridpoints_);
	}
	void correctEquatorialRightAscensionRangeFranzoesisch(std::vector<SG::candidate<T> >& points){
		for (unsigned int i = 0; i < points.size(); ++i){
			points[i].correctEquatorialRightAscensionRangeFranzoesisch();
		}
	}

	void correctEquatorialRightAscensionRangeFranzoesisch(){
		correctEquatorialRightAscensionRangeFranzoesisch(gridpoints_);
	}

	void correctEquatorialRightAscensionRangeMinusPi(std::vector<SG::candidate<T> >& points){
		for (unsigned int i = 0; i < points.size(); ++i){
			points[i].correctEquatorialRightAscensionRangeMinusPi();
		}
	}
	void correctEquatorialRightAscensionRangeMinusPi(){
		correctEquatorialRightAscensionRangeMinusPi(gridpoints_);
	}
	
	void correctEquatorialRightAscensionRangePlusPi(std::vector<SG::candidate<T> >& points){
		for (unsigned int i = 0; i < points.size(); ++i){
			points[i].azimuth_ += M_PI;
		}
	}

	void correctEquatorialRightAscensionRangePlusPi(){
		correctEquatorialRightAscensionRangePlusPi(gridpoints_);
	}

	void equatorialToGalactic(std::vector<SG::candidate<T> >& points){
		//assertRangeZeroTo2Pi(points);
		for (unsigned int i = 0; i < points.size(); ++i){
			points[i].equatorialToGalactic();
		}
		//assertRangeMpiToPi(points);
	}

	void equatorialToGalactic(){
		equatorialToGalactic(gridpoints_);	
	}

	void galacticToEquatorial(){
		//assertRangeMpiToPi(gridpoints_);
		for (unsigned int i = 0; i < gridpoints_.size(); ++i){
			gridpoints_.galacticToEquatorial();
		}
		//assertRangeZeroTo2Pi(gridpoints_);
	}
	void setTimeToDecAndEnergyToRa(){
		for (unsigned int i = 0; i < gridpoints_.size(); ++i){
			double oldRa = gridpoints_[i].azimuth_;
			double oldDec = gridpoints_[i].zenith_;

			gridpoints_[i].energy_ = oldRa;
			gridpoints_[i].time_ = oldDec;
		}
	}



	void setFineTEMPORARY(unsigned int f)
	{
		fine_ = f;
	}
	
	unsigned int fine() const
	{
		return fine_;
	}
	
	inline typename std::vector<SG::candidate<T> >::const_iterator begin() const
	{
		return gridpoints_.begin();
	}
	
	inline typename std::vector<SG::candidate<T> >::iterator begin()
	{
		return gridpoints_.begin();
	}
	
	inline typename std::vector<SG::candidate<T> >::const_iterator end() const
	{
		return gridpoints_.end();
	}
	
	inline typename std::vector<SG::candidate<T> >::iterator end()
	{
		return gridpoints_.end();
	}
	
	// Search for a connected path from cluster(center)A to cluster(center)B through the gridpoints
	// consideres only gridpoints with multiPurposeValue_ == value
	// (== within the width of a cluster or in overlapping, not independent parts)
	// return 0 if no path is found, 1 if a path was found
	// should be a const function, but not now ...
	int findPath(const SG::candidate<T>& clusterA, const SG::candidate<T>& clusterB, T value = 1.0) {//const{
		std::vector<SG::candidate<T> > gridpointsBkp = gridpoints_;
		T checked = value+1.0;
		
		// Find the gridpoints corresponding to the two clusters 
		typename std::vector<SG::candidate<T> >::iterator clustA = findGridpoint(clusterA.azimuth_, clusterA.zenith_);
		typename std::vector<SG::candidate<T> >::iterator clustB = findGridpoint(clusterB.azimuth_, clusterB.zenith_);

		// start at cluster A:
		typename std::vector< typename std::vector<SG::candidate<T> >::iterator > searchpoints;
		clustA->multiPurposeValue_ = checked;
		searchpoints.push_back(clustA);
		
		bool foundConnection = false;
		// search as long as we still have unchecked gridpoints
		while (searchpoints.size()){
			// check one gridpoint in each loop
			typename std::vector<SG::candidate<T> >::iterator currentPoint = *searchpoints.begin();

			// if the current point has the same position than the second cluster we have found a connection
			if (currentPoint->samePos(*clustB)){
				foundConnection = true;
				break;
			} else {
				// if this point isn't at the same position than the second cluster, maybe his neighbours are!
				typename std::vector<SG::candidate<T> >::iterator itN,itS,itW,itE; // four are enough, the other ones 
				itN = n(currentPoint);
				itS = s(currentPoint);
				itW = w(currentPoint);
				itE = e(currentPoint);
				
				// now check neighbour points if they haven't been tested already
				if (itW->multiPurposeValue_ == value){
					itW->multiPurposeValue_ = checked;
					searchpoints.push_back(itW);
				}
				if (itE->multiPurposeValue_ == value){
					itE->multiPurposeValue_ = checked;
					searchpoints.push_back(itE);
				}
				if (itN->multiPurposeValue_ == value){
					itN->multiPurposeValue_ = checked;
					searchpoints.push_back(itN);
				}
				if (itS->multiPurposeValue_ == value){
					itS->multiPurposeValue_ = checked;
					searchpoints.push_back(itS);
				}
				
				// remove the currently searched gridpoint from the search vector
				searchpoints.erase(searchpoints.begin());
// 				searchpoints.erase(currentPoint);
			}
		}
		
// 		std::cerr << "fertig    " << i << std::endl;
		
		// if we are here, either a connection has been found and the loop was broken
		// or all gridpoints reachable from clusterA with multiPurposeValue_==1 have been checked 
		// and their multiPurposeValue_ is set to 2
		// In both cases: we are done! :)

		// Clean up: reset the original state of the gridpoints
		gridpoints_ = gridpointsBkp;

		// Done : )
		return foundConnection;
	}
	
	// Set the multipurpose variable to the value segmentNumber for all gridpoints 
	// with a proximity > 0 (= surviving the segmentation as foreground) 
	// and connected to gridIt by other foreground segmented gridpoints
	// returns the number of gridpoints marked within this segment
	int connectSegment(const typename std::vector<SG::candidate<T> >::iterator gridIt, T segmentNumber, bool doAgain = false) {
		// the vector to contain a list of the yet-to-search points
		typename std::vector< typename std::vector<SG::candidate<T> >::iterator > searchpoints;
		
		// Find the gridpoints corresponding to the starting point iterator
		typename std::vector<SG::candidate<T> >::iterator startPoint = findGridpoint(gridIt->azimuth_, gridIt->zenith_);
		// check if the starting point is part of the foreground
		if (! (startPoint->proximity_ > 0.0)){
			std::cerr << "The gridPoint to find a segment around is not part of the foreground (the proximity_ is 0). No segments can be created in the background!" << std::endl;
			return 0;
		}
		
		// check if this gridpoint is already part of a segment
		if(startPoint->multiPurposeValue_ != 0.0){
			if (doAgain){
				std::cerr << "The gridPoint to find a segment around is already part of a segment!" << std::endl;
				std::cerr << "A new segment will be created containing the old. This works, but this is horribly inefficient!" << std::endl;
				std::cerr << "Try to avoid calling connectSegment on an already connected segment!" << std::endl;
			} else {
				return 0;
			}
		}
		
		// make the starting point the actual starting point for a segment
		startPoint->multiPurposeValue_ = segmentNumber;
		searchpoints.push_back(startPoint);
		
		// count how many points are marked for this label
		int numPoints = 0;
		// search as long as we still have unchecked gridpoints
		while (searchpoints.size()){
			++numPoints;
			// check one gridpoint in each loop
			typename std::vector<SG::candidate<T> >::iterator currentPoint = *searchpoints.begin();
		
// 			std::cout << searchpoints.size() << " " << findGridpointIndex(currentPoint) << std::endl;
			if (searchpoints.size() > size()){
				std::cerr << "ERROR: the number of remaining gridpoints cannot be greater than the number of gridpoints with multipurpose == 1.0" << std::endl;
				return -1;
			}
			
			// get the neighbours of the current searchpoint
			typename std::vector<SG::candidate<T> >::iterator itN,itS,itW,itE; // four are enough, the other ones 
			itN = n(currentPoint);
			itS = s(currentPoint);
			itW = w(currentPoint);
			itE = e(currentPoint);
			
			// now check neighbour points if 
			// they belong to the foreground == prox > 0
			// and haven't been marked already (by this segment (or any other which should not be possible))
			if (itW->proximity_ > 0.0 && itW->multiPurposeValue_ == 0.0){
				itW->multiPurposeValue_ = segmentNumber;
				searchpoints.push_back(itW);
			}
			if (itE->proximity_ > 0.0 && itE->multiPurposeValue_ == 0.0){
				itE->multiPurposeValue_ = segmentNumber;
				searchpoints.push_back(itE);
			}
			if (itN->proximity_ > 0.0 && itN->multiPurposeValue_ == 0.0){
				itN->multiPurposeValue_ = segmentNumber;
				searchpoints.push_back(itN);
			}
			if (itS->proximity_ > 0.0 && itS->multiPurposeValue_ == 0.0){
				itS->multiPurposeValue_ = segmentNumber;
				searchpoints.push_back(itS);
			}
			
			// remove the currently searched gridpoint from the search vector
			searchpoints.erase(searchpoints.begin());
// 			searchpoints.erase(currentPoint);
		}
		
// 		std::cerr << "Done connectSegment, #points: " << numPoints << std::endl;
		return numPoints;
	}
	
	void keepTop(double fraction)
	{
		if (fraction <= 0.0)
		{
			fraction = 0.0;
			return;
		}
		else if (fraction >= 1.0)
		{
			return;
		}
		if (size() == 0)
		{
			std::cerr << "There is nothing to segment in an empty grid." << std::endl;
			return;
		}

		std::vector<double> allProximities;
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it)
		{
			allProximities.push_back(it->proximity_);
		}
		
		std::sort(allProximities.begin(), allProximities.end());
		double thres = allProximities[ static_cast<unsigned int>( (1.0-fraction)*allProximities.size() ) ];
		
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it)
		{
			if (it->proximity_ < thres)
			{
				it->proximity_ = 0.0;
			}
		}
// 		std::cout << "DEBUG fraction=" << fraction << " allProximities.size()=" << allProximities.size() << " thres=" << thres << " c=" << c << std::endl;
	}
	
	unsigned int size() const
	{
		return gridpoints_.size();
	}
	
	inline SG::candidate<T>& operator[](unsigned int i)
	{
		return gridpoints_[i];
	}
	
	inline const SG::candidate<T>& operator[](unsigned int i) const
	{
		return gridpoints_[i];
	}
	
	typename std::vector<SG::candidate<T> >::iterator findGridpoint(typename std::vector<SG::candidate<T> >::const_iterator it)
	{
		return gridpoints_.begin()+findGridpointIndex(it->azimuth_, it->zenith_);
	}
	
	typename std::vector<SG::candidate<T> >::const_iterator findGridpoint(typename std::vector<SG::candidate<T> >::const_iterator it) const
	{
		return gridpoints_.begin()+findGridpointIndex(it->azimuth_, it->zenith_);
	}
	
	typename std::vector<SG::candidate<T> >::iterator findGridpoint(T az, T zen){
		return gridpoints_.begin()+findGridpointIndex(az, zen);
	}
	
	typename std::vector<SG::candidate<T> >::const_iterator findGridpoint(T az, T zen) const{
		return gridpoints_.begin()+findGridpointIndex(az, zen);
	}
	
	void highpassFilterCannyStandard(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		
	// 	// scharr x
		gNW = 3.0; gN = 0.0; gNE = -3.0;
		gW = 10.0; gC = 0.0; gE = -10.0;
		gSW = 3.0; gS = 0.0; gSE = -3.0;

		//sobel y
	// 	gNW = 0.0; gN = 0.0; gNE = 0.0;
	// 	gW = 1.0; gC = 0.0; gE = -1.0;
	// 	gSW = 0.0; gS = 0.0; gSE = 0.0;

		std::vector<SG::candidate<T> > gridpointsFilteredX;
		filter(gridpointsFilteredX, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		
	// 	// scharr y
		gNW = 3.0; gN = 10.0; gNE = 3.0;
		gW = 0.0; gC = 0.0; gE = 0.0;
		gSW = -3.0; gS = -10.0; gSE = -3.0;

		// sobel y
	// 	gNW = 0.0; gN = 1.0; gNE = 0.0;
	// 	gW = 0.0; gC = 0.0; gE = 0.0;
	// 	gSW = 0.0; gS = -1.0; gSE = 0.0;
		
		std::vector<SG::candidate<T> > gridpointsFilteredY;
		filter(gridpointsFilteredY, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		
		typename std::vector<SG::candidate<T> >::iterator itX = gridpointsFilteredX.begin();
		typename std::vector<SG::candidate<T> >::iterator itY = gridpointsFilteredY.begin();
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it){
			it->proximity_ = sqrt(itX->proximity_*itX->proximity_ + itY->proximity_*itY->proximity_);
			++itX;
			++itY;
		}
	}
	
	void highpassFilter(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		
	// 	// scharr x
		gNW = 3.0; gN = 0.0; gNE = -3.0;
		gW = 10.0; gC = 0.0; gE = -10.0;
		gSW = 3.0; gS = 0.0; gSE = -3.0;

		std::vector<SG::candidate<T> > gridpointsFilteredX;
		filter(gridpointsFilteredX, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		
	// 	// scharr y
		gNW = 3.0; gN = 10.0; gNE = 3.0;
		gW = 0.0; gC = 0.0; gE = 0.0;
		gSW = -3.0; gS = -10.0; gSE = -3.0;
		
		std::vector<SG::candidate<T> > gridpointsFilteredY;
		filter(gridpointsFilteredY, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		
		// scharr-like 3
		gNW = 10.0; gN = 3.0; gNE = 0.0;
		gW = 3.0; gC = 0.0; gE = -3.0;
		gSW = 0.0; gS = -3.0; gSE = -10.0;

		std::vector<SG::candidate<T> > gridpointsFiltered3;
		filter(gridpointsFiltered3, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		
		// scharr-like 4
		gNW = 0.0; gN = 3.0; gNE = 10.0;
		gW = -3.0; gC = 0.0; gE = 3.0;
		gSW = -10.0; gS = -3.0; gSE = 0.0;

		std::vector<SG::candidate<T> > gridpointsFiltered4;
		filter(gridpointsFiltered4, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		
		typename std::vector<SG::candidate<T> >::iterator itX = gridpointsFilteredX.begin();
		typename std::vector<SG::candidate<T> >::iterator itY = gridpointsFilteredY.begin();
		typename std::vector<SG::candidate<T> >::iterator it3 = gridpointsFiltered3.begin();
		typename std::vector<SG::candidate<T> >::iterator it4 = gridpointsFiltered4.begin();
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it){
			it->proximity_ = sqrt(itX->proximity_*itX->proximity_ + itY->proximity_*itY->proximity_ + it3->proximity_*it3->proximity_ + it4->proximity_*it4->proximity_);
			++itX;
			++itY;
			++it3;
			++it4;
		}
	}
	
	// high-pass filter the grid proximity values
	void laplaceFilter(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		
		// Laplace
		gNW = 1.0; gN = 1.0; gNE = 1.0;
		gW = 1.0; gC = -8.0; gE = 1.0;
		gSW = 1.0; gS = 1.0; gSE = 1.0;
		
		std::vector<SG::candidate<T> > gridpointsFiltered;
		filter(gridpointsFiltered, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		gridpoints_ = gridpointsFiltered;
	}
	
	// Gaussian filter the grid
	void gaussianFilter(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		
		T factor = 1.0/16.0;
		// gaussian
		gNW = 1.0*factor; gN = 2.0*factor; gNE = 1.0*factor;
		gW = 2.0*factor; gC = 4.0*factor; gE = 2.0*factor;
		gSW = 1.0*factor; gS = 2.0*factor; gSE = 1.0*factor;
		
		std::vector<SG::candidate<T> > gridpointsFiltered;
		filter(gridpointsFiltered, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		gridpoints_ = gridpointsFiltered;
	}
	
	// normalize proximity_
	void normalize(){
		T max = std::numeric_limits<T>::min();
		T min = std::numeric_limits<T>::max();
		for (typename std::vector<SG::candidate<T> >::const_iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it){
			if (it->proximity_ > max) max = it->proximity_;
			if (it->proximity_ < min) min = it->proximity_;
		}
		if (max-min == 0.0) {
			for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it)
				it->proximity_ = 0.0;
		} else {
			for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it)
				it->proximity_ = (it->proximity_-min)/(max-min);
		}
	}


	// smear filter the grid with a high emphasis of the centered value
	void smearCentered(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		
		T factor = 1.0/20.0;
// 		gNW = 1.0*factor; gN = 2.0*factor; gNE = 1.0*factor;
// 		gW = 2.0*factor; gC = 8.0*factor; gE = 2.0*factor;
// 		gSW = 1.0*factor; gS = 2.0*factor; gSE = 1.0*factor;
		gNW = 1.0*factor; gN = 1.0*factor; gNE = 1.0*factor;
		gW = 1.0*factor; gC = 12.0*factor; gE = 1.0*factor;
		gSW = 1.0*factor; gS = 1.0*factor; gSE = 1.0*factor;
		
		std::vector<SG::candidate<T> > gridpointsFiltered;	// TODO it seems to me that this filtering isn't even used! the function stores the result in gridpoints_, but this is then copied over by the original again
		filter(gridpointsFiltered, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		gridpoints_ = gridpointsFiltered;
	}
	
	
	void LoG(){	// edge detection using laplassian of gaussian LoG filter
		gaussianFilter();
		laplaceFilter();
	}
	
	// Lowpass filter the grid
	void lowpassFilter(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		T factor = 1.0/9.0;
		
		// LowPass
		gNW = factor; gN = factor; gNE = factor;
		gW = factor; gC = factor; gE = factor;
		gSW = factor; gS = factor; gSE = factor;
		
		std::vector<SG::candidate<T> > gridpointsFiltered;
		filter(gridpointsFiltered, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		gridpoints_ = gridpointsFiltered;
	}
	
	// Lowpass filter the grid
	void lowpassFilterLeftRight(){
		T gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE;
		T factor = 1.0/3.0;
		T factor0 = 0.0;
		
		// LowPass
		gNW = factor0; gN = factor0; gNE = factor0;
		gW = factor; gC = factor; gE = factor;
		gSW = factor0; gS = factor0; gSE = factor0;
		
		std::vector<SG::candidate<T> > gridpointsFiltered;
		filter(gridpointsFiltered, gNW, gN, gNE, gW, gC, gE, gSW, gS, gSE);
		gridpoints_ = gridpointsFiltered;
	}
	
	// This function returns iterators to the neighbours of this gridpoint
	// It assumes the geometry hasn't been changed (e.g. by a sort)
	// this is the slow version, but the other, faster one seems to become unmaintainable
	int getNeighbours(const typename std::vector<SG::candidate<T> >::const_iterator gC,
		typename std::vector<SG::candidate<T> >::iterator& gNW, typename std::vector<SG::candidate<T> >::iterator& gN, typename std::vector<SG::candidate<T> >::iterator& gNE, 
		typename std::vector<SG::candidate<T> >::iterator& gW,                                                         typename std::vector<SG::candidate<T> >::iterator& gE, 
		typename std::vector<SG::candidate<T> >::iterator& gSW, typename std::vector<SG::candidate<T> >::iterator& gS, typename std::vector<SG::candidate<T> >::iterator& gSE)
	{
		
		// tiny checks to prevent usage on resorted grid
		double minDistance = 1.5*M_PI/(18*fine_);	// more the distance they should have, but meant as a margin
// 		assert(distanceRad(gridpoints_.begin(), (gridpoints_.begin()+1)) < minDistance);
// 		assert(distanceRad(gridpoints_.begin()+0.5*gridpoints_.size(), (gridpoints_.begin()+0.5*gridpoints_.size()+1)) < minDistance);
		
		if (! (distanceRad(*gridpoints_.begin(), *(gridpoints_.begin()+1)) < minDistance)){
			std::cerr << "distance assert failed: " << distanceRad(*gridpoints_.begin(), *(gridpoints_.begin()+1)) << "<" << minDistance << std::endl;
			throw(-1);
		}
		if(! (distanceRad(*(gridpoints_.begin()+0.5*gridpoints_.size()), *(gridpoints_.begin()+0.5*gridpoints_.size()+1)) < minDistance)){
			std::cerr << "distance assert failed: " << distanceRad(*(gridpoints_.begin()+0.5*gridpoints_.size()), *(gridpoints_.begin()+0.5*gridpoints_.size()+1)) << "<" << minDistance << std::endl;
			throw(-1);
		}
		
		gN = n(gC);
		gNE = e(n(gC));
		gNW = w(n(gC));
		gS = s(gC);
		gSE = e(s(gC));
		gSW = w(s(gC));
		gE = e(gC);
		gW = w(gC);
		
		return 0;
	}
	
	unsigned int countMultipurposeValue(T value){
		unsigned int counter = 0;
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpoints_.begin(); it != gridpoints_.end(); ++it){
			if (it->multiPurposeValue_ == value)
				++counter;
		}
		return counter;
	}

	T median(std::vector<T> vec){
		std::sort(vec.begin(), vec.end());
		return vec.at(0.5*vec.size()+0.5);
	}
	
	T median(T nw, T n, T ne, T w, T c, T e, T sw, T s, T se){
		std::vector<T> vec;
		vec.push_back(nw);
		vec.push_back(n);
		vec.push_back(ne);
		vec.push_back(w);
		vec.push_back(c);
		vec.push_back(e);
		vec.push_back(sw);
		vec.push_back(s);
		vec.push_back(se);
		return median(vec);
	}
	
	T averageNoZeros(T nw, T n, T ne, T w, T c, T e, T sw, T s, T se){
		std::vector<T> vec;
		vec.push_back(nw);
		vec.push_back(n);
		vec.push_back(ne);
		vec.push_back(w);
		vec.push_back(c);
		vec.push_back(e);
		vec.push_back(sw);
		vec.push_back(s);
		vec.push_back(se);
		
		double aver = 0.0;
		int count = 0;
		for (unsigned int i = 0; i < vec.size(); ++i){
			aver += vec[i];
			if (vec[i] != 0.0){
				++count;
			}
		}
		if (count > 0){
			aver /= (double)count;
		}
		return aver;
	}
	
	void medianFilter(){
		std::vector<SG::candidate<T> > gridpointsFiltered;
		medianFilterInternal(gridpointsFiltered);
		gridpoints_ = gridpointsFiltered;
	}
	
	void fillZeros(){
		std::vector<SG::candidate<T> > gridpointsFiltered;
		fillZeroesInternal(gridpointsFiltered);
		gridpoints_ = gridpointsFiltered;
	}
	
	void bkpProximity(){
		for (unsigned int i = 0; i < gridpoints_.size(); ++i){
			gridpoints_[i].multiPurposeValue_ = gridpoints_[i].proximity_;
		}
	}
	
	void connectMaxes(){
		std::vector<double> proximities;
		for (unsigned int i = 0; i < gridpoints_.size(); ++i){
			proximities.push_back(gridpoints_[i].proximity_);
		}
		
		std::sort(proximities.begin(), proximities.end());
		double maxVal = 0.0;
		maxVal = proximities[proximities.size()-1];
		double connectionVal = 0.0;
		for (int i = proximities.size()-1; i > 0; --i){
			if (proximities[i] != maxVal){
				connectionVal = proximities[i];
				break;
			}
		}
		
		for (unsigned int i = 0; i < gridpoints_.size(); ++i){
			if (gridpoints_[i].proximity_ < connectionVal){
				gridpoints_[i].proximity_ = 0.0;
			}
		}

	}
	
	int writeHammerProjection(std::string outputFileName)
	{
		return SG::writeHammerProjection(gridpoints_, outputFileName);
	}
	
	void writeAsEquatorialHammerProj(std::string outputFolder, std::string fileName = "hammerProjEquatorial.txt")
	{
		correctEquatorialRightAscensionRangeFranzoesisch();
		correctEquatorialRightAscensionRangeMinusPi();
		invertEquatorialRightAscensionRange();
		writeHammerProjection(outputFolder+"/"+fileName);
	}

	void writeAsGalacticHammerProj(std::string outputFolder, std::string fileName = "hammerProjGalactic.txt")
	{
		correctEquatorialRightAscensionRangeFranzoesisch();
		equatorialToGalactic();
		correctEquatorialRightAscensionRangeFranzoesischInvers();
		invertEquatorialRightAscensionRange();
		writeHammerProjection(outputFolder+"/"+"hammerProjGalactic.txt");
	}

private:
	unsigned int fine_;
	
	typename std::vector<SG::candidate<T> >::iterator w(typename std::vector<SG::candidate<T> >::const_iterator it){
		unsigned int zenIndex, azIndex;
		findGridpointIndices(it, zenIndex, azIndex);
		return w(zenIndex, azIndex);
	}
	
	typename std::vector<SG::candidate<T> >::iterator e(typename std::vector<SG::candidate<T> >::const_iterator it){
		unsigned int zenIndex, azIndex;
		findGridpointIndices(it, zenIndex, azIndex);
		return e(zenIndex, azIndex);
	}
	
	typename std::vector<SG::candidate<T> >::iterator n(typename std::vector<SG::candidate<T> >::const_iterator it){
		unsigned int zenIndex, azIndex;
		findGridpointIndices(it, zenIndex, azIndex);
		return n(zenIndex, azIndex);
	}
	
	typename std::vector<SG::candidate<T> >::iterator s(typename std::vector<SG::candidate<T> >::const_iterator it){
		unsigned int zenIndex, azIndex;
		findGridpointIndices(it, zenIndex, azIndex);
		return s(zenIndex, azIndex);
	}
	
	typename std::vector<SG::candidate<T> >::iterator w(unsigned int zenIndex, unsigned int azIndex, unsigned int baseIndex = std::numeric_limits<unsigned int>::max() ){
		assert(zenIndex < fine_*18+1);
		unsigned int newAzIndex = 0;
		if (azIndex > 0){
			newAzIndex = azIndex-1;
		} else {
			double tempZenith = static_cast<double>(M_PI*zenIndex/(18*fine_)) - 0.5*M_PI;
			// newAzIndex = numAz-1
			// The line below was the old code before debugging
// 			newAzIndex = static_cast<unsigned int>(fine_*36*cos(tempZenith) + 0.5)-1;
			// The lines after are the new code?
			unsigned int numAz = static_cast<unsigned int>(fine_*36*cos(tempZenith) + 0.5);
// 			unsigned int numAz2 = static_cast<unsigned int>(fine_*36*cos(tempZenith) + 0.5);
			if (numAz == 0) numAz = 1;
			newAzIndex = numAz-1;
// 			if (azIndex == 0){
// 				std::cerr << "west: " << zenIndex << " " << azIndex << " " << numAz2 << " " << tempZenith << " " << fine_*36*cos(tempZenith) << " " << newAzIndex << std::endl;
// 			}
		}
		if (baseIndex > gridpoints_.size()){
			baseIndex = getBaseIndex(zenIndex);
		}
		return gridpoints_.begin()+baseIndex+newAzIndex;
	}
	
	typename std::vector<SG::candidate<T> >::iterator e(unsigned int zenIndex, unsigned int azIndex, unsigned int baseIndex = std::numeric_limits<unsigned int>::max() ){
		assert(zenIndex < fine_*18+1);
		unsigned int newAzIndex = 0;
		double tempZenith = static_cast<double>(M_PI*zenIndex/(18*fine_)) - 0.5*M_PI;
		unsigned int numAz = static_cast<unsigned int>(fine_*36*cos(tempZenith) + 0.5);
		if (numAz == 0) numAz = 1;
		assert(azIndex < numAz);
		if (azIndex < numAz-1){
			newAzIndex = azIndex+1;
		} else {
			newAzIndex = 0;
		}
		if (baseIndex > gridpoints_.size()){
			baseIndex = getBaseIndex(zenIndex);
		}
// 		if (baseIndex+newAzIndex == 165016){
// 			std::cerr << "east: " << tempZenith << " " << numAz << " " << azIndex << " " << zenIndex << " " << newAzIndex << " " << newAzIndex << std::endl;
// 		}
		return gridpoints_.begin()+baseIndex+newAzIndex;
	}
	
	typename std::vector<SG::candidate<T> >::iterator n(unsigned int zenIndex, unsigned int azIndex, unsigned int baseIndex = std::numeric_limits<unsigned int>::max() ){
		assert(zenIndex < fine_*18+1);
		if (baseIndex > gridpoints_.size()){
			baseIndex = getBaseIndex(zenIndex);
		}
		unsigned int newZenIndex = 0;
		unsigned int newAzIndex = 0;
		unsigned int newBaseIndex = 0;
		double tempZenith = static_cast<double>(M_PI*zenIndex/(18*fine_)) - 0.5*M_PI;
		unsigned int numAz = static_cast<unsigned int>(fine_*36*cos(tempZenith) + 0.5);
		if (numAz == 0) numAz = 1;
		double whereWeAre = static_cast<double>(azIndex)/static_cast<double>(numAz);
// 		double whereWeAre = static_cast<double>(azIndex)/static_cast<double>(numAz+1);
		assert (whereWeAre >= 0.0);
		assert (whereWeAre <= 1.0);
		
		if (zenIndex > 0){
			newZenIndex = zenIndex-1;
			
			double tempZenithN = static_cast<double>(M_PI*newZenIndex/(18*fine_)) - 0.5*M_PI;
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*cos(tempZenithN) + 0.5);
			if (numAzN == 0) numAzN = 1;
			if (numAzN > 1){
				newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
			} else {
				newAzIndex = 0;
			}
// 			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN));
			
			assert(baseIndex >= numAzN);
			newBaseIndex = baseIndex-numAzN;
		} else {
			// what is north of the northest point of a sphere?
			// NOTE: It seems unclear whether the function should return 
			// - the gridpoint itself (but this changes behaviour compared to all other points unnoticed)
			// - an error (logical, but helps no one)
			// - a point on the backside of the ring one south (good logic for filtering and 
			// if you don't want to care where you are on the sphere, but weired for other purposes and
			// not neccesarily consistent which point!)
			
			// Current solution: Try to treat it similar to any other point on the sphere:
			// give the northest point "virtual northern neighbours" in the south with 180deg shifted azimuth
			
			// this sets the "virtual north" to one south
			newZenIndex = 1;
			// this shifts the azimuth for 180deg
			whereWeAre += 0.5;
			// this keeps the azimuth between 0 and 360 deg
			if (whereWeAre >= 1.0)
				whereWeAre -= 1.0;
			double tempZenithN = static_cast<double>(M_PI*newZenIndex/(18*fine_)) - 0.5*M_PI;
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*cos(tempZenithN) + 0.5);
			if (numAzN == 0) numAzN = 1;
			if (numAzN > 1){
				newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
			} else {
				newAzIndex = 0;
			}
// 			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN));
			
			assert (baseIndex == 0);
			newBaseIndex = numAz;
		}
		
		return gridpoints_.begin()+newBaseIndex+newAzIndex;
	}
	
	typename std::vector<SG::candidate<T> >::iterator s(unsigned int zenIndex, unsigned int azIndex, unsigned int baseIndex = std::numeric_limits<unsigned int>::max() ){
		assert(zenIndex < fine_*18+1);
		if (baseIndex > gridpoints_.size()){
			baseIndex = getBaseIndex(zenIndex);
		}
		unsigned int newZenIndex = 0;
		unsigned int newAzIndex = 0;
		unsigned int newBaseIndex = 0;
		double tempZenith = static_cast<double>(M_PI*zenIndex/(18*fine_)) - 0.5*M_PI;
		unsigned int numAz = static_cast<unsigned int>(fine_*36*cos(tempZenith) + 0.5);
		if (numAz == 0) numAz = 1;
		double whereWeAre = static_cast<double>(azIndex)/static_cast<double>(numAz);
// 		double whereWeAre = static_cast<double>(azIndex)/static_cast<double>(numAz+1);
		assert (whereWeAre >= 0.0);
		assert (whereWeAre <= 1.0);
		unsigned int numAzS = 0;
		
		if (zenIndex < fine_*18){
			newZenIndex = zenIndex+1;
			
			double tempZenithS = static_cast<double>(M_PI*newZenIndex/(18*fine_)) - 0.5*M_PI;
			numAzS = static_cast<unsigned int>(fine_*36*cos(tempZenithS) + 0.5);
// 			unsigned int numAzS = static_cast<unsigned int>(fine_*36*cos(tempZenithS) + 0.5);
			if (numAzS == 0) numAzS = 1;
			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
// 			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS-1) +0.5);
// 			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS));
			if (newAzIndex >= numAzS){
				newAzIndex = numAzS-1;
			}
			
			assert(baseIndex < gridpoints_.size()-numAz);
			newBaseIndex = baseIndex+numAz;
		} else {
			// what is south of the southest point of a sphere?
			// NOTE: See north function
			// Current solution: Try to treat it similar to any other point on the sphere:
			// give the southest point "virtual southern neighbours" in the north with 180deg shifted azimuth
			
			// this sets the "virtual south" to one north
			newZenIndex = fine_*18-1;
			// this shifts the azimuth for 180deg
			whereWeAre += 0.5;
			// this keeps the aziuth between 0 and 360 deg
			if (whereWeAre >= 1.0)
				whereWeAre -= 1.0;
			double tempZenithS = static_cast<double>(M_PI*newZenIndex/(18*fine_)) - 0.5*M_PI;
			numAzS = static_cast<unsigned int>(fine_*36*cos(tempZenithS) + 0.5);
// 			unsigned int numAzS = static_cast<unsigned int>(fine_*36*cos(tempZenithS) + 0.5);
			if (numAzS == 0) numAzS = 1;
			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
// 			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS-1) +0.5);
// 			newAzIndex = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS));
			
			if (newAzIndex >= numAzS){
				newAzIndex = numAzS-1;
			}
			
			assert (baseIndex == gridpoints_.size()-numAz);
			newBaseIndex = baseIndex-numAzS;
		}
		
		return gridpoints_.begin()+newBaseIndex+newAzIndex;
	}
	
	
public:	
	unsigned int findGridpointIndex(typename std::vector<SG::candidate<T> >::const_iterator it) const {
		return findGridpointIndex(it->azimuth_, it->zenith_);
	}
	unsigned int findGridpointIndex(T az, T zen) const {
		unsigned int zenIndex = 0;
		unsigned int azIndex = 0;
		findGridpointIndices(az, zen, zenIndex, azIndex);
		unsigned int finalIndex = getBaseIndex(zenIndex) + azIndex;
		assert(finalIndex < gridpoints_.size());
		return finalIndex;
	}
private:
	
	inline void findGridpointIndices(typename std::vector<SG::candidate<T> >::const_iterator it, unsigned int& zenIndex, unsigned int& azIndex) const {
		findGridpointIndices(it->azimuth_, it->zenith_, zenIndex, azIndex);
	}
	
	inline void findGridpointIndices(T az, T zen, unsigned int& zenIndex, unsigned int& azIndex) const {
		// To find the indicies, we reverse the computations done when setting up the grid
		zenIndex = static_cast<unsigned int>( ((zen+0.5*M_PI)*18*fine_)/M_PI + 0.5 );
		assert(zenIndex < fine_*18+1);
		// the zenith was simple! Now the azimuth
		unsigned int numAz = static_cast<unsigned int>(fine_*36*cos(zen) + 0.5);
		if (numAz == 0)
			numAz = 1;
		azIndex = static_cast<unsigned int>( (az+M_PI)*numAz/(2.0*M_PI) +0.5 );
		// Now we know where we are in the grid, zenith and azimuth
		assert(azIndex < numAz);
	}
	
	inline unsigned int getBaseIndex(unsigned int zenIndex) const {
		assert(zenIndex < fine_*18+1);
		// To find the indicies, we reverse the computations done when setting up the grid
		unsigned int baseIndex = 0;
		for (unsigned int j = 0; j < zenIndex; ++j){
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor + 0.5);
			if (numAz == 0) 
				numAz = 1;
			baseIndex += numAz;
		}
		assert(baseIndex < gridpoints_.size());
		return baseIndex;
	}
	
	void filter(std::vector<SG::candidate<T> >& gridpointsFiltered,
		    T gNW, T gN, T gNE, 
		    T gW,  T gC, T gE, 
		    T gSW, T gS, T gSE){
		gridpointsFiltered = gridpoints_;
		
		// Should be unneccessary, but make sure:
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpointsFiltered.begin(); it != gridpointsFiltered.end(); ++it)
			it->proximity_ = 0.0;
		
		//    		nw	north	ne
		//  j |		west    center  east
		//    v		sw	south	se
		//
		//			i ->
		
		unsigned int nAzBaseIndex = 0;
		unsigned int cAzBaseIndex = 0;
		unsigned int sAzBaseIndex = 0;
		unsigned int cAzBaseIndexBkp = 0;
		{
			unsigned int j = 0;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			if (numAz == 0) numAz = 1;
			
			cAzBaseIndex += numAz;
			sAzBaseIndex += numAz;
			cAzBaseIndexBkp = numAz;
			
			{
				unsigned int j = 1;
				double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
				
				// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
				double densityFactor = cos(tempZenith);
				unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
				if (numAz == 0) numAz = 1;
				
				sAzBaseIndex += numAz;
			}
		}
		
		// The base indices for rows 0 (north), 1 (current center) and 2 (south) are set
		
		// center parts first
		for (unsigned int j = 1; j < fine_*18; ++j){
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithN = static_cast<double>(M_PI*(j-1)/(18*fine_)) - 0.5*M_PI;
			double tempZenithS = static_cast<double>(M_PI*(j+1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorN = cos(tempZenithN);
			double densityFactorS = cos(tempZenithS);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*densityFactorN+0.5);
			unsigned int numAzS = static_cast<unsigned int>(fine_*36*densityFactorS+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzN == 0) numAzN = 1;
			if (numAzS == 0) numAzS = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
				unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					assert(iN < numAzN);
	// 				if (iN == numAzN-1) iN = numAzN-2;
	// 				if (iN == 0) iN = 1;
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					assert(iS < numAzS);
	// 				if (iS == numAzS-1) iS = numAzS-2;
	// 				if (iS == 0) iS = 1;
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			{	// west border
				unsigned int i = 0;
	// 			double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
	// 			assert (whereWeAre >= 0.0);
	// 			assert (whereWeAre <= 1.0);
	// 			unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
	// 			unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				unsigned int iN = 0;	// we know where we are
				unsigned int iS = 0;
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN+numAzN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS+numAzS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			{	// east border
				unsigned int i = numAz-1;
	// 			double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
	// 			assert (whereWeAre >= 0.0);
	// 			assert (whereWeAre <= 1.0);
	// 			unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
	// 			unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				unsigned int iN = numAzN-1;	// we know where we are
				unsigned int iS = numAzS-1;
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN-numAzN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS-numAzS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			
			nAzBaseIndex += numAzN;
			cAzBaseIndex += numAz;
			sAzBaseIndex += numAzS;
		}

		{	// southern border = lower end of the sphere
			unsigned int j = fine_*18;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithN = static_cast<double>(M_PI*(j-1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorN = cos(tempZenithN);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*densityFactorN+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzN == 0) numAzN = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					assert(iN < numAzN);
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+1)%numAzN].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			{	// west border
				unsigned int i = 0;
				unsigned int iN = 0;	// we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN+numAzN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+numAzN+1)%numAz].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			{	// east border
				unsigned int i = numAz-1;
				unsigned int iN = numAzN-1; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN-numAzN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+1)%numAzN].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-numAzN-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
		}
		
		cAzBaseIndex = 0;
		sAzBaseIndex = cAzBaseIndexBkp;
		{	// norther boundary (the top of the sphere)
			unsigned int j = 0;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithS = static_cast<double>(M_PI*(j+1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorS = cos(tempZenithS);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzS = static_cast<unsigned int>(fine_*36*densityFactorS+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzS == 0) numAzS = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					assert(iS < numAzS);
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			{	// west border
				unsigned int i = 0;
				unsigned int iS = 0; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+numAzS+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS+numAzS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
			{	// east border
				unsigned int i = numAz-1;
				unsigned int iS = numAzS-1; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-numAzS-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS-numAzS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = gNW*nw + gN*n + gNE*ne + gW*w + gC*c + gE*e + gSW*sw + gS*s + gSE*se;
			}
		}
		// Done : )
		
	}	// end of function filter
	
	void medianFilterInternal(std::vector<SG::candidate<T> >& gridpointsFiltered){
		gridpointsFiltered = gridpoints_;
		
		// Should be unneccessary, but make sure:
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpointsFiltered.begin(); it != gridpointsFiltered.end(); ++it)
			it->proximity_ = 0.0;
		
		//    		nw	north	ne
		//  j |		west    center  east
		//    v		sw	south	se
		//
		//			i ->
		
		unsigned int nAzBaseIndex = 0;
		unsigned int cAzBaseIndex = 0;
		unsigned int sAzBaseIndex = 0;
		unsigned int cAzBaseIndexBkp = 0;
		{
			unsigned int j = 0;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			if (numAz == 0) numAz = 1;
			
			cAzBaseIndex += numAz;
			sAzBaseIndex += numAz;
			cAzBaseIndexBkp = numAz;
			
			{
				unsigned int j = 1;
				double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
				
				// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
				double densityFactor = cos(tempZenith);
				unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
				if (numAz == 0) numAz = 1;
				
				sAzBaseIndex += numAz;
			}
		}
		
		// The base indices for rows 0 (north), 1 (current center) and 2 (south) are set
		
		// center parts first
		for (unsigned int j = 1; j < fine_*18; ++j){
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithN = static_cast<double>(M_PI*(j-1)/(18*fine_)) - 0.5*M_PI;
			double tempZenithS = static_cast<double>(M_PI*(j+1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorN = cos(tempZenithN);
			double densityFactorS = cos(tempZenithS);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*densityFactorN+0.5);
			unsigned int numAzS = static_cast<unsigned int>(fine_*36*densityFactorS+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzN == 0) numAzN = 1;
			if (numAzS == 0) numAzS = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
				unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					assert(iN < numAzN);
	// 				if (iN == numAzN-1) iN = numAzN-2;
	// 				if (iN == 0) iN = 1;
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					assert(iS < numAzS);
	// 				if (iS == numAzS-1) iS = numAzS-2;
	// 				if (iS == 0) iS = 1;
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// west border
				unsigned int i = 0;
	// 			double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
	// 			assert (whereWeAre >= 0.0);
	// 			assert (whereWeAre <= 1.0);
	// 			unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
	// 			unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				unsigned int iN = 0;	// we know where we are
				unsigned int iS = 0;
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN+numAzN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS+numAzS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// east border
				unsigned int i = numAz-1;
	// 			double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
	// 			assert (whereWeAre >= 0.0);
	// 			assert (whereWeAre <= 1.0);
	// 			unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
	// 			unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				unsigned int iN = numAzN-1;	// we know where we are
				unsigned int iS = numAzS-1;
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN-numAzN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS-numAzS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			
			nAzBaseIndex += numAzN;
			cAzBaseIndex += numAz;
			sAzBaseIndex += numAzS;
		}

		{	// southern border = lower end of the sphere
			unsigned int j = fine_*18;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithN = static_cast<double>(M_PI*(j-1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorN = cos(tempZenithN);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*densityFactorN+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzN == 0) numAzN = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					assert(iN < numAzN);
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+1)%numAzN].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// west border
				unsigned int i = 0;
				unsigned int iN = 0;	// we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN+numAzN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+numAzN+1)%numAz].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// east border
				unsigned int i = numAz-1;
				unsigned int iN = numAzN-1; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN-numAzN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+1)%numAzN].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-numAzN-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
		}
		
		cAzBaseIndex = 0;
		sAzBaseIndex = cAzBaseIndexBkp;
		{	// norther boundary (the top of the sphere)
			unsigned int j = 0;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithS = static_cast<double>(M_PI*(j+1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorS = cos(tempZenithS);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzS = static_cast<unsigned int>(fine_*36*densityFactorS+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzS == 0) numAzS = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					assert(iS < numAzS);
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// west border
				unsigned int i = 0;
				unsigned int iS = 0; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+numAzS+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS+numAzS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// east border
				unsigned int i = numAz-1;
				unsigned int iS = numAzS-1; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-numAzS-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS-numAzS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = median(nw,n,ne,w,c,e,sw,s,se);
			}
		}
		// Done : )
		
	}	// end of function filter
	
	
	
	
	
	void fillZeroesInternal(std::vector<SG::candidate<T> >& gridpointsFiltered){
		gridpointsFiltered = gridpoints_;
		
		// Should be unneccessary, but make sure:
		for (typename std::vector<SG::candidate<T> >::iterator it = gridpointsFiltered.begin(); it != gridpointsFiltered.end(); ++it)
			it->proximity_ = 0.0;
		
		//    		nw	north	ne
		//  j |		west    center  east
		//    v		sw	south	se
		//
		//			i ->
		
		unsigned int nAzBaseIndex = 0;
		unsigned int cAzBaseIndex = 0;
		unsigned int sAzBaseIndex = 0;
		unsigned int cAzBaseIndexBkp = 0;
		{
			unsigned int j = 0;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			if (numAz == 0) numAz = 1;
			
			cAzBaseIndex += numAz;
			sAzBaseIndex += numAz;
			cAzBaseIndexBkp = numAz;
			
			{
				unsigned int j = 1;
				double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
				
				// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
				double densityFactor = cos(tempZenith);
				unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
				if (numAz == 0) numAz = 1;
				
				sAzBaseIndex += numAz;
			}
		}
		
		// The base indices for rows 0 (north), 1 (current center) and 2 (south) are set
		
		// center parts first
		for (unsigned int j = 1; j < fine_*18; ++j){
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithN = static_cast<double>(M_PI*(j-1)/(18*fine_)) - 0.5*M_PI;
			double tempZenithS = static_cast<double>(M_PI*(j+1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorN = cos(tempZenithN);
			double densityFactorS = cos(tempZenithS);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*densityFactorN+0.5);
			unsigned int numAzS = static_cast<unsigned int>(fine_*36*densityFactorS+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzN == 0) numAzN = 1;
			if (numAzS == 0) numAzS = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
				unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					assert(iN < numAzN);
	// 				if (iN == numAzN-1) iN = numAzN-2;
	// 				if (iN == 0) iN = 1;
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					assert(iS < numAzS);
	// 				if (iS == numAzS-1) iS = numAzS-2;
	// 				if (iS == 0) iS = 1;
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// west border
				unsigned int i = 0;
	// 			double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
	// 			assert (whereWeAre >= 0.0);
	// 			assert (whereWeAre <= 1.0);
	// 			unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
	// 			unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				unsigned int iN = 0;	// we know where we are
				unsigned int iS = 0;
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN+numAzN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS+numAzS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// east border
				unsigned int i = numAz-1;
	// 			double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
	// 			assert (whereWeAre >= 0.0);
	// 			assert (whereWeAre <= 1.0);
	// 			unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
	// 			unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				unsigned int iN = numAzN-1;	// we know where we are
				unsigned int iS = numAzS-1;
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN-numAzN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS-numAzS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			
			nAzBaseIndex += numAzN;
			cAzBaseIndex += numAz;
			sAzBaseIndex += numAzS;
		}

		{	// southern border = lower end of the sphere
			unsigned int j = fine_*18;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithN = static_cast<double>(M_PI*(j-1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorN = cos(tempZenithN);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzN = static_cast<unsigned int>(fine_*36*densityFactorN+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzN == 0) numAzN = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iN = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzN) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					assert(iN < numAzN);
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+1)%numAzN].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// west border
				unsigned int i = 0;
				unsigned int iN = 0;	// we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN+numAzN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+numAzN+1)%numAz].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// east border
				unsigned int i = numAz-1;
				unsigned int iN = numAzN-1; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzN >= 3){
					nw = gridpoints_[nAzBaseIndex+iN-1].proximity_;
					n = gridpoints_[nAzBaseIndex+iN].proximity_;
					ne = gridpoints_[nAzBaseIndex+iN-numAzN+1].proximity_;
				} else {
					nw = gridpoints_[nAzBaseIndex].proximity_;
					n = gridpoints_[nAzBaseIndex].proximity_;
					ne = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAzN >= 3){
					sw = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)+1)%numAzN].proximity_;
					s = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5))%numAzN].proximity_;
					se = gridpoints_[nAzBaseIndex+(iN+static_cast<unsigned int>(0.5*numAzN+0.5)-numAzN-1)%numAzN].proximity_;
				} else {
					sw = gridpoints_[nAzBaseIndex].proximity_;
					s = gridpoints_[nAzBaseIndex].proximity_;
					se = gridpoints_[nAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
		}
		
		cAzBaseIndex = 0;
		sAzBaseIndex = cAzBaseIndexBkp;
		{	// norther boundary (the top of the sphere)
			unsigned int j = 0;
			double tempZenith = static_cast<double>(M_PI*j/(18*fine_)) - 0.5*M_PI;
			double tempZenithS = static_cast<double>(M_PI*(j+1)/(18*fine_)) - 0.5*M_PI;
			
			// the density factor is neccessary because the grid points become closer in azimuth when the zenith is closer to the poles
			double densityFactor = cos(tempZenith);
			double densityFactorS = cos(tempZenithS);
			unsigned int numAz = static_cast<unsigned int>(fine_*36*densityFactor+0.5);
			unsigned int numAzS = static_cast<unsigned int>(fine_*36*densityFactorS+0.5);
			if (numAz == 0) numAz = 1;
			if (numAzS == 0) numAzS = 1;
			
			// middle entries first
			for (unsigned int i = 1; i < numAz-1; ++i){
				double whereWeAre = static_cast<double>(i)/static_cast<double>(numAz);
				assert (whereWeAre >= 0.0);
				assert (whereWeAre <= 1.0);
				unsigned int iS = static_cast<unsigned int>( whereWeAre*static_cast<double>(numAzS) +0.5);
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					assert(iS < numAzS);
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// west border
				unsigned int i = 0;
				unsigned int iS = 0; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+numAzS+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS+numAzS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i+numAz-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
			{	// east border
				unsigned int i = numAz-1;
				unsigned int iS = numAzS-1; // we know where we are
				T nw,n,ne,sw,s,se,w,c,e;
				
				if (numAzS >= 3){
					nw = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)+1)%numAzS].proximity_;
					n = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5))%numAzS].proximity_;
					ne = gridpoints_[sAzBaseIndex+(iS+static_cast<unsigned int>(0.5*numAzS+0.5)-numAzS-1)%numAzS].proximity_;
				} else {
					nw = gridpoints_[sAzBaseIndex].proximity_;
					n = gridpoints_[sAzBaseIndex].proximity_;
					ne = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAzS >= 3){
					sw = gridpoints_[sAzBaseIndex+iS-1].proximity_;
					s = gridpoints_[sAzBaseIndex+iS].proximity_;
					se = gridpoints_[sAzBaseIndex+iS-numAzS+1].proximity_;
				} else {
					sw = gridpoints_[sAzBaseIndex].proximity_;
					s = gridpoints_[sAzBaseIndex].proximity_;
					se = gridpoints_[sAzBaseIndex].proximity_;
				}
				if (numAz >= 3){
					w = gridpoints_[cAzBaseIndex+i-1].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i-numAz+1].proximity_;
				} else {
					w = gridpoints_[cAzBaseIndex+i].proximity_;
					c = gridpoints_[cAzBaseIndex+i].proximity_;
					e = gridpoints_[cAzBaseIndex+i].proximity_;
				}
				gridpointsFiltered[cAzBaseIndex+i].proximity_ = averageNoZeros(nw,n,ne,w,c,e,sw,s,se);
			}
		}
		// Done : )
		
	}	// end of function filter

};	// end of class sphericalGrid


template <typename T>
void renormalizeGrid(SG::sphericalGrid<T>& grid, T allSumBefore)
{
	T allSumAfter = 0.0;
	for (unsigned int i = 0; i < grid.size(); ++i)
	{
		allSumAfter += grid[i].proximity_;
	}
	T preserveFactor = allSumBefore/allSumAfter;
	for (unsigned int i = 0; i < grid.size(); ++i)
	{
		grid[i].proximity_ *= preserveFactor;
	}
}

template <typename T>
int save(SG::sphericalGrid<T> const & grid, std::string filename)
{
	std::ofstream raus;
	raus.open(filename.c_str());
	if (! raus.is_open()){
		std::cerr << "Failed to write spherical grid to file " << filename << "." << std::endl;
		return -1;
	}
	raus << grid.size() << " " << grid.fine() << std::endl;
	if (! raus.good()){
		std::cerr << "Failed to write spherical grid size to file " << filename << "." << std::endl;
		return -2;
	}
	for (unsigned int i = 0; i < grid.size(); ++i){
		raus << grid.gridpoints_[i] << std::endl;
		if (! raus.good()){
			std::cerr << "Failed to write spherical grid entry " << i << " to file " << filename << "." << std::endl;
			return -3;
		}
	}
	if (! raus.good()){
		std::cerr << "Failed to write spherical grid entry " << grid.size()-1 << " to file " << filename << "." << std::endl;
		return -3;
	}
	raus.close();
	return 0;
}

template <typename T>
int load(std::string filename, SG::sphericalGrid<T> & grid)
{
	std::ifstream rein;
	rein.open(filename.c_str());
	if (! rein.is_open()){
		std::cerr << "Failed to read spherical grid from file " << filename << "." << std::endl;
		return -1;
	}
	unsigned int num = 0;
	rein >> num;
	if (! rein.good()){
		std::cerr << "Failed to read spherical grid size from file " << filename << "." << std::endl;
		return -2;
	}
	std::cout << num << std::endl;
	unsigned int fine = 0;
	rein >> fine;
	if (! rein.good()){
		std::cerr << "Failed to read spherical grid fine from file " << filename << "." << std::endl;
		return -2;
	}
	std::cout << fine << std::endl;
	grid.setFineTEMPORARY(fine);

	grid.gridpoints_.clear();
	while (rein.good()){
		SG::candidate<T> gridpoint;
		rein >> gridpoint;
		if (rein.good()){
			grid.gridpoints_.push_back(gridpoint);
		}
	}
	rein.close();

	if (grid.gridpoints_.size() != num){
		std::cerr << "Successfully read " << grid.gridpoints_.size() << " gridpoints from file " << filename << " but expected " << num << " from file header!" << std::endl;
		return -3;
	}
	return 0;
}

// copied from the onoff project
/*
void writePGM(grid const & g, std::string const & fileName)
{
//	cout << "DEBUG: writePGM start" << endl;

	const int signal = 255;
//	const double delta = 0.5*(M_PI / 180.0);
	const unsigned int numX = 720;
	const unsigned int numY = 360;

	double deltaAz = (2.0*M_PI)/(numX-1);
	double azMin = -1.0*M_PI;
	double deltaZen = (1.0*M_PI)/(numY-1);
	double zenMin = -0.5*M_PI;

	// get maximum sig value of any point in g
	double maxValD = maxSignificance(g);

	// compose image
	std::vector<int> imageData(numX*numY, 0);

	// check all gridpoints
	for (unsigned int i = 0; i < g.size(); ++i)
	{
		// only do something for signal gridpoints
		if (g.get(i).sig_ > 0.0)
		{
			// find corresponding indices for this gridpoint
			double zen = g.get(i).zenith_;
			double az = g.get(i).azimuth_;
			int x = (az - azMin)/(deltaAz);
			int y = numY-( (zen - zenMin)/(deltaZen) );

			assert( x < static_cast<int>(numX) );
			assert( y < static_cast<int>(numY) );
			assert( x >= 0 );
			assert( y >= 0 );

//			imageData[y*numX+x] = signal;

			double s = g.get(i).sig_;
			assert (s <= maxValD);
			imageData[y*numX+x] = static_cast<unsigned int>( static_cast<double>(signal)*(s/maxValD) );
		}
	}

	// open stream
	std::ofstream raus(fileName.c_str());
	raus << "P2" << endl;
	if (! raus.good() ) throw "Writing to file failed";

	raus << numX << " " << numY << endl;
	raus << signal << endl;

	for (unsigned int i = 0; i < numY; ++i)
	{
		for (unsigned int j = 0; j < numX; ++j)
		{
			raus << imageData[i*numX+j] << " ";
		}
		raus << endl;
	}

	raus.close();

//	cout << "DEBUG: writePGM end" << endl;
}

//unsigned int numX = 720, unsigned int numY = 360
void readPGM(std::string const fileName, grid & g)
{
	cout << "DEBUG: readPGM start" << endl;
	g.clear();

	// open stream
	std::ifstream rein(fileName.c_str());
	std::string temp;
	rein >> temp;
	if (! rein.good() )
	{
		cerr << "Reading failed from file " << fileName << endl;
		throw "Reading from file failed";
	}
	if (temp != "P2")
	{
		cerr << "File " << fileName << " isn't a suited PGM file (it doesn't start with P2)" << endl;
		throw "File isn't a suited PGM file (it doesn't start with P2)";
	}

	// TODO handle comments starting with # in the pgm file (gimp uses comments) ... getline seems to be necessary?

	int signal = 1;
	unsigned int numX = 720;
	unsigned int numY = 360;

	rein >> numX;
	rein >> numY;
	rein >> signal;

	cout << numX << " " << numY << " " << signal << endl;

	double deltaAz = (2.0*M_PI)/(numX-1);
	double azMin = -1.0*M_PI;
	double deltaZen = (1.0*M_PI)/(numY-1);
	double zenMin = -0.5*M_PI;

	// compose image
	std::vector<int> imageData(numX*numY, 0);
	for (unsigned int i = 0; i < numY; ++i)
	{
		for (unsigned int j = 0; j < numX; ++j)
		{
			int temp;
			rein >> temp;
			if (! rein.good() )
			{
				cerr << "Reading failed from file " << fileName << " at " << i << " / " << j << endl;
				throw "Reading from file failed";
			}
			g.add( point(-1.0*(zenMin+i*deltaZen), azMin+j*deltaAz, temp) );
		}
	}

	rein.close();

	cout << "DEBUG: readPGM end" << endl;
}
*/

template <typename T>
int convertFromImage(std::string filename, SG::sphericalGrid<T> & grid)
{
	std::ifstream rein;
	rein.open(filename.c_str());
	if (! rein.is_open()){
		std::cerr << "Failed to read spherical grid from file " << filename << "." << std::endl;
		return -1;
	}
	unsigned int num = 0;
	rein >> num;
	if (! rein.good()){
		std::cerr << "Failed to read spherical grid size from file " << filename << "." << std::endl;
		return -2;
	}
	return 0;
}

} // end of namespace SG

#endif 
