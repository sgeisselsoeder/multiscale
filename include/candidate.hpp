/*
 * candidate.h
 *
 *  Created on: 21.09.2012
 *  Author: sgeisselsoeder
 */

#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
// TODO: not working without astro package
// #include "astro/astro.h"

namespace SG{

const double minAzimuth = -1.0*M_PI;
const double maxAzimuth = 1.0*M_PI;
const double minZenith = -0.5*M_PI;
const double maxZenith = +0.5*M_PI;
const double twoPi = 2.0*M_PI;

template <typename T>
class candidate{
private:
	void init(){
		azimuth_ = 0.0;
		zenith_ = 0.0;
		time_ = 0.0;
		energy_ = 0.0;
		proximity_ = 0.0;
		multiPurposeValue_ = 0.0;
	}

public:
	inline int set(T a, T z, T t = 0.0, T e = 0.0, T p = 0.0){
		if (z > maxZenith || z < minZenith){
			std::cerr << "The zenith value (" << z << ") must be within the range from " << minZenith << " (down from above) to " << maxZenith << " (up from below)." << std::endl;
// 			std::cerr << "Setting this zenith to 0.0 - BUT THIS IS A REAL ERROR!" << std::endl;
			std::cerr << "Dropping this event - BUT THIS IS A REAL ERROR!" << std::endl;
			std::cerr << "I'M SERIOUS ABOUT THIS! CHECK ALL ZENITH VALUES! THIS IS NOT A \"tiny warning-error\". YOU ARE DOING SOMETHING HORRIBLY WRONG!" << std::endl;
			std::cerr << "YOUR RESULTS WILL BE WRONG! - THE COLORFUL PICTURES => WRONG - THE PRESENTATION YOU WILL GIVE => WRONG - SO CHECK THE ZENITH VALUES!!!" << std::endl;
			return -2;
// 			z = 0.0;
		}
		zenith_ = z;

		azimuth_ = a;

		time_ = t;
		energy_ = e;
		proximity_ = p;
		
		return 0;
	}

	
	// Azimuth, -pi to pi, 0.0 is anywhere
	T azimuth_;
	// Zenith from pi/2 (exactly downgoing) to -pi/2 (exactly upgoing); therefore the horizon is at 0
	T zenith_;
	// The time of the event
	T time_;
	// The Energy of the event
	T energy_;
	// Integrated proximity of (other) candidates
	T proximity_;
	T multiPurposeValue_;

	candidate(){
		init();
	}

	candidate(T a, T z, T t = 0.0, T e = 0.0, T p = 0.0){
		set(a,z,t,e,p);
	}

	~candidate(){
		
	}
	
	inline T x() const {
		const T r = 1.0;
		//return r*sin(zenith_)*cos(azimuth_);
		return r*sin(zenith_+0.5*M_PI)*cos(azimuth_);
	}
	
	inline T y() const {
		const T r = 1.0;
		//return r*sin(zenith_)*sin(azimuth_);
		return r*sin(zenith_+0.5*M_PI)*sin(azimuth_);
	}
	
	inline T z() const {
		const T r = 1.0;
		//return r*cos(zenith_);
		return r*cos(zenith_+0.5*M_PI);
	}
	
	inline T xP(T maxProx) const {
		T r = 1.0+0.1*proximity_/maxProx;
		return r*sin(zenith_+0.5*M_PI)*cos(azimuth_);
	}

	inline T yP(T maxProx) const {
		T r = 1.0+0.1*proximity_/maxProx;
		return r*sin(zenith_+0.5*M_PI)*sin(azimuth_);
	}

	inline T zP(T maxProx) const {
		T r = 1.0+0.1*proximity_/maxProx;
		return r*cos(zenith_+0.5*M_PI);
	}


	void equatorialToGalactic(){
		// TODO: not working without astro package
/*
		Astro::EqCoords eqApparent;
		eqApparent.SetRaDecRad(azimuth_, zenith_);
		Astro::LongLatCoords galCoords = Astro::EquatorialJ2000ToGalactic1958(eqApparent);
		zenith_ = galCoords.GetLatRad();        // aka zenith
		azimuth_ = galCoords.GetLongRad();      // aka azimuth
*/  
      }

	void galacticToEquatorial(){
		// TODO: not working without astro package
/*
		Astro::LongLatCoords galCoords;
		galCoords.SetLongLatRad(azimuth_, zenith_);
		Astro::EqCoords eqApparent = Astro::Galactic1958ToEquatorialJ2000(galCoords);
		zenith_ = eqApparent.GetDecRad(); // declination in rad
		azimuth_ = eqApparent.GetRaRad(); // right ascension in rad
*/  
      }
	
	void correctEquatorialRightAscensionRangeFranzoesisch()
	{
		if (azimuth_ < 0.0){
			azimuth_ += 2.0*M_PI;
		}
    }

	void correctEquatorialRightAscensionRangeMinusPi()
	{
		azimuth_ -= M_PI;
	}

	void invertEquatorialRightAscensionRange()
	{
		azimuth_ *= -1.0;
	}

	void correctEquatorialRightAscensionRangeFranzoesischInvers()
	{
		if (azimuth_ > M_PI){
			azimuth_ -= 2.0*M_PI;
		}
	}

	// implementing the hammer projection
	// http://en.wikipedia.org/wiki/Hammer_projection
	// this is also known as hammer-aitov projection. 
	// Beware: projections and coordinate systems fuck with your head

	// (2.0*sqrt(2)*cos(zen)*sin(az/2)) / (sqrt(1.0+cos(zen)*cos(az/2)))
	inline T hammerProjectionX() const {
		double zenith = zenith_;
		double azimuth = azimuth_;
		
		return 2.0*sqrt(2.0)*cos(zenith)*sin(0.5*azimuth) / (sqrt(1.0+cos(zenith)*cos(0.5*azimuth)));
	}
	// sqrt(2) * sin(zen) / (sqrt(1.0+cos(zen)*cos(az/2)))
	inline T hammerProjectionY() const {
		double zenith = zenith_;
		double azimuth = azimuth_;

		return sqrt(2.0)*sin(zenith) / (sqrt(1.0+cos(zenith)*cos(0.5*azimuth)));
	}
	
	bool operator<(const candidate<T>& b) const {
		return compareProximity( *this, b);
	}
	
	bool samePos(const candidate<T>& b) const {
		if ((b.zenith_ == zenith_) && (b.azimuth_ == azimuth_)){
			return true;
		} else {
			return false;
		}
	}
	
	std::ostream& outputToStream(std::ostream& raus) const{
		raus << azimuth_ << " " << zenith_ <<  " ";
		raus.precision(20);
		raus << time_ << " ";
		raus.precision(6);
		raus << energy_ << " " << proximity_ << " " << multiPurposeValue_;
		return raus;
	}
	
	std::istream& inputFromStream(std::istream& rein){
		rein >> azimuth_ >> zenith_ >> time_ >> energy_ >> proximity_ >> multiPurposeValue_;
		return rein;
	}
};

typedef candidate<double> candidate_t;

template <typename T>
std::ostream& operator<<(std::ostream& lhs, const candidate<T>& rhs)
{
	return rhs.outputToStream(lhs);
}

template <typename T>
std::istream& operator>>(std::istream& lhs, candidate<T>& rhs)
{
	return rhs.inputFromStream(lhs);
}

template <typename T>
void readCandidatesFromFile(std::string const & fileName, std::vector<SG::candidate<T> >& candidates)
{
	std::cout << "Reading candidates from file " << fileName << ". " << std::endl;
	candidates.clear();
	std::ifstream rein(fileName.c_str());
	if (! rein.is_open() || ! rein.good())
	{
		std::cerr << "Failed to open file " << fileName << " to read candidates." << std::endl;
		throw "ball";
	}
	while (rein.good())
	{
		SG::candidate<double> c;
		rein >> c;
		if (rein.good())
		{
			candidates.push_back(c);
		}
	}
	rein.close();
	std::cout << "Read " << candidates.size() << " candidates." << std::endl;
}

template <typename T>
bool compareProximity(const candidate<T>& a, const candidate<T>& b){
	return (a.proximity_ > b.proximity_);
}

template <typename T>
T quickDist(const candidate<T>& a, const candidate<T>& b){
	assert(a.zenith_ >= -0.5*M_PI);
	assert(a.zenith_ <= 0.5*M_PI);
	assert(b.zenith_ >= -0.5*M_PI);
	assert(b.zenith_ <= 0.5*M_PI);

	double zenithDiff = fabs(b.zenith_ - a.zenith_);

	double azimuthDiff = fabs(b.azimuth_ - a.azimuth_);
	double azimuthDiff2 = fabs(b.azimuth_+2.0*M_PI - a.azimuth_);
	double azimuthDiff3 = fabs(b.azimuth_-2.0*M_PI - a.azimuth_);

	return std::min(std::min(zenithDiff, azimuthDiff), std::min(azimuthDiff2, azimuthDiff3));
}

template <typename T>
inline T distanceRad(const candidate<T>& a, const candidate<T>& b){
	return greatCircleDistance(a,b);
}

// The great circle distance is used to compute the distance between two points on a sphere
// in (BOGENMASS), breite = zenith, laenge = azimuth
// this implements the numerically more stable computation
// C = 2 * arcsin( sqrt( sin^2( (zenith_B - zenith_A)/2) + cos(zenith_A)*cos(zenith_B)*sin^2( (azimuth_B - azimuth_A)/2) )
// from http://en.wikipedia.org/wiki/Great-circle_distance#Formulas
template <typename T>
inline T greatCircleDistance(const candidate<T>& a, const candidate<T>& b){
	assert(a.zenith_ >= -0.5*M_PI);
	assert(a.zenith_ <= 0.5*M_PI);
	assert(b.zenith_ >= -0.5*M_PI);
	assert(b.zenith_ <= 0.5*M_PI);

	T zenithDiff = sin(0.5*(b.zenith_ - a.zenith_));
	zenithDiff *= zenithDiff;
	T azimuthDiff = sin(0.5*(b.azimuth_ - a.azimuth_));
	azimuthDiff *= azimuthDiff;

	return 2.0*asin(sqrt(zenithDiff + cos(a.zenith_)*cos(b.zenith_)*azimuthDiff));
}


template <typename T>
void findRange(const std::vector<SG::candidate<T> >& points, T& minAz, T& maxAz, T& minZen, T& maxZen){
	minZen = 9999999999.9;
	maxZen = -9999999999.9;
	minAz = 9999999999.9;
	maxAz = -9999999999.9;
	for (unsigned int i = 0; i < points.size(); ++i){
		double az = points[i].azimuth_;
		double zen = points[i].zenith_;
		if (az < minAz) minAz = az;
		if (az > maxAz) maxAz = az;
		if (zen < minZen) minZen = zen;
		if (zen > maxZen) maxZen = zen;
	}
}

template <typename T>
void findRange(const std::vector<SG::candidate<T> >& points){
	double minZen = 9999999999.9;
	double maxZen = -9999999999.9;
	double minAz = 9999999999.9;
	double maxAz = -9999999999.9;
	//std::cout << "Analyzing " << points.size() << " points:" << std::endl;
	findRange(points, minAz, maxAz, minZen, maxZen);	
	std::cout << "Azimuth from " << minAz << " to " << maxAz << std::endl;
	std::cout << "Zenith from " << minZen << " to " << maxZen << std::endl;
}

template <typename T>
int isInRange(const std::vector<SG::candidate<T> >& points, T minAzVal, T maxAzVal, T minZenVal, T maxZenVal){
	const double eps = 0.001;
	double minZen = 9999999999.9;
	double maxZen = -9999999999.9; 
	double minAz = 9999999999.9; 
	double maxAz = -9999999999.9;
//	std::cout << "Analyzing " << points.size() << " points:" << std::endl;
	findRange(points, minAz, maxAz, minZen, maxZen);	
	bool errorFlag = false;
	if (minZen < minZenVal-eps) {
		std::cerr << "minZen = " << minZen << std::endl;
		errorFlag = true;
	}
	if (maxZen > maxZenVal+eps) {
		std::cerr << "maxZen = " << maxZen << std::endl;
		errorFlag = true;
	}
	if (minAz < minAzVal-eps) {
		std::cerr << "minAz = " << minAz << std::endl;
		errorFlag = true;
	}
	if (maxAz > maxAzVal+eps) {
		std::cerr << "maxAz = " << maxAz << std::endl;
		errorFlag = true;
	}
	if (errorFlag) 
		return 0;
	else
		return 1;
}

template <typename T>
int isInRangeMpiToPi(const std::vector<SG::candidate<T> >& points){
	double minZen = -0.5*M_PI;
	double maxZen = +0.5*M_PI; 
	double minAz = -1.0*M_PI; 
	double maxAz = 1.0*M_PI;
	
	return isInRange(points, minAz, maxAz, minZen, maxZen);
}

template <typename T>
int isInRangeZeroTo2Pi(const std::vector<SG::candidate<T> >& points){
	double minZen = -0.5*M_PI;
	double maxZen = +0.5*M_PI; 
	double minAz = 0.0; 
	double maxAz = 2.0*M_PI;
	
	return isInRange(points, minAz, maxAz, minZen, maxZen);
}
template <typename T>
void assertRangeMpiToPi(const std::vector<SG::candidate<T> >& points){
	int ret = isInRangeMpiToPi(points);
	if (ret == 0) throw (-1);
}

template <typename T>
void assertRangeZeroTo2Pi(const std::vector<SG::candidate<T> >& points){
	int ret = isInRangeZeroTo2Pi(points);
	if (ret == 0) throw (-1);
}

template <typename T>
int writeHammerProjection(const std::vector<SG::candidate<T> >& points, std::string outputFileName){
	assertRangeMpiToPi(points);
	if (outputFileName == ""){
		std::cerr << "writeHammerProjection: Failed to open file with no name. Please specify a outputFileName!" << std::endl;
		return -1;
	}
	std::ofstream raus(outputFileName.c_str());
	if (! raus.is_open()){
		std::cerr << "writeHammerProjection: Failed to open file with fileName " << outputFileName << " !" << std::endl;
		return -2;
	}
	if (! raus.good()){
		std::cerr << "writeHammerProjection: Cannot write to file with fileName " << outputFileName << " !" << std::endl;
		return -3;
	}
	for (unsigned int i = 0; i < points.size(); ++i){
		raus << points[i].hammerProjectionX() << " " << points[i].hammerProjectionY() << " " << points[i].proximity_ << std::endl;
		if (! raus.good()){
			std::cerr << "writeHammerProjection: Failed to write to file with fileName " << outputFileName << " !" << std::endl;
			return -4;
		}
	}

	raus.close();
	return 0;
}

template <typename T>
void writeEquatorialProjection(std::vector<SG::candidate<T> > eventsForOutput, std::string fileName)
{
	for (unsigned int i = 0; i < eventsForOutput.size(); ++i){
		eventsForOutput[i].correctEquatorialRightAscensionRangeFranzoesisch();
		eventsForOutput[i].correctEquatorialRightAscensionRangeMinusPi();
		eventsForOutput[i].invertEquatorialRightAscensionRange();
	}
	writeHammerProjection(eventsForOutput, fileName);
}

template <typename T>
void writeGalacticProjection(std::vector<SG::candidate<T> > eventsForOutput, std::string fileName)
{
	for (unsigned int i = 0; i < eventsForOutput.size(); ++i){
		eventsForOutput[i].equatorialToGalactic();
		eventsForOutput[i].correctEquatorialRightAscensionRangeFranzoesischInvers();
		eventsForOutput[i].invertEquatorialRightAscensionRange();
	}
	writeHammerProjection(eventsForOutput, fileName);
}

} // end of namespace SG

#endif /* CANDIDATE_H_ */
