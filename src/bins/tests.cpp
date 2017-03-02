#include "candidate.hpp"
#include "signalFirstInternal.h"
#include <vector>
#include "config.h"
#include <iostream>
#include "sgcu/sgcu.hpp"
#include "utils.hpp"



int main()
{
/*
	std::vector<double> radialSymmetricValues = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.07749e-07,2.15497e-07,0,0,0,0,0,0,0,5.37246e-06,1.34312e-05,2.07377e-05,2.37463e-05,2.56804e-05,3.3248e-05,3.83051e-05,3.89507e-05,4.23452e-05,4.49311e-05,4.79071e-05,4.8015e-05,4.77579e-05,4.68509e-05,5.18281e-05,5.43928e-05,5.1214e-05,4.92214e-05,5.48399e-05,5.37186e-05,4.85568e-05,4.85161e-05,4.80376e-05,5.19409e-05,5.04401e-05,4.88611e-05,4.90409e-05,4.97336e-05,5.11329e-05,4.95037e-05,4.41918e-05,4.92718e-05,4.91906e-05,4.14848e-05,4.81697e-05,4.98865e-05,4.96933e-05,4.84796e-05,3.99882e-05,4.20429e-05,4.35384e-05,4.23548e-05,4.54741e-05,3.90281e-05,4.0998e-05,3.86236e-05,3.60563e-05,4.25499e-05,4.73822e-05,3.83864e-05,4.0755e-05,4.27284e-05,3.94589e-05,4.27661e-05,4.45864e-05,3.74105e-05,3.61273e-05,4.04539e-05,3.99097e-05,3.1587e-05,2.98792e-05,3.73091e-05,3.94967e-05,3.57837e-05,3.84918e-05,3.53042e-05,3.39426e-05,3.13322e-05,3.60252e-05,3.37894e-05,3.88129e-05,3.3665e-05,3.29659e-05,3.18005e-05,3.41745e-05,3.333e-05,3.58846e-05,3.73763e-05,3.7991e-05,3.63163e-05,3.06931e-05,3.29537e-05,3.27691e-05,3.10028e-05,3.2536e-05,3.045e-05,3.11215e-05,3.37909e-05,3.86817e-05,4.00363e-05,3.55659e-05,3.03326e-05,3.11001e-05,3.11949e-05,3.30882e-05,3.26327e-05,3.2101e-05,3.00047e-05,3.24533e-05,3.48797e-05,3.275e-05,3.24318e-05,3.1082e-05,3.31558e-05,3.33553e-05,3.09101e-05,3.64484e-05,3.47051e-05,3.18233e-05,3.12785e-05,3.57783e-05,3.22348e-05,2.93712e-05,3.1569e-05,3.63135e-05,3.51834e-05,3.15935e-05,3.23034e-05,3.02364e-05,3.58203e-05,3.78221e-05,4.03162e-05,4.43344e-05,3.64328e-05,3.43329e-05,3.0992e-05,3.50147e-05,3.97663e-05,3.77024e-05,3.59366e-05,3.53734e-05,3.59649e-05,3.90778e-05,3.50427e-05,3.52403e-05,3.3356e-05,3.46703e-05,3.3993e-05,3.47836e-05,4.20587e-05,3.74238e-05,3.22348e-05,2.90939e-05,3.3592e-05,3.3803e-05,3.35123e-05,3.1532e-05,2.88217e-05,3.08672e-05,2.17585e-05,2.12126e-05,2.7507e-05,2.48986e-05,2.21038e-05,3.29511e-05,3.52103e-05,3.24927e-05,2.73996e-05,2.48467e-05,2.06696e-05,2.41374e-05,2.08036e-05,1.64312e-05,2.55445e-05,1.702e-05,1.15222e-05,1.84617e-05,1.98126e-05,1.65042e-05,2.35454e-05,2.57878e-05,2.48668e-05,2.47563e-05,1.75826e-05,1.6287e-05,1.2478e-05,6.18908e-06,8.14352e-06,2.38041e-05,2.57878e-05,7.73635e-05};
	sgcu::saveVectorToFile(radialSymmetricValues, "vec0");

	// smoothing away the last two bins on both ends
	radialSymmetricValues[0] = radialSymmetricValues[2];
	radialSymmetricValues[1] = radialSymmetricValues[2];
	radialSymmetricValues[radialSymmetricValues.size()-1] = radialSymmetricValues[radialSymmetricValues.size()-3];
	radialSymmetricValues[radialSymmetricValues.size()-2] = radialSymmetricValues[radialSymmetricValues.size()-3];
	sgcu::saveVectorToFile(radialSymmetricValues, "vec1");
	

	sgcu::smoothMedian(radialSymmetricValues);
	sgcu::saveVectorToFile(radialSymmetricValues, "vec2");

	//sgcu::smoothMedian5(radialSymmetricValues);
	sgcu::saveVectorToFile(radialSymmetricValues, "vec3");


	unsigned int numSmoothNoCloseToZeros = 20;	// for IceCube, large dataset
	// cout << "debug: smoothRadialSymmetricValues: number of smoothing operations : " << numSmoothNoZeros << endl;
	for (unsigned int i = 0; i < numSmoothNoCloseToZeros; ++i)
	{
		smoothNoCloseToZerosKillingEnds(radialSymmetricValues);
	}
	sgcu::saveVectorToFile(radialSymmetricValues, "vec4");


	unsigned int numSmoothNoZeros = 60;	// for ANTARES
	numSmoothNoZeros = 15;	// for IceCube, IC40
	numSmoothNoZeros = 2;	// for IceCube, large dataset
	// cout << "debug: smoothRadialSymmetricValues: number of smoothing operations : " << numSmoothNoZeros << endl;
	for (unsigned int i = 0; i < numSmoothNoZeros; ++i)
	{
		smoothNoZerosKillingEnds(radialSymmetricValues);
	}
	sgcu::saveVectorToFile(radialSymmetricValues, "vec5");

*/

//
//	time_t seed = 123;
//	if (useRandomSeed)
//	{
//		seed = time(NULL);
//	}
//	srand(seed);
///*
//	std::string commando = "mkdir -p ";
//	std::string outputPath = "output/scenarioIC40"+boost::lexical_cast<std::string>(seed);
//	std::cout << outputPath << std::endl;
//	system(std::string(commando + outputPath).c_str());
//*/
//
//	// set up an object to handle all computations with a simple interface to this module
//	SG::signalFirstInternal intern;
//	intern.init(fine);
//
//	const double oneDeg = M_PI/180.0;
//
//	std::vector<SG::candidate_t> otherCandidates;
//	//readCandidatesFromFile(candidateFile, allCandidates);
////	readCandidatesFromFile("input/IC40/IceCube40_converted.txt", otherCandidates);
//	//readCandidatesFromFile("input/IC40/IC40With80Events.txt", otherCandidates);
//	readCandidatesFromFile("input/ICFull/IceCubeDatasetLocalCoords.txt", otherCandidates);
//
//	std::vector<SG::candidate_t> otherCandidatesOrg = otherCandidates;
//	// use randomized events
//	intern.getRandomEventsFromEvents(otherCandidatesOrg, otherCandidates);
//
//
//
//
//	// add source events
////	intern.addSourceGaussian(otherCandidates, 260.0*oneDeg, 20.0*oneDeg, 10.0*oneDeg, 50);
////	intern.dumpEvents(otherCandidates, "input/ICFull/with100Events20x20.txt");
//	intern.addSource(otherCandidates, 200.0*oneDeg, 35.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 300);
//	intern.dumpEvents(otherCandidates, "input/ICFull/with300Events20x20.txt");
//
//	// get a new expectation for these events
//	intern.computeNormalizationsFromDataStep1(otherCandidates, fine, false);
//
//	// Calculate the scenario
//	intern.setOutputFolder("output/scenarioICFull/");
////	intern.setOutputFolder(outputPath + "/");
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(otherCandidates, fine, configDoDebugOutput);
//
//return 0;
////	/*
//	std::vector<SG::candidate_t> allCandidates;
//	readCandidatesFromFile(candidateFile, allCandidates);
//	std::vector<SG::candidate_t> allCandidatesOrg = allCandidates;
//
//	// Calculate scenario with a diffuse hypothesis and a perfectly matching diffuse source
//	intern.setOutputFolder("output/scenarioLargeFullHit2/");
//	// allCandidates = allCandidatesOrg;
//	// use randomized events
//	intern.getRandomEventsFromEvents(allCandidatesOrg, allCandidates);
//	// add source events
//	intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 80);
//	//intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 50);
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(allCandidates, fine, configDoDebugOutput);
//return 0;
//	// Calculate scenario with a diffuse hypothesis and a perfectly matching diffuse source
//	intern.setOutputFolder("output/scenarioLargeFullHit/");
//	// allCandidates = allCandidatesOrg;
//	// use randomized events
//	intern.getRandomEventsFromEvents(allCandidatesOrg, allCandidates);
//	// add source events
//	// intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 80);
//	intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 50);
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(allCandidates, fine, configDoDebugOutput);
//
//	// Calculate scenario with a diffuse hypothesis and a contained source half the size
//	intern.setOutputFolder("output/scenarioLargeHalfHit/");
//	// use randomized events
//	intern.getRandomEventsFromEvents(allCandidatesOrg, allCandidates);
//	// add source events
//	// intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 5.0*oneDeg, 10.0*oneDeg, 60);
//	intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 5.0*oneDeg, 10.0*oneDeg, 30);
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(allCandidates, fine, configDoDebugOutput);
//
//	// Calculate scenario with a diffuse hypothesis and a source of the estimated size, contained just half
//	intern.setOutputFolder("output/scenarioLargeHalfHitHalfMiss/");
//	// use randomized events
//	intern.getRandomEventsFromEvents(allCandidatesOrg, allCandidates);
//	// add source events
//	// intern.addSource(allCandidates, 270.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 80);
//	intern.addSource(allCandidates, 270.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 50);
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(allCandidates, fine, configDoDebugOutput);
//
//	// Calculate scenario with a diffuse hypothesis and a three contained weak point sources
//	intern.setOutputFolder("output/scenarioLarge3PSHit/");
//	// use randomized events
//	intern.getRandomEventsFromEvents(allCandidatesOrg, allCandidates);
//	// add source events
//	//intern.addSource(allCandidates, 250.0*oneDeg, 10.0*oneDeg, 1.5*oneDeg, 30);
//	//intern.addSource(allCandidates, 260.0*oneDeg, 20.0*oneDeg, 1.5*oneDeg, 20);
//	//intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 1.0*oneDeg, 20);
//	intern.addSource(allCandidates, 250.0*oneDeg, 10.0*oneDeg, 1.5*oneDeg, 15);
//	intern.addSource(allCandidates, 260.0*oneDeg, 20.0*oneDeg, 1.5*oneDeg, 10);
//	intern.addSource(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 1.0*oneDeg, 10);
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(allCandidates, fine, configDoDebugOutput);
//
//
//	// Calculate scenario with a diffuse hypothesis and a three contained weak point sources
//	intern.setOutputFolder("output/scenarioLargeGaussian/");
//	// use randomized events
//	intern.getRandomEventsFromEvents(allCandidatesOrg, allCandidates);
//	// add source events
//	intern.addSourceGaussian(allCandidates, 260.0*oneDeg, 10.0*oneDeg, 10.0*oneDeg, 40);
//	// Evaluate the clustering/density of the candidates here
//	intern.evaluateEventsToSkymap(allCandidates, fine, configDoDebugOutput);
//
//	//	*/

	return 0;
}
