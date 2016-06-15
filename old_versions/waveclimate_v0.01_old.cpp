//==============================================================
//
// waveclimate_vX.XX.cpp
//
// cpp file for the waveclimate objects
//
// objects storing wave climate PDFs or data with functions
// that can sample to get offshore waves to feed to the 
// coastline object
//
// build in functionality to be able to modify wave climates
// through time
//
// CWT has code to handle 360 degree distributions
//
// Martin Hurst, Started June 2013
//
// Also:
// Andrew Barkwith
// Chris Thomas
// Mike Ellis
//
// British Geological Survey
//
//=====================================================================

#include <math.h>
#include <vector>
#include "waveclimate_v0.01.hpp"

using namespace std;

#ifndef waveclimate_CPP
#define wavecliamte_CPP

//=====================================================================
//
// Single Wave Climate Object
//
//---------------------------------------------------------------------

	/*********************************************\
	| WaveClimate Object Initialisation Functions |
	\*********************************************/
	
void SingleWaveClimate::Initialise()
{
	cout 	<< "SingleWaveClimate.Initialise: Error, " 	
			<< "initialised an empty SingleWaveClimate object" << endl;
	exit(EXIT_FAILURE);
}

//void SingleWaveClimate::Initialise(double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight)
//{
//	double Dir_0 = OffshoreWaveDirection;
//	double Period = WavePeriod;
//	double Height_0 = OffshoreWaveHeight;
//}

////=====================================================================
////
//// PDF Wave Climate Object
////
////---------------------------------------------------------------------

//void GaussianWaveClimate::Initialise()
//{
//	cout 	<< "SingleWaveClimate.Initialise: Error, "
//			<< "initialised an empty SingleWaveClimate object" << endl;
//	exit(EXIT_FAILURE);
//}

//void GaussianWaveClimate::Initialise(OffshoreWaveHeight, OffshoreWavePeriod, OffshoreMeanWaveDirection, OffshoreStDWaveDirection)
//{
//	double Height_0 = OffshoreWaveHeight;
//	double Period = OffshoreWavePeriod;
//	double MeanDir_0 = OffshoreMeanWaveDirection;
//	double StDDir_0 = OffshoreStDWaveDirection;
//}

//double GaussianWaveClimate::GetWaveDirection()
//{
//	/* Get Wave Direction returns a wave direction sampled randomly 
//		from a normal distribution described by a mean and standard
//		deviation of offshore wave conditions

//		Random normal distribution from Box-Muller transform 
//		sqrt(-2*log((double)rand()/RAND_MAX))*sin(2*PI*((double)rand()/RAND_MAX)) */
//		
//		double OffshoreWaveDirection;
//		
//		//Get two random numbers
//		double rand1 = (double)rand()/RAND_MAX;
//		double rand2 = (double)rand()/RAND_MAX;
//		
//		//Generate random direction from distribution using random numbers
//		OffshoreWaveDirection = OffshoreMeanWaveDirection + 
//										OffshoreStDWaveDirection*sqrt(-2*log(rand1))*cos(2*PI*(rand2));
//		
//		return OffshoreWaveDirection;
//}


//=====================================================================
#endif 
