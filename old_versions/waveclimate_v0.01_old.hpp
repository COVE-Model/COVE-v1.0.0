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
//==============================================================

#include <math.h>
#include <vector>

using namespace std;

#ifndef waveclimate_HPP
#define wavecliamte_HPP

class SingleWaveClimate {
	
	private:

	//data members
	double Dir_0;
	double Period;
	double Height_0;
	
	public:
	
	SingleWaveClimate()	
	{ 
		Initialise(); 
	}
	SingleWaveClimate(double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight) 
	{ 
		Initialise(OffshoreWavePeriod, OffshoreWaveDirection, OffshoreWaveHeight); 
	}
	void Initialise();	
	void Initialise(double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight);
	
	/****************************************\
	| Get Function to return private members |
	\****************************************/

	//double get_OffshoreWavePeriod() const			{return Period;}		//get Wave Period
	//double get_OffshoreWaveDirection() const		{return Dir_0;}		//get Wave Direction
	//double get_OffshoreWaveHeight() const			{return Height_0;}	//get Wave Height
	
};

//class GaussianWaveClimate {
//	
//	//data members
//	double Height_0;
//	double Period;
//	double MeanDir_0;
//	double StDDir_0;
//};

#endif

