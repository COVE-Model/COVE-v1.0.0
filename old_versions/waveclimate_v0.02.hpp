/*==============================================================

 waveclimate_vX.XX.cpp

 cpp file for the waveclimate objects

 objects storing wave climate PDFs or data with functions
 that can sample to get offshore waves to feed to the 
 coastline object

 build in functionality to be able to modify wave climates
 through time

 CWT has code to handle 360 degree distributions

 Martin Hurst, Started June 2013

 Also:
 Andrew Barkwith
 Chris Thomas
 Mike Ellis

 British Geological Survey

==============================================================

 Updates detailed in the waveclimate_vX.cpp file

==============================================================
*/

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <cstring>

using namespace std;

#ifndef waveclimate_HPP
#define waveclimate_HPP


class Wave
{
	private:
	
	//data members
	double Dir;
	double Period;
	double Height;
	
	public:
		
	/********************************************\
	| SingleWaveClimate Initialisation Functions |
	\********************************************/
	Wave()	
	{ 
		Initialise(); 
	}
	Wave(double WvPrd, double WvHgt, double WvDir)
	{
		Initialise( WvPrd, WvHgt, WvDir );
	}
	
	void Initialise();
	void Initialise( double WvPrd,  double WvHgt, double WvDir );
	
	void AssignWaveDirection( double WvDir );
	void AssignWaveHeight( double WvHgt );
	void AssignWavePeriod( double WvPrd );
	
	double Get_WavePeriod();		//get Wave Period	
	double Get_WaveDirection();		//get Wave Direction
	double Get_WaveHeight();		//get Wave Height
};


class SingleWaveClimate
{
	private:

	//data members
	double Dir;
	double Period;
	double Height;
	
	public:
	
	/********************************************\
	| SingleWaveClimate Initialisation Functions |
	\********************************************/
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

	/*****************************************\
	| Get Functions to return private members |
	\*****************************************/

	double Get_WavePeriod();		//get Wave Period	
	double Get_WaveDirection();		//get Wave Direction
	double Get_WaveHeight();		//get Wave Height
	
};

class GaussianWaveClimate {
	
	private:
	
	//data members
	double Period_Mean;
	double Period_StD;
	double Dir_Mean;
	double Dir_StD;
	double Height_Mean;
	double Height_StD;
	
	public:
	
	/*****************************************************\
	| GaussianWaveClimate Object Initialisation Functions |
	\*****************************************************/
	
	GaussianWaveClimate()	
	{ 
		Initialise(); 
	}
	GaussianWaveClimate(double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveDirection, double OffshoreStDWaveDirection,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight)
	{
		Initialise(	OffshoreMeanWavePeriod, OffshoreStDWavePeriod,
						OffshoreMeanWaveDirection, OffshoreStDWaveDirection,
						OffshoreMeanWaveHeight, OffshoreStDWaveHeight);
	}
	void Initialise();
	void Initialise(	double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveDirection, double OffshoreStDWaveDirection,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight);											
	
	Wave Get_Wave();
	
};

class RealWaveClimate {
	
	private:
	
	//data members
	vector<double> Period;
	vector<double> Dir;
	vector<double> Height;
	int NoWaves, WaveCounter;
	
	public:
	
	/*****************************************************\
	|   RealWaveClimate Object Initialisation Functions   |
	\*****************************************************/
	
	RealWaveClimate()	
	{ 
		Initialise(); 
	}
	RealWaveClimate(string WaveFileName)
	{
		Initialise(	WaveFileName);
	}
	void Initialise();
	void Initialise(string WaveFileName);											
	
	Wave Get_Wave();
	
};

#endif
