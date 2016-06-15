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
		
	/*************************************\
	| SingleWave Initialisation Functions |
	\*************************************/
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
	
	double Get_WavePeriod();			//get Wave Period	
	double Get_WaveDirection();		//get Wave Direction
	double Get_WaveHeight();			//get Wave Height
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

	double Get_WavePeriod();			//get Wave Period	
	double Get_WaveDirection();		//get Wave Direction
	double Get_WaveHeight();			//get Wave Height
	
};

class GaussianWaveClimate 
{
	
	//BimodalWaveClimate will need to be a friend since it uses two GaussianWaveClimates
	friend class BimodalWaveClimate;
	
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
	
	//Return functions to allow BimodalWaveClimate access to private members
	//Could be replaced subsequently by making them "friends"!?
	
};

class RealWaveClimate 
{
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

class UAWaveClimate 
{
	private:
	
	//data members
	double Period_Mean;
	double Period_StD;
	double Height_Mean;
	double Height_StD;
	double U;
	double A;
	double CoastTrend;
	
	public:
	
	/*****************************************************\
	|    UAWaveClimate Object Initialisation Functions    |
	\*****************************************************/
	
	UAWaveClimate()	
	{ 
		Initialise(); 
	}
	UAWaveClimate(	double input_U, double input_A, double Trend, 
						double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
						double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight)
	{
		Initialise(	input_U, input_A,Trend, 
						OffshoreMeanWavePeriod, OffshoreStDWavePeriod, 
						OffshoreMeanWaveHeight, OffshoreStDWaveHeight);
	}
	void Initialise();
	void Initialise( double input_U, double input_A, double Trend, 
							double OffshoreMeanWavePeriod, double OffshoreStDWavePeriod,
							double OffshoreMeanWaveHeight, double OffshoreStDWaveHeight);
	
	Wave Get_Wave();											
};


class BimodalWaveClimate 
{
	private:
	
	//data members
	GaussianWaveClimate WaveMode1;
	GaussianWaveClimate WaveMode2;
	
	public:
	
	/*****************************************************\
	| BimodalWaveClimate Object Initialisation Functions |
	\*****************************************************/
	
	BimodalWaveClimate()	
	{ 
		Initialise(); 
	}
	BimodalWaveClimate(double FractionWaveDirection1,
							double OffshoreMeanWavePeriod1, double OffshoreStDWavePeriod1,
							double OffshoreMeanWaveDirection1, double OffshoreStDWaveDirection1,
							double OffshoreMeanWaveHeight1, double OffshoreStDWaveHeight1,
							double OffshoreMeanWavePeriod2, double OffshoreStDWavePeriod2,
							double OffshoreMeanWaveDirection2, double OffshoreStDWaveDirection2,
							double OffshoreMeanWaveHeight2, double OffshoreStDWaveHeight2)
	{
		Initialise(	FractionWaveDirection1,
						OffshoreMeanWavePeriod1, OffshoreStDWavePeriod1,
						OffshoreMeanWaveDirection1, OffshoreStDWaveDirection1,
						OffshoreMeanWaveHeight1, OffshoreStDWaveHeight1,
						OffshoreMeanWavePeriod2, OffshoreStDWavePeriod2,
						OffshoreMeanWaveDirection2, OffshoreStDWaveDirection2,
						OffshoreMeanWaveHeight2, OffshoreStDWaveHeight2);
	}
	void Initialise();
	void Initialise(double FractionWaveDirection1,
							double OffshoreMeanWavePeriod1, double OffshoreStDWavePeriod1,
							double OffshoreMeanWaveDirection1, double OffshoreStDWaveDirection1,
							double OffshoreMeanWaveHeight1, double OffshoreStDWaveHeight1,
							double OffshoreMeanWavePeriod2, double OffshoreStDWavePeriod2,
							double OffshoreMeanWaveDirection2, double OffshoreStDWaveDirection2,
							double OffshoreMeanWaveHeight2, double OffshoreStDWaveHeight2);											
	
	Wave Get_Wave();
	
};
#endif
