//==============================================================
//
// coastline.hpp
//
// cpp file for the coastline object
// coastline is a vector based object defining the location
// and paramaters of a soft-sediment coast line.
//
// Martin Hurst, started May 2013
//
// Also:
// Andrew Barkwith
// Chris Thomas
// Mike Ellis
//
// British Geological Survey
//
//==============================================================
//
// Updates detailed in the coastline_vX.cpp file
//
//==============================================================

#include <cmath>
#include <cstdlib>
#include <vector>
#include <cstring>
#include "waveclimate_v0.02.hpp"

using namespace std;

#ifndef coastline_HPP
#define coastline_HPP

	//will need shelf slope/shoreface slope?
	//pin platform/shelf inflection point at mean sea-level
	//keep track of beach width?
	//
	//
	//function members
	//get_aspect				(private)
	//get_curvature			(private)
	//reinterpolate			(private)
	//
	//
	// Shadow functions
	// - Feed transformed wave data object as a friend or subclass?
	// - Refraction? Diffraction?
	//
	// Compute sediment flux functions
	// - Feed transformed wave data object as a friend or subclass?
	// - allow various different equations
	// -

class Coastline
{
	private:
	//data members
	int NoNodes;								//Number of nodes along coastline
	vector<double> X;							//position in x (m) (private)
	vector<double> Y;							//position in y (m) (private)
	vector<double> Orientation2Cell;		//shoreline orientation bewteen i, i+1 (private)
	vector<double> Orientation3Cell;		//orientation/azimuth of shoreline across i-1, i+1 (private)
	vector<double> CellWidth;				//Width of individual cells parallel to orientation
	vector<double> BreakingWaveHeight;		//Heights of breaking waves
	vector<double>	BreakingWaveAngle;		//Angle of breaking waves
	vector<double> Shadows;					//Shadow Zone Vector
	vector<double> ShadowZoneWaveDirection;
	vector<double> ShadowZoneWaveHeight;
	vector<double> LongshoreFlux;			//Volume of sediment transported alongshore in m3/day
	vector<double> VolumeChange;			//Volume change in each cell during a particular timestep (m3/day)
	vector<double> PositionChange;			//Magnitude of change in shoreline position
	vector<double> EmptyVector;				//For reseting vectors
	
	vector<int> NodesToRemove;				//For removing nodes that are getting squashed
	
	int StartBoundary;						//Boundary type periodic 1 /noflux 2 /continous flux 3 /sinks 4 (private)
	int EndBoundary;						//Boundary type periodic 1 /noflux 2 /continousflux 3 /sinks 4 (private)
	int MeanNodeSpacing;					//average spacing between nodes (m) (private)
	int DesiredNodeSpacing;					//for adding and removing nodes?
	int ShadowFlag;
	double Trend;							//redundant?
	
	double OffshoreWavePeriod, OffshoreWaveHeight, OffshoreWaveDirection;
	double ClosureDepth;					// } This needs to 
	double LostFluxFraction;				// } go somewhere else?
	double ShorefaceSlope;					// }
	
	//Use equals operator for assignment
	// CODE GOES HERE //

	void CalculateMorphology();
	void CheckNodeSpacing();
	void CalculateMeanNodeSpacing();
	void IntersectionAnalysis();
	void TransformWave(int i);
	void TransformWaves();
	void GetShadows();
	void RefractDiffractShadowZone();

	public:

	//temporarily public
	
	
	/*******************************************\
	| Coastline Object Initialisation Functions |
	\*******************************************/
	Coastline()
	{
		Initialise();
	}
	Coastline(string xyfilename)
	{
		Initialise(xyfilename);
	}
	Coastline(string xyfilename, int StartTime)
	{
		Initialise(xyfilename, StartTime);
	}
	Coastline(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary)
	{
		Initialise(MeanNodeSpacing, CoastLength, Trend, StartBoundary, EndBoundary);
	}

	void Initialise();
	void Initialise(string xyfilename);
	void Initialise(string xyfilename, int StartTime);
	void Initialise(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary);

	void ReadCoast(string InputFileName, double Time);
	void WriteCoast(string OutputFileName, double Time);

	/****************************************\
	| Get Function to return private members |
	\****************************************/

	int get_NoNodes() const						{return NoNodes;}			//get Number of nodes along coastline
	vector<double> get_X() const				{return X;}					//get X vector position in x (m)
	vector<double> get_Y() const				{return Y;}					//get Y vectir position in y (m)
	//vector<double> get_Orientation() const 	 	{return Orientation;}		//get orientation of shoreline between i-1 and i+1
	//vector<double> get_FluxOrientation() const 	{return FluxOrientation;}	//get shoreline orientation bewteen i and i+1
	//vector<double> get_Curvature() const		{return Curvature;}			//get Curvature across i-1, i+1
	vector<double> get_Shadows() const			{return Shadows;}			//get Shadows	
	int get_StartBoundary() const				{return StartBoundary;}		//get Boundary type start
	int get_EndBoundary() const					{return EndBoundary;}		//get Boundary type end
	int get_MeanNodeSpacing() const				{return MeanNodeSpacing;}	//get average spacing between nodes (m)

	
	/*Flux Type 	0 = Simple Diffusion
					1 = CERC equation (default)
					2 = Kamphius equation
					3 = Bailard equation
					4 = Deigaard equation
	*/
	
	void TransportSediment(double &TimeDelta, double MaxTimeStep, Wave TheWave, int FluxType=1);
	void CalculateFlux(int i, int FluxType);
	void SimpleDiffusion(double DiffCoeff, double TimeDelta);
};

#endif
