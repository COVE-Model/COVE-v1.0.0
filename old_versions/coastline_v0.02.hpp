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
	vector<double> Orientation;			//orientation/azimuth of shoreline across i-1, i+1 (private)
	vector<double> FluxOrientation;		//shoreline orientation bewteen i, i+1 (private)
	vector<double> Curvature;				//Curvatureature across i-1, i+1 (private)
	vector<double> Area;						//Area of coastal cell (varies with curvature)
	vector<double> CellWidth;				//Width of individual cells perpendicular to shoreline trend (for now)
	vector<double> BreakingWaveHeight;	//Heights of breaking waves
	vector<double>	BreakingWaveAngle;	//Angle of breaking waves
	vector<double> Shadows;					//Shadow Zone Vector
	vector<double> LongshoreFlux;			//Volume of sediment transported alongshore in m3/day
	vector<double> VolumeChange;			//Volume change in each cell during a particular timestep (m3/day)
	vector<double> PositionChange;		//Magnitude of change in shoreline position
	vector<double> EmptyVector;
	int StartBoundary;						//Boundary type periodic 1 /noflux 2 /continous flux 3 /sinks 4 (private)
	int EndBoundary;							//Boundary type periodic 1 /noflux 2 /continousflux 3 /sinks 4 (private)
	int MeanNodeSpacing;						//average spacing between nodes (m) (private)
	double Trend;
	double ClosureDepth;
	double LostFluxFraction;
	
	//Use equals operator for assignment
	// CODE GOES HERE //

	void CalculateMorphology();
	void CalculateMeanNodeSpacing();
	void GetArea();
	
	public:
	
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
	
	void ReadCoast(string InputFileName, int Time);
	void WriteCoast(string OutputFileName, int Time);
	
	/****************************************\
	| Get Function to return private members |
	\****************************************/

	int get_NoNodes() const							 {return NoNodes;}			//get Number of nodes along coastline
	vector<double> get_X() const					 {return X;}					//get X vector position in x (m)
	vector<double> get_Y() const					 {return Y;}					//get Y vectir position in y (m) 
	vector<double> get_Orientation() const 	 {return Orientation;}		//get orientation of shoreline between i-1 and i+1
	vector<double> get_FluxOrientation() const {return FluxOrientation;}	//get shoreline orientation bewteen i and i+1
	vector<double> get_Curvature() const		 {return Curvature;}			//get Curvature across i-1, i+1
	int get_StartBoundary() const					 {return StartBoundary;}	//get Boundary type start
	int get_EndBoundary() const					 {return EndBoundary;}		//get Boundary type end
	int get_MeanNodeSpacing() const				 {return MeanNodeSpacing;}	//get average spacing between nodes (m)
		
	void SimpleDiffusion(double DiffCoeff, double TimeDelta);
	void CERCDiffusion(double &TimeDelta, double MaxTimeStep, double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight);
	void TransformWaves(double T, double H_0, double Theta_0);		//Transform wave of period T and offshore height H_0
	void GetShadows(double OffshoreWaveDirection);
	
	void KamphiusDiffusion();
};

#endif 
