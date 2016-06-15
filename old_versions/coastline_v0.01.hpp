//==============================================================
//
// coastline.cpp
//
// cpp file for the coastline object
// coastline is a vector based object defining the location 
// and paramaters of a soft-sediment coast line.
//
// Martin Hurst, May 2013
//
//==============================================================


#include <math.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#ifndef coastline_CPP
#define coastline_CPP

class coastline {
	//data members
	//vector x;					//position in x (private)
	//vector y;					//position in y (private)
	//vector orentation;		//orientation of shoreline gradient across i-1, i+1 (private)
	//vector slope;			//shoreline gradient bewteen i, i+1 (private)
	//vector curvature;		//curvature across i-1, i+1 (private)
	//int boundary_type;		//periodic/noflux/continousflux/sinks (private)
	//node_spacing; 			//average spacing between nodes (private)
	//
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
	//Create Functions
	//
	//Calculate Morphology Functions
	//	- slope, aspect, curvature
	//
	// Shadow functions
	// - Feed transformed wave data object as a friend or subclass?
	// - Refraction? Diffraction?
	//
	// Compute sediment flux functions
	// - Feed transformed wave data object as a friend or subclass?
	// - allow various different equations
	// - 
}

class Coastline 
{
	private:
	//data members
	int NoNodes;							//Number of nodes along coastline
	vector<double> X;						//position in x (m) (private)
	vector<double> Y;						//position in y (m) (private)
	vector<double> Orientation;		//orientation/azimuth of shoreline across i-1, i+1 (private)
	vector<double> Slope;				//shoreline orientation bewteen i, i+1 (private)
	vector<double> Curvature;			//Curvatureature across i-1, i+1 (private)
	int StartBoundary;					//Boundary type periodic 1 /noflux 2 /continous flux 3 /sinks 4 (private)
	int EndBoundary;						//Boundary type periodic 1 /noflux 2 /continousflux 3 /sinks 4 (private)
	int MeanNodeSpacing;					//average spacing between nodes (m) (private)

	//Use equals operator for assignment
	// CODE GOES HERE //


	public:
	
	/********************************************
	* Coastline Object Initialisation Functions *
	********************************************/
	
	void Initialise();
	
	void Initialise(string xyfilename);
	void Initialise(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary);
	
	private:
	void CalculateMorphology();
	int CalculateMeanNodeSpacing();
		
	public:
	void SimpleDiffusion(double DiffCoeff, double TimeDelta);

};

#endif 
