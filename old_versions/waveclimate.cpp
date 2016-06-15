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

#endif 
