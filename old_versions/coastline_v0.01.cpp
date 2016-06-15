//========================================================================
//
// coastline.cpp
//
// cpp file for the coastline object
// coastline is a vector based object defining the location 
// and paramaters of a soft-sediment coast line.
//
// Martin Hurst, Started May 2013
//
// Also:
// Andrew Barkwith
// Chris Thomas
// Mike Ellis
// 
// British Geological Survey
//
//
//=========================================================================
//
//	18/6/13 	Initialise() functions added: initialises coastline from XY file
//				CalculateMorphology() added: gets Orientation, Slope, Curvatureature
//
//=========================================================================


#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "coastline_v0.01.hpp"

using namespace std;

#ifndef coastline_CPP
#define coastline_CPP
	
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
	
	void Initialise()
	{
		cout << "Coastline.Initialise: Error, initialised an empty Coastline object" << endl;
		exit(EXIT_FAILURE);
	}
	
	void Initialise(string xyfilename)
	{
		/*	Read coastline.XY text file which has the following structure (| is delimeter)
		
						x		  | 		y
				StartBoundary | EndBoundary
						x0		  |		y0
						x1		  |		y1
						...	  |		...
						Xn		  |		yn		*/
		
		cout << "Coastline: Initialising Coastline from XY file: " << xyfilename << endl;
		
		//Declare temporary variables for reading data
		double tempX, tempY;
		char temp[5];
		NoNodes = 0;
		
		//Open the file in a filestream and read data into X and Y
		ifstream ReadCoastFile;
		ReadCoastFile.open(xyfilename.c_str());
		if (ReadCoastFile.is_open())
		{
			//handle header lines
			if (ReadCoastFile.good()) ReadCoastFile >> temp >> temp;
			if (ReadCoastFile.good()) ReadCoastFile >> StartBoundary >> EndBoundary;
						
			//read XY data
			while (ReadCoastFile.good())
			{
				ReadCoastFile >> tempX >> tempY;
				X.push_back(tempX);
				Y.push_back(tempY);
				NoNodes += 1;
			}
		}
		//Populate empty vectors
		vector<double> EmptyVector(NoNodes, -9999);
		Orientation = EmptyVector;
		Slope = EmptyVector;
		Curvature = EmptyVector;
		
		//Calculate Orientation, Slope and Curvature
		CalculateMorphology();
		CalculateMeanNodeSpacing();
	}
	void Initialise(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary)
	{
	
		/*	Initialises the coast as a straight segment with low amplitude noise
			Coast has a fixed length and an orientation/trend. Coordinates of the
			first node in the arrays is [0][0].
			Boundary conditions must be specified */
		
		cout 	<< "Coastline: Initialising Coastline as straight segment"
				<< "\n\t MeanNodeSpacing = " << MeanNodeSpacing
				<< "\n\t CoastLength = " << CoastLength
				<< "\n\t Trend = " << Trend 
				<< "\n\t Boundary Conditions = " << StartBoundary << ", " << EndBoundary << endl;
		
		//declare temporary variables for setting up coastline object
		double tempX, tempY;
		NoNodes = 0;
		double Distance = 0;
		
		//initialise coast
		while (Distance<CoastLength)
		{
			if (Distance==0)
			{
				tempX = 0; tempY = 0;
				X.push_back(tempX); 
				Y.push_back(tempY);
				NoNodes += 1;
			}
			else
			{
				//sort out random addition
				tempX = tempX + MeanNodeSpacing*cos(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
				tempY = tempY - MeanNodeSpacing*sin(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
				NoNodes += 1;
			}
			Distance += MeanNodeSpacing;
		}
		
		//Populate empty vectors
		vector<double> EmptyVector(NoNodes, -9999);
		Orientation = EmptyVector;
		Slope = EmptyVector;
		Curvature = EmptyVector;
		
		//Calculate Orientation, Slope and Curvature
		CalculateMorphology();
		
	}
	
	private:
	void CalculateMorphology()
	{
		/************************************************************
		Calculates Orientation, Slope and Curvatureature members
		of Coastline object from X and Y positions
		
		Orientation is the slope over total distance between i-1 and i+1
		Slope (orientation between i and i+1) is dy/dx
		Curvature is change in Slope over total distance between i-1 and i+1

		first handle periodic boundary
		
		!! ** WILL NEED TO ADD LOGIC FOR DIFFERENT BOUNDARIES** !!
		
		************************************************************/
		
		//declare some temporary variables
		double dX, dY;
		for (int i=0; i<NoNodes; ++i) 
		{
			//Handle StartBoundary
			if (i==0) 
			{
				if (StartBoundary == 1) 
				{
					Slope[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]);
				}
				else 
				{
					cout << "Coastline.CalculateMorphology: Error, non-periodic StartBoundary" << endl;
					exit(EXIT_FAILURE);
				}
			}
			//Deal with EndBoundary
			else if (i==NoNodes-1)
				if (EndBoundary == 1)
				{
					dX = X[i]-X[i-1]+X[1]-X[0];
					dY = Y[i]-Y[i-1]+Y[1]-Y[0];
					Slope[i] = Slope[0];
					Curvature[i] = (Slope[i]-Slope[i-1])/sqrt(dX*dX + dY*dY);
					Curvature[0] = Curvature[i];
					if (dX > 0) Orientation[i] = (M_PI*0.5 - atan((dY/dX)));
					else if (dX < 0) Orientation[i] = (M_PI*1.5 - atan((dY/dX)));
					Orientation[0] = Orientation[i];
				}
				else
				{
					cout << "Coastline.CalculateMorphology: Error, non-periodic EndBoundary" << endl;
					exit(EXIT_FAILURE);
				}
			//Finally do all the bits inbetween
			else 
			{
				dX = (X[i+1]-X[i-1]);
				dY = (Y[i+1]-Y[i-1]);
				Slope[i] = (Y[i+1]-Y[i])/(X[i+1]-X[i]);
				Curvature[i] = (Slope[i]-Slope[i-1])/sqrt(dX*dX + dY*dY);
				if (dX > 0) Orientation[i] = (M_PI*0.5 - atan((dY/dX)));
				else if (dX < 0) Orientation[i] = (M_PI*1.5 - atan((dY/dX)));
			}
		}
	}

	int CalculateMeanNodeSpacing()
	{
		//Returns an integer of the mean node spacing from X and Y data
		vector<double> Distance(NoNodes);
		double dx, dy, MeanDistanceDouble;
		double TotalDistance = 0;
		int MeanDistanceInt;
		
		for (int i=0; i<NoNodes; ++i)
		{
			//For the last cell Distance is the same as for cell i=0
			if (i==NoNodes-1) Distance[i] = Distance[0];
			
			//Calculate distance between i and i+1
			dx = fabs(X[i+1]-X[i]);
			dy = fabs(Y[i+1]-Y[i]);
			Distance[i] = sqrt(dx*dx + dy*dy);
			TotalDistance += Distance[i];
		}

		//Calculate MeanNodeSpacing
		MeanDistanceDouble = TotalDistance/NoNodes;
		
		//Convert to int
		MeanDistanceInt = (int) MeanDistanceDouble;
		return MeanDistanceInt;
	}
	
	
	
	public:
	void SimpleDiffusion(double DiffCoeff, double TimeDelta)
	{
		//do diffusion
		//Diffuse coastline as a function of curvature
		//Coastline adjusts perpendicular to orientation
		//Calls CalculateMorphology() at the end to update
		
		//temp variables
		double alpha, Change;
		
		for (int i=0; i<NoNodes; ++i) 
		{
			Change = DiffCoeff*Curvature[i]*TimeDelta;
			
			//update X and Y
			if (Orientation[i] < M_PI/2) 
			{
				alpha = Orientation[i];
				X[i] += Change*cos(alpha);
				Y[i] -= Change*sin(alpha);
			}
			else if (Orientation[i] < M_PI)
			{
				alpha = M_PI-Orientation[i];
				X[i] -= Change*cos(alpha);
				Y[i] -= Change*sin(alpha);
			}
			else if (Orientation[i] < 3*M_PI/2)
			{
				alpha = Orientation[i]-M_PI;
				X[i] -= Change*cos(alpha);
				Y[i] += Change*sin(alpha);
			}
			else
			{
				alpha = Orientation[i]-M_PI;
				X[i] += Change*cos(alpha);
				Y[i] += Change*sin(alpha);
			}
		}
		
		//Update Coastal Morphology
		CalculateMorphology();
	}
	
	//will need shelf slope/shoreface slope?
	//M_PIn platform/shelf inflection point at mean sea-level
	//keep track of beach width?
	//
	//
	//function members
	//get_aspect				(private)
	//get_Curvatureature			(private)
	//reinterpolate			(private)
	//
	//
	//Create Functions
	//
	//Calculate Morphology Functions
	//	- slope, aspect, Curvatureature
	//
	// Shadow functions
	// - Feed transformed wave data object as a friend or subclass?
	// - Refraction? Diffraction?
	//
	// Compute sediment flux functions
	// - Feed transformed wave data object as a friend or subclass?
	// - allow various different equations
	// - 
};

//int main()
//{
//	string filename = "coastline0.txt";
//	Coastline Coast;
//	Coast.Initialise(filename);
//}

#endif 
