//========================================================================
//
// coastline.cpp
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
//
//=========================================================================
//
// 4/4/14 v0.28 Added cubic and quadratic equation solver functions for finding the
//				change in position of the coast, cubic solution is a new approach to
//				handling triangular and complex cell geometries
//
// 12/3/14 v0.27 Modified handling of triangular cells so that cells only advance/retreat
//					over their effective shoreface depth (i.e. the depth to the point of
//					the triangle). The triangle retreats to become a trapezium in the case
//					of coastal recession, and advances as a triangle (overriding adjacent cells
//					in the case of coastal aggredation.
//
// 7/3/14  v0.26 Many changes and debugging to get the model working
//				 Modified treatment of quadratic solved to get position change 
//				 to account for the wedge-shaped cell geometry
//
// 20/2/14 v0.22 New Shadow finding algorithm without problems for hooked coasts
//
// 19/2/14 v0.21 Bugs ironed out in wave transformation code
//
// 19/2/14 v0.20 Multiple back ups due to network down time...
//
// 19/2/14 v0.19 Multiple back ups due to network down time...
//
// 13/2/14 v0.18 Improved approach to diffraction in the lee of a headland
//				following Goda et al 1978 and Kraus 1984. 
//
// 13/2/14 v0.17 Longshore variation in wave height affects sediment flux following
//				Ozasa and Brampton (1980). NOT WORKING YET
//
// 13/2/14 v0.16 Modifying the code to allow fixed sediment inputs at the up and downdrift
//				boundaries. A function of waves i.e. the amount
//				that waves can transport in given the orientation of the boundary
//
// 12/2/14 v0.15 Reverted to v0.12 after breaking v0.14!
//				Transition form alpha_0 < 45 to alpha_0 > 45 and vica versa
//				now uses alpha_0 = 45 to maximise transport aka Ashton et al 2006.
//				Fixed random seed generation to generate truly random waves at each
//				time step, this has smoothed out much of the modelling issues for
//				high angle waves.
//				
// 11/2/14 v0.14 Reverted to v0.12 after breaking v0.13!
//				Transition form alpha_0 < 45 to alpha_0 > 45 and vica versa
//				now uses alpha_0 = 45 to maximise transport aka Ashton et al 2006.
//
// 8/2/14 v0.13 Transition form alpha_0 < 45 to alpha_0 > 45 and vica versa
//				now uses alpha_0 = 45 to maximise transport aka Ashton et al 2006.
//
// 7/2/14 v0.12 Wave refraction/diffraction finalised and debugged
//				Code clear up
//
// 4/2/14 v0.11 Bug fixing v0.10 had a memory leak that would not allow it to add
//				new nodes to the vectors. Reverted to v0.09 and updated
//				Shadow Zone now has 5 int codes (see description)
//
// 3/2/14 v0.10 Added simple rules for wave refraction/diffraction in the shadow zone
//				Offshore angle modified by 1.5x(omega-theta)
//				Offshore height modified by cos(omega-theta)
//				Offshore Waves > 10 degrees from coast trend only refracted
//
// 29/1/13 v_0.09 Seperate out sediment flux calculations and moving of sediment
//					to allow different sediment flux laws to be used
//					New function to calculate longshore flux at a particular node i
//					with 5 different flux laws in it
//
// 28/1/13 v0.08 Added function to check for coastline intersecting itself and flag
//
// 21/1/13 v_0.07 Shadows algorithm debugged.
//
// 15/1/13 	v_0.06 Shadows Algorithm improved to better handle complex coastlines
//
// 27/11/13 v 0.05 Non-rectilinear cell shapes added and node migration perpendicular
//					to the orientation rather than the trend
//
// 22/11/13 v 0.04 CERC diffusion added, rectilinear style
//
//	7/8/13	WriteCoast() and ReadCoast() functions added
//
// 5/8/13	Intialise() functions implemented and tested
//				CaluclateMeanNodeSpacing() added: returns mean node spacing int
//				SimpleDiffusion() added: performs simple diffusion no waves
//
//	18/6/13 	Initialise() functions added: initialises coastline from XY file
//				CalculateMorphology() added: gets Orientation, Slope, Curvature
//
//=========================================================================

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "coastline_v0.28.hpp"
#include "waveclimate_v0.02.hpp"
#include "global_variables.hpp"

using namespace std;

#ifndef coastline_CPP
#define coastline_CPP

/********************************************
* Coastline Object Initialisation Functions *
********************************************/

void Coastline::Initialise()
{
	cout << "\nCoastline.Initialise: Error, initialised an empty Coastline object" << endl;
	exit(EXIT_FAILURE);
}

void Coastline::Initialise(string xyfilename)
{

	/*	Read coastline.XY text file */
	cout << "\nCoastline.Initialise: Initialising Coastline from XY file: " << xyfilename << endl;

	int Time=0;
	ReadCoast(xyfilename, Time);

	//Populate empty vectors
	vector<double> EmptyVector(NoNodes, -9999);
	Distance = EmptyVector;
	Orientation = EmptyVector;
	FluxOrientation = EmptyVector;
	RegionalOrientation = EmptyVector;
	Curvature = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	CellWidth = EmptyVector;
	LongshoreFlux = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	Shadows = EmptyVector;
	ShadowZoneWaveDirection = EmptyVector;
	ShadowZoneWaveHeight = EmptyVector;
	e1 = EmptyVector;
	e2 = EmptyVector;
	
	XMidPoints = EmptyVector;
	YMidPoints = EmptyVector;
	ShorefaceDepth = EmptyVector;
	X0 = EmptyVector;
	Y0 = EmptyVector;
	
	Volume = EmptyVector;
	vector<int> Zeros(NoNodes, 0);
	Fixed = Zeros;
	FlagTriangles = Zeros;
	EvolveFlag = Zeros;
		
	ClosureDepth = 10.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	ShorefaceSlope = 0.02;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	BeachWidth = 20.0;			//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	
	//Calculate Orientation, Slope and Curvature
	CalculateMorphology();
	CalculateMeanNodeSpacing();
	DesiredNodeSpacing = MeanNodeSpacing;

}

void Coastline::Initialise(string xyfilename, int StartTime)
{

	/*	Read coastline.XY text file */
	cout << "\nCoastline.Initialise: Initialising Coastline from XY file: " << xyfilename << " at Time: " << StartTime << endl;
	ReadCoast(xyfilename, StartTime);
	
	//Populate empty vectors
	vector<double> EmptyVector(NoNodes, -9999);
	Distance = EmptyVector;
	Orientation = EmptyVector;
	FluxOrientation = EmptyVector;
	RegionalOrientation = EmptyVector;
	Curvature = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	CellWidth = EmptyVector;
	LongshoreFlux = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	Shadows = EmptyVector;
	ShadowZoneWaveDirection = EmptyVector;
	ShadowZoneWaveHeight = EmptyVector;
	e1 = EmptyVector;
	e2 = EmptyVector;
	Volume = EmptyVector;
	
	XMidPoints = EmptyVector;
	YMidPoints = EmptyVector;
	ShorefaceDepth = EmptyVector;
	X0 = EmptyVector;
	Y0 = EmptyVector;
	
	vector<int> Zeros(NoNodes, 0);
	Fixed = Zeros;
	FlagTriangles = Zeros;
	EvolveFlag = Zeros;
	
	ClosureDepth = 10.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	ShorefaceSlope = 0.02;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	BeachWidth = 20.0;			//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
		
	//Calculate Orientation, Slope and Curvature
	CalculateMorphology();
	CalculateMeanNodeSpacing();
	DesiredNodeSpacing = MeanNodeSpacing;

}

void Coastline::Initialise(int NodeSpacing, double CoastLength, double Trend, int StartBoundaryInput, int EndBoundaryInput)
{

	/*	Initialises the coast as a straight segment with low amplitude noise
		Coast has a fixed length and an orientation/trend. Coordinates of the
		first node in the arrays is [0][0].
		Boundary conditions must be specified */

	//setup boundary conditons
	StartBoundary = StartBoundaryInput;
	EndBoundary = EndBoundaryInput;
	
	cout 	<< "\nCoastline: Initialising Coastline as straight segment"
			<< "\n\t MeanNodeSpacing = " << NodeSpacing
			<< "\n\t CoastLength = " << CoastLength
			//<< "\n\t Trend = " << Trend
			<< "\n\t Boundary Conditions = " << StartBoundary << ", " << EndBoundary << endl << endl;

	//declare temporary variables for setting up coastline object
	double tempX, tempY;
	NoNodes = 0;
	double Dist = 0;
	MeanNodeSpacing = (double)NodeSpacing;
	DesiredNodeSpacing = MeanNodeSpacing;
	
	//initialise coast
	while (Dist<CoastLength)
	{
		if (Dist==0)
		{
			tempX = 0; tempY = 0;
			X.push_back(tempX);
			Y.push_back(tempY);
			NoNodes += 1;
			
			tempX = tempX + 4*MeanNodeSpacing*sin(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
			tempY = tempY + 4*MeanNodeSpacing*cos(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
			X.push_back(tempX);
			Y.push_back(tempY);
			NoNodes += 1;
			Dist += 1.0;
		}
		else
		{
			//sort out random addition
			tempX = tempX + MeanNodeSpacing*sin(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
			tempY = tempY + MeanNodeSpacing*cos(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
			X.push_back(tempX);
			Y.push_back(tempY);
			NoNodes += 1;
			Dist += MeanNodeSpacing;
		}
		
	}
	tempX = tempX + 4*MeanNodeSpacing*sin(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
	tempY = tempY + 4*MeanNodeSpacing*cos(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
	X.push_back(tempX);
	Y.push_back(tempY);
	NoNodes += 1;
	
	//Populate empty vectors
	vector<double> EmptyVector(NoNodes, -9999);
	Distance = EmptyVector;
	Orientation = EmptyVector;
	FluxOrientation = EmptyVector;
	RegionalOrientation = EmptyVector;
	Curvature = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	CellWidth = EmptyVector;
	LongshoreFlux = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	Shadows = EmptyVector;
	ShadowZoneWaveDirection = EmptyVector;
	ShadowZoneWaveHeight = EmptyVector;
	e1 = EmptyVector;
	e2 = EmptyVector;
	Volume = EmptyVector;
	
	XMidPoints = EmptyVector;
	YMidPoints = EmptyVector;
	ShorefaceDepth = EmptyVector;
	X0 = EmptyVector;
	Y0 = EmptyVector;
	
	vector<int> Zeros(NoNodes, 0);
	Fixed = Zeros;
	FlagTriangles = Zeros;
	EvolveFlag = Zeros;
	
	ClosureDepth = 10.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	ShorefaceSlope = 0.02;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	BeachWidth = 20.0;			//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	
	//Calculate Orientation, Slope and Curvature
	CalculateMorphology();
}

void Coastline::ReadCoast(string InputFileName, double Time)
{
	/*	Reads a Coastline object from a text file for a specific time (default is 0)
		File format is
								Start Boundary | End Boundary
								Time | X[0] | X[1] | X[2] =====> X[NoNodes]
								Time | Y[0] | Y[1] | Y[2] =====> Y[NoNodes]	*/

	//Params
	int FindTime = 1;
	int TempTime;
	double TempX, TempY;
	string LineString, line;
	NoNodes = 0;

	//Open the file in a filestream, check the file exists and read data into X and Y
	ifstream ReadCoastFile;
	ReadCoastFile.open(InputFileName.c_str());
	if (ReadCoastFile.is_open())
	{
		ReadCoastFile >> StartBoundary >> EndBoundary;

		//Write Boundary Conditions to screen
		if 			(StartBoundary == 1 && EndBoundary == 1) cout << "Coastline.ReadCoast: Boundaries are periodic" << endl;
		else if 	(StartBoundary == 2 && EndBoundary == 2) cout << "Coastline.ReadCoast: Boundaries are fixed" << endl;
		else if 	(StartBoundary == 3 && EndBoundary == 3) cout << "Coastline.ReadCoast: Boundaries are prescribed" << endl;
		else
		{
			cout 	<< "Coastline.ReadCoast: Error in boundary conditions:"
					<< "\n\tBoundaries must be specified as either"
					<< "\n\t\t1 = periodic, \n\t\t2 = fixed, \n\t\t3 = prescribed"
					<< "\n\tCurrently both start and end boundary must be the same type." << endl << endl;
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		cout << "Coastline.ReadCoast: Error, the file " << InputFileName << " has not been read correctly." << endl;
		exit(EXIT_FAILURE);
	}

	//Loop through lines looking for the desired time
	while (FindTime == 1)
	{
		if (ReadCoastFile.is_open())
		{
			ReadCoastFile >> TempTime;

			//search for first instance of correct time
			if (Time == TempTime)
			{
				getline(ReadCoastFile, line);
				stringstream iss;
				iss << line;

				while (iss >> TempX)
				{
					X.push_back(TempX);
					NoNodes += 1;
				}

				iss << "";
				iss.clear();

				ReadCoastFile >> TempTime;
				getline(ReadCoastFile, line);

				iss << line;

				while (iss >> TempY) Y.push_back(TempY);

				FindTime = 0;
			}

			//otherwise discard lines
			else
			{
				getline(ReadCoastFile, line);
				getline(ReadCoastFile, line);
			}
		}
		else
		{
			cout << "Coastline.ReadCoast: Error, the file " << InputFileName << " has not been read correctly." << endl;
			exit(EXIT_FAILURE);
		}
	}
}

void Coastline::WriteCoast(string OutputFileName, double Time)
{
	/*	Writes a coastline object X and Y coordinates to file for a given time
		If the file already exists the data will be appended else a new file is
		created.

		File format is 	Time | X[0] | X[1] | X[2] =====> X[NoNodes]
								Time | Y[0] | Y[1] | Y[2] =====> Y[NoNodes]			*/

	//Print to screen
	cout.flush();
	cout << "Coastline: Time is " << setprecision(1) << fixed << Time << " years\r";

	//test if output file already exists
	int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WriteCoastFile;
	if (FileExists == 0)
	{
		WriteCoastFile.open(OutputFileName.c_str());
		if (WriteCoastFile.is_open()) WriteCoastFile << StartBoundary << " " << EndBoundary << endl;
	}
	WriteCoastFile.close();

	//open output filestream again to append coastline data
	WriteCoastFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteCoastFile.is_open())
	{
		//write X
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i<NoNodes; ++i) WriteCoastFile << setprecision(10) << " " << X[i];
		WriteCoastFile << endl;

		//write Y
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i< NoNodes; ++i) WriteCoastFile << setprecision(10) << " " << Y[i];
		WriteCoastFile << endl;
	}
	
	else
	{
		//report errors
		cout << "Coastline.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
}

void Coastline::CalculateMorphology()
{
	/************************************************************
	Calculates FluxOrientation [i:i+1], Orientation [i-1:i+1],
	Slope and Curvature members of Coastline object from X and Y positions

	Orientation is the direction of vector i-1 to i+1
	FluxOrientation is the direction of vector i to i+1
	Curvature is change in FluxOrientation over total distance between i-1 and i+1

	Currently handles periodic and fixed boundaries (21/11/13)

	!! ** WILL NEED TO ADD LOGIC FOR DIFFERENT BOUNDARIES** !!

	************************************************************/

	//declare some temporary variables
	double dX, dY;
	double SpatialGradient, SpatialGradientLast, Wb, D_sf;
		
	//Reset all vectors (could builda function for this)
	vector<double> EmptyVector(NoNodes,-9999);
	Distance = EmptyVector;
	Orientation = EmptyVector;
	FluxOrientation = EmptyVector;
	RegionalOrientation = EmptyVector;
	Curvature = EmptyVector;
	CellWidth = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	LongshoreFlux = EmptyVector;
	Shadows = EmptyVector;
	ShadowZoneWaveDirection = EmptyVector;
	ShadowZoneWaveHeight = EmptyVector;
	e1 = EmptyVector;
	e2 = EmptyVector;
	Volume = EmptyVector;
	
	XMidPoints = EmptyVector;
	YMidPoints = EmptyVector;
	ShorefaceDepth = EmptyVector;
	X0 = EmptyVector;
	Y0 = EmptyVector;
	
	vector<int> Zeros(NoNodes,0);
	FlagTriangles = Zeros;
	EvolveFlag = Zeros;
	TriangleFlag = 0;
	TotalVolume = 0;
	
	for (int i=0; i<NoNodes; ++i)
	{
		//Handle StartBoundary
		if (i==0)
		{
			//First get flux orientation (2-node)
			dX = X[1]-X[0];
			dY = Y[1]-Y[0];
			Distance[i] = 0;
			Distance[i+1] = sqrt(dX*dX + dY*dY);
			
			SpatialGradient = dY/dX;
			//convert to azimuths
			if (dX == 0 && dY < 0) FluxOrientation[i] = 180.;
			else if (dX == 0 && dY > 0) FluxOrientation[i] = 0.;
			else if (dX > 0) FluxOrientation[i] = (180./M_PI)*(M_PI*0.5 - atan(SpatialGradient));
			else if (dX < 0) FluxOrientation[i] = (180./M_PI)*(M_PI*1.5 - atan(SpatialGradient));
			else cout << "This is impossible!" << endl;

			//Get MidPoints
			XMidPoints[i] = (X[i+1]+X[i])/2.;
			YMidPoints[i] = (Y[i+1]+Y[i])/2.;
			
			//Next get curvature and orientation (3-node)
			//PERIODIC BOUNDARY
			if (StartBoundary == 1)
			{
				SpatialGradientLast = Y[NoNodes-1]-Y[NoNodes-2]/X[NoNodes-1]-X[NoNodes-2];
				dX = X[NoNodes-1]-X[NoNodes-2]+X[1]-X[0];
				dY = Y[NoNodes-1]-Y[NoNodes-2]+Y[1]-Y[0];
				Curvature[i] = (SpatialGradient-SpatialGradientLast)/sqrt(dX*dX + dY*dY);
			}
			//FIXED BOUNDARY
			else if (StartBoundary == 2)
			{
				dX = 2*dX;
				Curvature[i] = 0;
			}
			//BOUNDARY ERROR
			else
			{
				cout << "Coastline.CalculateMorphology: Error, did not recognise StartBoundary" << endl << endl;
				exit(EXIT_FAILURE);
			}
			dY = 2*dY;
			
			//convert to azimuths
			if (dX == 0 && dY < 0) Orientation[i] = 180.;
			else if (dX == 0 && dY > 0) Orientation[i] = 0.;
			else if (dX > 0) Orientation[i] = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
			else if (dX < 0) Orientation[i] = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
			else cout << "c'est impossible!" << endl;

			//get regional orientation
			RegionalOrientation[i] = Orientation[i];
			
			//calculate cell width
			CellWidth[i] = (sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i]));
		}

		//Deal with EndBoundary
		else if (i==NoNodes-1)
		{
			//Get MidPoints
				XMidPoints[i] = X[i] + (X[i]-X[i-1]);
				YMidPoints[i] = Y[i] + (Y[i]-Y[i-1]);
				
			//PERIODIC BOUNDARY
			if (EndBoundary == 1)
			{
				//First get flux orientation (2-node)
				FluxOrientation[i] = FluxOrientation[0];
				Orientation[i] = Orientation[0];
				RegionalOrientation[i] = RegionalOrientation[0];
				Curvature[i] = Curvature[0];
				CellWidth[i] = CellWidth[0];
			}
			//FIXED BOUNDARY
			else if (EndBoundary == 2)
			{
				FluxOrientation[i] = FluxOrientation[i-1];
				Orientation[i] = Orientation[i-1];
				RegionalOrientation[i] = RegionalOrientation[i-1];
				Curvature[i] = 0;
				CellWidth[i] = (sqrt(pow(X[i]-X[i-1],2.0) + pow(Y[i]-Y[i-1],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i]));
			}
			//BOUNDARY ERROR
			else
			{
				cout << "Coastline.CalculateMorphology: Error, did not recognise EndBoundary" << endl << endl;
				exit(EXIT_FAILURE);
			}
		}
		//Finally do all the bits inbetween
		else
		{
			//First get flux orientation (2-node)
			dX = X[i+1]-X[i];
			dY = Y[i+1]-Y[i];
			Distance[i+1] = Distance[i] + sqrt(dX*dX + dY*dY);
			
			//don't allow divide by zero, just make spatial gradient very large
			if (dX != 0) SpatialGradient = dY/dX;
			
			//convert to azimuths
			if (dX == 0 && dY < 0) FluxOrientation[i] = 180.;
			else if (dX == 0 && dY > 0) FluxOrientation[i] = 0.;
			else if (dX > 0) FluxOrientation[i] = (180./M_PI)*(M_PI*0.5 - atan(SpatialGradient));
			else if (dX < 0) FluxOrientation[i] = (180./M_PI)*(M_PI*1.5 - atan(SpatialGradient));
			else cout << "es imposible!" << endl;

			//Next get curvature and orientation (3-node)
			dX = (X[i+1]-X[i-1]);
			dY = (Y[i+1]-Y[i-1]);
			Curvature[i] = (tan((M_PI/180.)*(FluxOrientation[i]-FluxOrientation[i-1])))/sqrt(dX*dX + dY*dY);
			//convert to azimuths
			if (dX == 0 && dY < 0) Orientation[i] = 180.;
			else if (dX == 0 && dY > 0) Orientation[i] = 0.;
			else if (dX > 0) Orientation[i] = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
			else if (dX < 0) Orientation[i] = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
			else cout << "es ist unmÃ¶glich" << endl;

			//Get Regional Orientation
			dX = (X[i+2]-X[i-1]);
			dY = (Y[i+2]-Y[i-1]);
			//convert to azimuths
			if (dX == 0 && dY < 0) RegionalOrientation[i] = 180.;
			else if (dX == 0 && dY > 0) RegionalOrientation[i] = 0.;
			else if (dX > 0) RegionalOrientation[i] = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
			else if (dX < 0) RegionalOrientation[i] = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
			else cout << "ei bod yn amhosibl" << endl;
			
			if (i==NoNodes-2) RegionalOrientation[i] = RegionalOrientation[i-1];
			
			//Get MidPoints
			XMidPoints[i] = (X[i+1]+X[i])/2.;
			YMidPoints[i] = (Y[i+1]+Y[i])/2.;
			
			//Cell Width parallel to orientation AT the node position
			CellWidth[i] = (0.5*sqrt(pow(X[i]-X[i-1],2.0) + pow(Y[i]-Y[i-1],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i-1]));
			CellWidth[i] += (0.5*sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i]));
						
			//Get cell width at bottom of the shoreface
			e1[i] = Orientation[i]-FluxOrientation[i-1];
			e2[i] = FluxOrientation[i] - Orientation[i];
			Wb = CellWidth[i] + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));	
			
			//Flag triangles and calculate mass of sand assuming a fixed beach width
			if (Wb < 0) 
			{
				//Flag triangles
				FlagTriangles[i] = 1;
				TriangleFlag = 1;
				
				//Get effective shoreface depth
				D_sf = -(CellWidth[i]*ShorefaceSlope)/(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
				
				//Calculate volume
				Volume[i] = D_sf*CellWidth[i]*BeachWidth + 
							D_sf*D_sf*BeachWidth*(tan((M_PI/180.)*e1[i])+ tan((M_PI/180.)*e2[i]))/ShorefaceSlope + 
							D_sf*0.5*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]))*BeachWidth*BeachWidth;
				TotalVolume += Volume[i];
			}
			else
			{
				Volume[i] = ClosureDepth*CellWidth[i]*BeachWidth + 
							ClosureDepth*ClosureDepth*BeachWidth*(tan((M_PI/180.)*e1[i])+ tan((M_PI/180.)*e2[i]))/ShorefaceSlope + 
							ClosureDepth*0.5*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]))*BeachWidth*BeachWidth;
				TotalVolume += Volume[i];
			}
		}
	}
}

void Coastline::CheckNodeSpacing()
{
	//declare variables
	int SpacingFlag = 0;
	double Distance,MeanOrientation,X0,Y0;//,M1,M2; //, VolumeOld, VolumeNew, VolChange, a, b, c;
	//need to recalculate morphology and start check again after each new node is added
	while (SpacingFlag == 0)
	{
		for (int i=0; i<NoNodes-1; ++i)
		{
			if (i<1) {}
			else if (i>NoNodes-3) {}
//			else if (FlagTriangles[i] == 1)
//			{
//				//delete triangular node
//				X.erase(X.begin()+i); 
//				Y.erase(Y.begin()+i); 
//				Fixed.erase(Fixed.begin()+i+1);
//				NoNodes -= 1;
//				SpacingFlag = 0;
//				
//				//check volumes and update coastal morphology
//				VolumeOld = TotalVolume;
//				CalculateMorphology();
//				VolumeNew = TotalVolume;
//				VolChange = VolumeNew-VolumeOld;
//				
//				//remove required volume from adjacent nodes
//				if (i == 2)
//				{
//					//only change volume of next node down coast (cell i)
//					e1[i] = Orientation[i]-FluxOrientation[i-1];
//					e2[i] = FluxOrientation[i] - Orientation[i];
//					
//					//quadratic coefficients
//					a = tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]);
//					b = 2.*CellWidth[i] + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
//					c = -2.*VolChange/ClosureDepth;
//					
//					//Move the node
//					if (a==0) PositionChange[i] = VolumeChange[i]/(CellWidth[i]*ClosureDepth);
//					else PositionChange[i] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//					
//					//Decompose Position change into X and Y components
//					//Going to need a bit more logic here to handle different orientations?
//					//current setup x change is negative, y change is positive
//					X[i] -= PositionChange[i]*cos((M_PI/180.)*Orientation[i]);
//					Y[i] += PositionChange[i]*sin((M_PI/180.)*Orientation[i]);
//				}
//				else if (i == NoNodes-2)
//				{
//					//only change volume of next cell up the coast
//					e1[i-1] = Orientation[i-1]-FluxOrientation[i-2];
//					e2[i-1] = FluxOrientation[i-1] - Orientation[i-1];
//					
//					//quadratic coefficients
//					a = tan((M_PI/180.)*e1[i-1]) + tan((M_PI/180.)*e2[i-1]);
//					b = 2.*CellWidth[i-1] + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*e1[i-1]) + tan((M_PI/180.)*e2[i-1]));
//					c = -2.*VolChange/ClosureDepth;
//					
//					//Move the node
//					if (a==0) PositionChange[i-1] = VolumeChange[i-1]/(CellWidth[i-1]*ClosureDepth);
//					else PositionChange[i-1] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//					
//					//Decompose Position change into X and Y components
//					//Going to need a bit more logic here to handle different orientations?
//					//current setup x change is negative, y change is positive
//					X[i-1] -= PositionChange[i-1]*cos((M_PI/180.)*Orientation[i-1]);
//					Y[i-1] += PositionChange[i-1]*sin((M_PI/180.)*Orientation[i-1]);
//				}
//				else
//				{
//					//pass half the volume to each cell
//					VolChange = VolChange/2.;
//					break 
//					for (int j=i-1; j<i+1; ++j)
//					{ 
//						e1[j] = Orientation[j]-FluxOrientation[j-1];
//						e2[j] = FluxOrientation[j] - Orientation[j];
//						
//						//quadratic coefficients
//						a = tan((M_PI/180.)*e1[j]) + tan((M_PI/180.)*e2[j]);
//						b = 2.*CellWidth[j] + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*e1[j]) + tan((M_PI/180.)*e2[j]));
//						c = -2.*VolChange/ClosureDepth;
//					
//						//Move the node
//						if (a==0) PositionChange[j] = VolumeChange[j]/(CellWidth[j]*ClosureDepth);
//						else PositionChange[j] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//						
//						//Decompose Position change into X and Y components
//						//Going to need a bit more logic here to handle different orientations?
//						//current setup x change is negative, y change is positive
//						X[j] -= PositionChange[j]*cos((M_PI/180.)*Orientation[j]);
//						Y[j] += PositionChange[j]*sin((M_PI/180.)*Orientation[j]);
//					}
//				}
//				
//				break;
//			}
			else
			{
				Distance = sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));
				
				//if distance between two cells is too small, 
				//replace with a linearly interpolated point between
				if (Distance < 0.6*DesiredNodeSpacing)
				{
					//if cell is too narrow remove its adjacent nodes
					if (i == 1 || i == NoNodes-3) 
					{
						X.erase(X.begin()+i+1); 
						Y.erase(Y.begin()+i+1); 
						Fixed.erase(Fixed.begin()+i+1);
						NoNodes -= 1;
						SpacingFlag = 0;
					}
					else
					{
						X[i] = ((X[i+1]+X[i])/2.);
						Y[i] = ((Y[i+1]+Y[i])/2.);
						X.erase(X.begin()+i+1); 
						Y.erase(Y.begin()+i+1); 
						Fixed.erase(Fixed.begin()+i+1);
						NoNodes -= 1;
						SpacingFlag = 0;
					}
					break;
				}
				else if 
				(Distance > 1.5*DesiredNodeSpacing) 
				{
//					//if cell too far apart add new node
//					//just do with a straightline for now
					if (i == 1 || i == NoNodes-3 || e1[i+1] > 0)
					{
						X.insert(X.begin()+i+1,(X[i]+X[i+1])/2.);
						Y.insert(Y.begin()+i+1,(Y[i]+Y[i+1])/2.);
						Fixed.insert(Fixed.begin()+i+1,0.0);
						NoNodes += 1;
						SpacingFlag = 0;
						break;
					}
					else
					{
						MeanOrientation = (Orientation[i] + Orientation[i+1])/2.;
						X0 = ((X[i]+X[i+1])/2.) + 1.*cos((M_PI/180.)*MeanOrientation);
						Y0 = ((Y[i]+Y[i+1])/2.) - 1.*sin((M_PI/180.)*MeanOrientation);
						X.insert(X.begin()+i+1,X0);
						Y.insert(Y.begin()+i+1,Y0);
						Fixed.insert(Fixed.begin()+i+1,0.0);
						NoNodes += 1;
						SpacingFlag = 0;
						break;
					
						
//						//if cell too far apart add new node
//						//interpolate based on orientation of cell i+1 and i?
//						//Find intersection of lines along orientation direction
//						//get distance the length of the shore across the shadow zone
//						//get gradients of intersecting lines
//						M1 = 1./tan((M_PI/180.)*Orientation[i]);
//						M2 = 1./tan((M_PI/180.)*Orientation[i+1]);
//						
//						//find point of intersection
//						X0 = (Y[i+1] - M2*X[i+1] + M1*X[i] - Y[i])/(M1 - M2);
//						Y0 = M1*X0 + Y[i] - M1*X[i];
//						
//						X.insert(X.begin()+i+1,X0);
//						Y.insert(Y.begin()+i+1,Y0);
//						Fixed.insert(Fixed.begin()+i+1,0.0);
//						NoNodes += 1;
//						SpacingFlag = 0;
//						break;
					}
				}
				else SpacingFlag = 1;
			}
		}
		CalculateMorphology();
	}
}

void Coastline::CalculateMeanNodeSpacing()
{
	//Returns an integer of the mean node spacing from X and Y data
	double Distance;
	double dx, dy;
	double TotalDistance = 0;

	for (int i=0; i<NoNodes; ++i)
	{
		//For the last cell Distance is the same as for cell i=0
		if (i==NoNodes-1)
		{
			//Calculate distance between 0 and 1
			dx = fabs(X[1]-X[0]);
			dy = fabs(Y[1]-Y[0]);
			Distance = sqrt(dx*dx + dy*dy);
		}
		else
		{
			//Calculate distance between i and i+1
			dx = fabs(X[i+1]-X[i]);
			dy = fabs(Y[i+1]-Y[i]);
			Distance = sqrt(dx*dx + dy*dy);
		}
		TotalDistance += Distance;
	}

	//Calculate MeanNodeSpacing
	MeanNodeSpacing = (int)(TotalDistance/NoNodes);
}

void Coastline::TransformWaves()
{
	//run transform wave on every node along the coast
	for (int i=0; i<NoNodes; ++i) TransformWave(i);
}

void Coastline::TransformWave(int i)
{
	//Transforms waves assume shore-parallel contours and fixed shoreface gradient
	//Wave base depth calculated based on wave conditions?

	// C_0 is offshore wave celerity, L_0 is offshore wave length,
	// C is wave speed, H is wave height, L is wavelength, alpha is wave angle to the coast
	// h is depth, n is a shoaling factor, k is wave number,
	//	Ks is shoaling coefficient, Kr is refraction coefficient
	// See Komar (1998)

	double C_0, L_0, Alpha_0, Alpha_0_Last, Alpha_0_Next, H_0, Theta_0;
	double C, H, L, Alpha;
	double h, n, k, Ks, Kr;

	int BreakCondition = 0;

	if (Shadows[i] == 1)
	{
		H_0 = ShadowZoneWaveHeight[i];
		Theta_0 = ShadowZoneWaveDirection[i];
	}
	else 
	{
		H_0 = OffshoreWaveHeight;
		Theta_0 = OffshoreWaveDirection;
	}
	
	C_0 = (g*OffshoreWavePeriod)/(2.0*M_PI);			//Deep water wave speed (m/s)
	L_0 = C_0*OffshoreWavePeriod;						//Deep water wave length (m)

	BreakCondition = 0;	//for testing wave breaking
	h = 3.*H_0;		// water depth at wave base (metres) calculate this later based on period and height

		//determine alpha_0 angle between coast and wave crest approach
		if (Shadows[i] == 1)
	   	{
	   		if (ShadowZoneWaveDirection[i] <= FluxOrientation[i]) Alpha_0 = FluxOrientation[i]-ShadowZoneWaveDirection[i]-90.;
    		else if (ShadowZoneWaveDirection[i] > FluxOrientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-ShadowZoneWaveDirection[i];
	    	else Alpha_0= 270.-(ShadowZoneWaveDirection[i]-FluxOrientation[i]);
	   	}
	   	else if ((Shadows[i] == 2) && (i > 2))
	   	{
	   		if (Theta_0 <= FluxOrientation[i-1]) Alpha_0 = FluxOrientation[i-1]-Theta_0-90.;
    		else if (Theta_0 > FluxOrientation[i]+270.) Alpha_0 = (FluxOrientation[i-1]+270.)-Theta_0;
	   		else Alpha_0= 270.-(Theta_0-FluxOrientation[i-1]);
	   	}
	   	else if ((Shadows[i] == 4) && (i < NoNodes-3))
	   	{
	   		//determine alpha_0 angle between coast and wave crest approach
			if (Theta_0 <= FluxOrientation[i+1]) Alpha_0 = FluxOrientation[i+1]-Theta_0-90.;
    		else if (Theta_0 > FluxOrientation[i+1]+270.) Alpha_0 = (FluxOrientation[i+1]+270.)-Theta_0;
	   		else Alpha_0= 270.-(Theta_0-FluxOrientation[i+1]);
	   	}
	   	else
	   	{
	   		if (Theta_0 <= FluxOrientation[i+1]) Alpha_0 = FluxOrientation[i]-Theta_0-90.;
    		else if (Theta_0 > FluxOrientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-Theta_0;
	    	else Alpha_0= 270.-(Theta_0-FluxOrientation[i]);
	    }
	    
	    if (i > 0)
		{
	   		//need some additional logic here for inside shadow zone?
	   		//Alpha_0_Next based on ShadowZoneWaveDirection, not Theta_0
	   		if (Shadows[i-1] == 1)
	   		{
	   			if (ShadowZoneWaveDirection[i-1] <= Orientation[i-1]) Alpha_0_Last = FluxOrientation[i-1]-ShadowZoneWaveDirection[i-1]-90.;
	    		else if (ShadowZoneWaveDirection[i-1] > Orientation[i-1]+270.) Alpha_0_Last = (FluxOrientation[i-1]+270.)-ShadowZoneWaveDirection[i-1];
		   		else Alpha_0_Last = 270.-(ShadowZoneWaveDirection[i]-FluxOrientation[i-1]);
		   	}
		   	else
		   	{
			   	if (Theta_0 <= Orientation[i-1]) Alpha_0_Last = FluxOrientation[i-1]-Theta_0-90.;
	    		else if (Theta_0 > Orientation[i-1]+270.) Alpha_0_Last = (FluxOrientation[i-1]+270.)-Theta_0;
	   			else Alpha_0_Last = 270.-(Theta_0-FluxOrientation[i-1]);
	   		}
	   	}
	   	else Alpha_0_Last = Alpha_0;
	   	
	   	if (i<NoNodes-1)
	   	{
	   		//need some additional logic here for inside shadow zone?
	   		//Alpha_0_Next based on ShadowZoneWaveDirection, not Theta_0
	   		if (Shadows[i+1] == 1)
	   		{
	   			if (ShadowZoneWaveDirection[i+1] <= Orientation[i+1]) Alpha_0_Next = FluxOrientation[i+1]-ShadowZoneWaveDirection[i+1]-90.;
	    		else if (ShadowZoneWaveDirection[i+1] > Orientation[i+1]+270.) Alpha_0_Next = (FluxOrientation[i+1]+270.)-ShadowZoneWaveDirection[i+1];
		   		else Alpha_0_Next = 270.-(ShadowZoneWaveDirection[i+1]-FluxOrientation[i+1]);
		   	}
		   	else
		   	{
		   		if (Theta_0 <= Orientation[i+1]) Alpha_0_Next = FluxOrientation[i+1]-Theta_0-90.;
	    		else if (Theta_0 > Orientation[i+1]+270.) Alpha_0_Next = (FluxOrientation[i+1]+270.)-Theta_0;
		   		else Alpha_0_Next = 270.-(Theta_0-FluxOrientation[i+1]);
		   	}
	   	}
	   	
	   	//if high angle waves then use flux orientation of the previous node (updrift)
	   	//first if wave is offshore no transport
	   	//this bit is messy could do with a reorganise
	   	if (Alpha_0 <= -90. || Alpha_0 >= 90.) Alpha_0 = 0;
	   	//next if alpha_0 and alpha_0_last have the same sign
	   	else if ((Alpha_0_Last > 0) && (Alpha_0 > 0))
	   	{
	    	if ((Alpha_0_Last < 45.) && (Alpha_0 > 45.)) Alpha_0 = 45.;
	    	//else if ((Alpha_0_Last > 45.) && (Alpha_0 < 45.)) Alpha_0 = 45.;
	    	else if (Alpha_0 > 45.) Alpha_0 = Alpha_0_Last;
	    }
	    //next if alpha_0 and alpha_0_next have the same sign
	    else if ((Alpha_0_Next < 0) && (Alpha_0 < 0))
	    {
	    	if ((Alpha_0_Next > -45.) && (Alpha_0 < -45.)) Alpha_0 = -45.;
	    	//else if ((Alpha_0_Next < -45.) && (Alpha_0 > -45.)) Alpha_0 = -45.;
	    	else if (Alpha_0 < -45.) Alpha_0 = Alpha_0_Next;
	    }
	    //else for high angle waves use updrift orientation
	    else if (Alpha_0 > 45.0 && Alpha_0_Last > 0)
	    {
	    	if (Theta_0 <= Orientation[i]) Alpha_0 = FluxOrientation[i-1]-Theta_0-90.;
	    	else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (FluxOrientation[i-1]+270.)-Theta_0;
		   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i-1]);
	    }
		else if (Alpha_0 < -45.0 && Alpha_0_Next < 0)
		{
			if (Theta_0 <= Orientation[i]) Alpha_0 = FluxOrientation[i+1]-Theta_0-90.;
	    	else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (FluxOrientation[i+1]+270.)-Theta_0;
		   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i+1]);
		}
	//}
	
	if (Alpha_0 == 0) 
	{
		BreakingWaveHeight[i] = 0;
		BreakingWaveAngle[i] = 0;
	}
	else
	{
		while (BreakCondition == 0)
		{
			//Calculate new wave height
			L = L_0*sqrt(tanh((2.*M_PI*h)/L_0));		//Wavelength (m) in intermediate-shallow waters
			//L = L_0*pow(tanh(pow(pow(2.*M_PI/T,2.)*h/g,.75)),2.0/3.0); //Fenton & McKee (1990) formulation
		   	C = C_0*tanh((2*M_PI*h)/L);					//Wave speed (m/s) set by water depth, L and C_0
			k = 2*M_PI/L;								//Wave number (1/m)
	   	    n = ((2*h*k)/(sinh(2*h*k)) + 1)/2;			//Shoaling factor
	   	    Ks = sqrt(C_0/(n*C*2));							//Shoaling coefficient
	   	    Alpha = (180./M_PI)*asin((C/C_0)*sin((M_PI/180.)*Alpha_0));		//update theta
	   	    Kr = sqrt(cos((M_PI/180.)*Alpha_0)/cos((M_PI/180.)*Alpha));		//refraction coefficient
	   	    H = H_0*Ks*Kr;										//calculate new wave height

			//test if wave breaks
			if (H > h*0.78)
			{
				BreakCondition = 1;
				BreakingWaveHeight[i] = H;
				BreakingWaveAngle[i] = Alpha;
			}

			//update water depth
			h -= 0.05;

			//Catch negative water depths!
			//Must be Alpha -> 90 so assume zero transport and set wave height and angle to 0
			if (h < 0)
			{
				BreakingWaveHeight[i] = 0;
				BreakingWaveAngle[i] = 0;
				BreakCondition = 1; // dont bother with refraction for this cell
			}
		}
	}
}

void Coastline::GetShadows()
{
	/*description goes here
	
	Loops through coast from start to end and end to start to check 
	whether the coastline casts shadows given the offshore wave direction
	There are five codes to the Shadows vector:
		#0: No shadow
		#1: In Shadow
		#2: Casts a shadow forward along the coast
		#3: Casts a shadow backward along the coast
		#4: Special case for cell just downdrift of backward shadow (#3)

	These Codes are required for handling sediment transport andc
	for the refraction/diffraction code
	
	*/
	
	//declare temporary variables
	double XProd, X1, Y1, X2, Y2, dX12, dY12, X3, Y3, X4, Y4, dX34, dY34, dX31, dY31, S, T;
	double ShadowAngle;
	int i, ShadowStart, ShadowEnd, CasterFlag, XProdPos;
	
	//reset shadows
	vector<double> Zeros(NoNodes, 0);
	vector<double> Empty(NoNodes, -9999);
	Shadows = Zeros;
	ShadowZoneWaveDirection = Empty;
	ShadowZoneWaveHeight = Empty;
	ShadowFlag = 0;
	
	//get incloming wave direction and project to shadow angle
	ShadowAngle = OffshoreWaveDirection+180.;
	if (ShadowAngle >= 360.) ShadowAngle -= 360.;

	//loop through coast
	for (i=1; i<NoNodes-1; ++i)
	{
	    if 	(	(OffshoreWaveDirection < FluxOrientation[i]-180.) 
	    	||	(OffshoreWaveDirection > FluxOrientation[i]) 
	    	||	(OffshoreWaveDirection < FluxOrientation[i-1]-180.)) 
	    {
	    	Shadows[i] = 1;
	    	ShadowFlag = 1;
	    }

	    else
	    {
			//set point for casting shadow    
        	X1 = X[i];
	        Y1 = Y[i];
        
			//create a hypothetical point for offshore vector in direction of OffshoreWaveDirection
	        X2 = X[i] + 100000.*sin((M_PI/180.)*OffshoreWaveDirection);
    	    Y2 = Y[i] + 100000.*cos((M_PI/180.)*OffshoreWaveDirection);
           	dX12 = X2-X1;
	        dY12 = Y2-Y1;
        
        	for (int j=0; j<NoNodes-1; ++j)
        	{
	            if ((j < i-1) || (j > i))
	            {
	            	//Get points for coastal segment and diffs
                	X3 = X[j];
                	Y3 = Y[j];
                	X4 = X[j+1];
                	Y4 = Y[j+1];
                	dX34 = X4-X3;
                	dY34 = Y4-Y3;
                	
                	//Find the cross product of the two vectors
	                XProd = dX12*dY34 - dX34*dY12;
        
					//Xproduct == 0 for parallel lines        
                	if (XProd != 0)
                	{
                    	if (XProd > 0) XProdPos = 1;
	                    else XProdPos = 0;
        				
                    	//assign third test segment
	                    dX31 = X1-X3;
    	                dY31 = Y1-Y3;
                    
                    	//get cross products
	                    S = dX12*dY31 - dY12*dX31;
	                    T = dX34*dY31 - dY34*dX31;
                    
                    	//logic for collision occurence
                    	if ((S < 0) == XProdPos) {}
                    	else if ((T < 0) == XProdPos) {}
                    	else if (((S > XProd) == XProdPos) or ((T > XProd) == XProdPos)) {}
                        else 
                       	{
                       		Shadows[i] = 1; 
                       		ShadowFlag = 1; 
                       		break;
                       	}
					}
				}
			}
		}
	}
	
	//if no shadows break out of function
	if (ShadowFlag == 0) return;
	
	//Otherwise find which nodes cast the shadows 
	i=0;
		  
	while (i < NoNodes-1)
	{
    	++i;

		if (Shadows[i] == 1)
		{
        	ShadowStart = i;
        	ShadowEnd = i;
        	while (Shadows[ShadowEnd] == 1) ++ShadowEnd;
        	if (ShadowEnd >= NoNodes) 
        	{
        		cout << "ShadowEnd >= NoNodes" << endl;
        		ShadowEnd = NoNodes-2;
        	}
			
			//set point for casting shadow    
			X1 = X[ShadowStart];
			Y1 = Y[ShadowStart];
		    
	    	//create a hypothetical point for vector in direction of ShadowAngle
		    X2 = X[ShadowStart] + 100000.*sin((M_PI/180.)*ShadowAngle);
		    Y2 = Y[ShadowStart] + 100000.*cos((M_PI/180.)*ShadowAngle);
	    	dX12 = X2-X1;
		    dY12 = Y2-Y1;
	    
	    	CasterFlag = 0;
	    
	    	for (int j=ShadowStart+1; j<ShadowEnd+1; ++j)
	    	{
	        	//Get points for coastal segment and diffs
	        	X3 = X[j];
	        	Y3 = Y[j];
	        	X4 = X[j+1];
	        	Y4 = Y[j+1];
	        	dX34 = X4-X3;
	        	dY34 = Y4-Y3;
	        
	        	//Find the cross product of the two vectors
		        XProd = dX12*dY34 - dX34*dY12;
	
	        	//Xproduct == 0 for parallel lines
	        	if (XProd != 0)
	        	{
					if (XProd > 0) XProdPos = 1;
		            else XProdPos = 0;
					
					//assign third test segment
		            dX31 = X1-X3;
		            dY31 = Y1-Y3;
		            
		            //get cross products
		            S = dX12*dY31 - dY12*dX31;
		            T = dX34*dY31 - dY34*dX31;
		            
		            //logic for collision occurence
		            if ((S < 0) == XProdPos) {}
		            else if ((T < 0) == XProdPos) {}
		            else if (((S > XProd) == XProdPos) or ((T > XProd) == XProdPos)) {}
		            else 
		            {
		            	CasterFlag = 1;
		            	break;
		            }
				}
	    	}
            
    		//if shadowstart casts the shadow set Shadows to 2
        	//otherwise shadowend must be the caster
    		if (CasterFlag == 1) Shadows[ShadowStart] = 2;
    		else 
    		{
        	    //loop back and do intersection analysis from ShadowEnd
        	    //set point for casting shadow    
			    X1 = X[ShadowEnd];
			    Y1 = Y[ShadowEnd];
		    
		    	//create a hypothetical point for vector in direction of ShadowAngle
			    X2 = X[ShadowEnd] + 100000.*sin((M_PI/180.)*ShadowAngle);
			    Y2 = Y[ShadowEnd] + 100000.*cos((M_PI/180.)*ShadowAngle);
		    	dX12 = X2-X1;
			    dY12 = Y2-Y1;
		    
		    	CasterFlag = 0;
		    
		    	for (int j=ShadowEnd-2; j>ShadowStart-2; --j)
		    	{
		        	if (j == 0)
		        	{
		        		CasterFlag = 1;
		        		break;
		        	}
		        	
		        	//Get points for coastal segment and diffs
		        	X3 = X[j];
		        	Y3 = Y[j];
		        	X4 = X[j+1];
		        	Y4 = Y[j+1];
		        	dX34 = X4-X3;
		        	dY34 = Y4-Y3;
		        
		        	//Find the cross product of the two vectors
			        XProd = dX12*dY34 - dX34*dY12;
		
		        	//Xproduct == 0 for parallel lines        
		        	if (XProd != 0)
		        	{
						if (XProd > 0) XProdPos = 1;
			            else XProdPos = 0;
						
						//assign third test segment
			            dX31 = X1-X3;
			            dY31 = Y1-Y3;
			            
			            //get cross products
			            S = dX12*dY31 - dY12*dX31;
			            T = dX34*dY31 - dY34*dX31;
			            
			            //logic for collision occurence
			            if ((S < 0) == XProdPos) {}
			            else if ((T < 0) == XProdPos) {}
			            else if (((S > XProd) == XProdPos) or ((T > XProd) == XProdPos)) {}
			            else 
			            {
			            	CasterFlag = 1; 
			            	break;
			            }
					}
		    	}
		    	
		    	if (CasterFlag == 1)
		    	{
        	    	Shadows[ShadowEnd] = 3;
        	    	if (ShadowEnd > 1)
        	    	{
        	    		Shadows[ShadowEnd-1] = 4;
        	    		Shadows[ShadowEnd-2] = 1;
        	    	}
        	    }
        	    else Shadows[ShadowStart] = 2;
        	}    
        
        	//continue from after shadowed section
        	i = ShadowEnd;
        }
	}
	//check here
	cout << "";
}

void Coastline::RefractDiffractShadowZone()
{
	/*
	
	description goes here
	
	Function to compute wave angles and heights in the shadow zone
	based on simple rule set
	
	*/
	
	//declare temporary variables
	double ShadowAngle, SpatialGradient, ShadowZoneDistance, dX, dY, X0, Y0, M1, M2, Dist, AngleDiff, TempWaveHeight;
	int i, j, k, ShadowStart, ShadowEnd;
	
	//get incloming wave direction and project to shadow angle
	ShadowAngle = OffshoreWaveDirection+180.;
	if (ShadowAngle >= 360.) ShadowAngle -= 360.;

	//loop through
	for (i=1; i<NoNodes-1; ++i)
	{
		//catch a shadow casting cell
		if (Shadows[i] == 2)
		{
			//find end of shadow zone
			ShadowStart=i;
			ShadowEnd = ShadowStart+1;
			while ((Shadows[ShadowEnd] == 1) && (ShadowEnd<NoNodes-1)) ++ShadowEnd;
			
			//if shadow cast merges with an up coast shadow (shouldn't happen)
			if (Shadows[ShadowEnd] == 4) 
			{
				cout << endl << "Shadows Merge" << endl;
				Shadows[ShadowEnd] = 1;
				Shadows[ShadowEnd+1] = 1;
				ShadowEnd += 1;
			}
			
			//setup iterator to loop back across shadow
			ShadowEnd -= 1;
			j = ShadowEnd;
									
			//get distance the length of the shore across the shadow zone
			//get gradients of intersecting lines
			M1 = 1./tan((M_PI/180.)*ShadowAngle);
			M2 = 1./tan((M_PI/180.)*FluxOrientation[ShadowEnd]);
			//find point of intersection
			X0 = (Y[ShadowEnd] - M2*X[ShadowEnd] + M1*X[ShadowStart] - Y[ShadowStart])/(M1 - M2);
			Y0 = M1*X0 + Y[ShadowStart] - M1*X[ShadowStart];
			//Find distance difference
			dX = X0-X[ShadowEnd];
			dY = Y0-Y[ShadowEnd];
			Dist = sqrt(dX*dX + dY*dY);
			ShadowZoneDistance = Distance[ShadowEnd]-Distance[ShadowStart]+Dist;
			
			//loop back through and set BreakingWaveHeight and Direction in the shadow
			while (j > ShadowStart)
			{
				//calculate shadow zone angle
				dX = X[j]-X[ShadowStart];
				dY = Y[j]-Y[ShadowStart];
				SpatialGradient = dY/dX;
							
				//convert to azimuths and determine refracted diffracted wave properties
				if (dX == 0 && dY < 0) AngleDiff = ShadowAngle-180.;
				else if (dX == 0 && dY > 0) AngleDiff = ShadowAngle-360;
				else if (dX > 0) AngleDiff = (180./M_PI)*(M_PI*0.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				else if (dX < 0) AngleDiff = (180./M_PI)*(M_PI*1.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*AngleDiff;
				//if (1.5*AngleDiff > 90.) ShadowZoneWaveHeight[j] = 0;
				//else
				ShadowZoneWaveHeight[j] = OffshoreWaveHeight*0.5*(1-sin((M_PI/180.)*(AngleDiff)));
				--j;
			}
			
			//Reduce wave heights outside the shadow zone over the length that they were redued in the shadow zone
			//This should preserve energy.
			
			if (Shadows[ShadowStart] == 1) cout << "What the fuck!?" << endl << endl << endl;
			
			Dist = (Distance[ShadowEnd]-Distance[ShadowEnd-1])-Dist;
			k = ShadowEnd+1;
			while(Dist <= ShadowZoneDistance)
			{
				ShadowZoneWaveHeight[k] = OffshoreWaveHeight*(1.-0.5*(1.-sin((M_PI/180.)*(90.*(Dist/ShadowZoneDistance)))));
				ShadowZoneWaveDirection[k] = OffshoreWaveDirection;
				Shadows[k] = 1;
				k+=1;
				if (Shadows[k] != 0 || k == NoNodes-1 || k == 0) Dist = ShadowZoneDistance+1;
				else Dist += (Distance[k]-Distance[k-1]);
			}
		}
		
		else if (Shadows[i] == 3)
		{
			//find end of shadow zone
			ShadowStart = i;
			ShadowEnd = i-1;
			while ((Shadows[ShadowEnd] == 1 || Shadows[ShadowEnd] == 4) && (ShadowEnd > 0)) --ShadowEnd;
						
			//get distance the length of the shore across the shadow zone
			//get gradients of intersecting lines1
			M1 = 1./tan((M_PI/180.)*ShadowAngle);
			M2 = 1./tan((M_PI/180.)*FluxOrientation[j]);
			//find point of intersection
			X0 = (Y[ShadowEnd] - M2*X[ShadowEnd] + M1*X[ShadowStart] - Y[ShadowStart])/(M1 - M2);
			Y0 = M1*X0 + Y[ShadowStart] - M1*X[ShadowStart];
			//Find distance difference
			dX = X0-X[ShadowEnd];
			dY = Y0-Y[ShadowEnd];
			Dist = sqrt(dX*dX + dY*dY);
			ShadowZoneDistance = Distance[ShadowStart]-Distance[ShadowEnd]-Dist;

			//setup iterator to loop back across shadow
			j = ShadowEnd;

			//loop back through and set BreakingWaveHeight and Direction in the shadow
			while (j < ShadowStart-1)
			{
				//calculate shadow zone angle to node SUPPLYING sediment
				dX = X[j+1]-X[ShadowStart];
				dY = Y[j+1]-Y[ShadowStart];
				SpatialGradient = dY/dX;
							
				//convert to azimuths and determine refracted diffracted wave properties
				if (dX == 0 && dY < 0) AngleDiff = ShadowAngle-180.;
				else if (dX == 0 && dY > 0) AngleDiff = ShadowAngle-360;
				else if (dX > 0) AngleDiff = (180./M_PI)*(M_PI*0.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				else if (dX < 0) AngleDiff = (180./M_PI)*(M_PI*1.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				//if (1.5*AngleDiff < -90.) TempWaveHeight = 0;
				//else
				TempWaveHeight = OffshoreWaveHeight*0.5*(1-sin((M_PI/180.)*(-AngleDiff)));

				if (TempWaveHeight < ShadowZoneWaveHeight[j] || ShadowZoneWaveHeight[j] == -9999)
				{
					ShadowZoneWaveHeight[j] = TempWaveHeight;
					ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*AngleDiff;
				}
				++j;
			}
			
			//Reduce wave heights outside the shadow zone over the length that they were redued in the shadow zone
			//This should preserve energy.
			k = ShadowEnd;
			while(Dist < ShadowZoneDistance)
			{
				ShadowZoneWaveHeight[k] = OffshoreWaveHeight*(1.-0.5*(1.-sin((M_PI/180.)*(90.*(Dist/ShadowZoneDistance)))));
				ShadowZoneWaveDirection[k] = OffshoreWaveDirection;
				Shadows[k+1] = 1;
				k-=1;
				if (Shadows[k+1] != 0 || k<2 || k >NoNodes-3) Dist = ShadowZoneDistance+1;
				else Dist += (Distance[k+1]-Distance[k]);
			}
		}
	}
}

void Coastline::IntersectionAnalysis()
{
	/*description goes here
	
	Function to loop through the coastline vector and check for interesctions due 
	to shoreline meeting up or coast eating itself.
	
	Finds the cross product between two line segments and the line connecting
	the start of the two segments in order to determine if they intersect
	
	No logic to deal with this yet, just a detection flag
	
	Martin Hurst, 28/1/2014
	
	*/
	
	//declare temporary variables
	double X1, X2, X3, X4, Y1, Y2, Y3, Y4, dX12, dY12, dX34, dY34, dX31, dY31, XProd, S, T;
	int XProdPos = 0;
	
	for (int i=0; i<NoNodes-2; ++i)
	{
		//assign first segment
		X1 = X[i];
		Y1 = Y[i];
		X2 = X[i+1];
		Y2 = Y[i+1];
		dX12 = X2-X1;
		dY12 = Y2-Y1;
		
		for (int j=i+2; j<NoNodes-1; ++j)
		{
			//assign second segment
			X3 = X[j];
			Y3 = Y[j];
			X4 = X[j+1];
			Y4 = Y[j+1];
			dX34 = X4-X3;
			dY34 = Y4-Y3;
			
			//Find the cross product of the two vectors
			XProd = dX12*dY34 - dX34*dY12;
			
			if (XProd == 0) {} 	//lines are colinear/parallel
			else
			{
				if (XProd > 0) XProdPos = 1;
				else XProdPos = 0;
				
				//assign third test segment
				dX31 = X1-X3;
				dY31 = Y1-Y3;
				
				//get cross products
				S = dX12*dY31 - dY12*dX31;
				T = dX34*dY31 - dY34*dX31;
				
				//logic for collision occurence
				if ((S < 0) == XProdPos) {}
				else if ((T < 0) == XProdPos) {}
				else if (((S > XProd) == XProdPos) || ((T > XProd) == XProdPos)) {}
				else
				{
					if ((Fixed[i] == 0) && (Fixed[j] == 0))
					{	
						cout << "Intersection in coastline detected!!" << endl;
						cout << "i = " << i << "; j = " << j << "; NoNodes = " << NoNodes << endl;
					}

					if (i==0) Fixed[j] = 1;
					else if (j==NoNodes-2) Fixed[i] = 1;
					else
					{
						X.erase(X.begin()+i+1, X.begin()+j+1);
						Y.erase(Y.begin()+i+1, Y.begin()+j+1);
						Fixed.erase(Fixed.begin()+i+1, Fixed.begin()+j+1);
						NoNodes -= (j-i);
					}
				}
			}
		}
	}
	CheckNodeSpacing();
	CalculateMorphology();
}

void Coastline::SetupQueue()
{
	/* DESCRIPTION GOES HERE
	
	Function to find shoreface triangles and create priority queue for computing sediment
	fluxes and updating the coastline_CPP
	
	Expand this to populate a vector of effective shoreface depths?
	
	Martin Hurst, March 2014
	
	*/
	
	CoastNode CurrentNode;
	double Wb, D_sf;
	
	//Loop through coast and look for triangles
	for (int i=2; i<NoNodes-2; ++i)
	{
		//Get cell width at bottom of the shoreface
		Wb = CellWidth[i] + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));	
			
		//Flag triangles and calculate mass of sand assuming a fixed beach width
		if (Wb < 0) 
		{
			//Flag triangles
			FlagTriangles[i] = 1;
							
			//Get effective shoreface depth
			D_sf = -(CellWidth[i]*ShorefaceSlope)/(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
			
			CurrentNode.ShorefaceDepth = D_sf;
			CurrentNode.i = i;
			CoastlineQueue.push(CurrentNode);
		}
		else 
		{
			CurrentNode.ShorefaceDepth = ClosureDepth;
			CurrentNode.i = i;
			CoastlineQueue.push(CurrentNode);
		}
	}
}

void Coastline::BuildCellGeometries()
{
	/* DESCRIPTION GOES HERE
	
	Function to build pseudo-bathymetry in order to get depths over which to calculate sediment flux
	
	*/
	
	//First setup a priority queue based on ShorefaceDepth
	CoastNode CurrentNode;
	double Wb, D_sf, ShorefaceLength;
	double O1, O2, M1, M2, dX, dY, O;
	int i;
	vector<double> MeshOrientation = Orientation;
	vector<double> Tracker(NoNodes,1);
	
	//Loop through coast and look for triangles
	for (int i=2; i<NoNodes-2; ++i)
	{
		//Get cell width at bottom of the shoreface
		Wb = CellWidth[i] + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));	
			
		//Flag triangles and calculate mass of sand assuming a fixed beach width
		if (Wb < 0) 
		{
			//Flag triangles
			FlagTriangles[i] = 1;
							
			//Get effective shoreface depth
			D_sf = -(CellWidth[i]*ShorefaceSlope)/(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
			
			CurrentNode.ShorefaceDepth = D_sf;
			CurrentNode.i = i;
			CoastlineQueue.push(CurrentNode);
		}
		else 
		{
			CurrentNode.ShorefaceDepth = ClosureDepth;
			CurrentNode.i = i;
			CoastlineQueue.push(CurrentNode);
		}
	}
	
	//loop through the queue and calculate new shoreface depth, adding adjacent nodes
	//to the queue as we go
	while (!CoastlineQueue.empty())
	{
		//pop a node out of the priority queue
		CoastNode CurrentNode = CoastlineQueue.top();
		i = CurrentNode.i;
		D_sf = CurrentNode.ShorefaceDepth;
		CoastlineQueue.pop();
		
		if (i<2 || i>NoNodes-3 || D_sf<=0 || Tracker[i]==0) continue;
		else if (Tracker[i] == 1 || Tracker[i] == 4)
		{
			if (D_sf < ClosureDepth)
			{
				//Find intersection point between cell boundaries
				O1 = Orientation[i]-90.-e1[i];
        		O2 = Orientation[i]-90.+e2[i];
        		M1 = 1./tan((M_PI/180.)*O1);
        		M2 = 1./tan((M_PI/180.)*O2);
            	X0[i] = (YMidPoints[i] - M2*XMidPoints[i] + M1*XMidPoints[i-1] - YMidPoints[i-1])/(M1 - M2);
        		Y0[i] = M1*X0[i] + YMidPoints[i-1] - M1*XMidPoints[i-1];
        
        		//get Dsf
        		//may need to multiply by cos angle diff
        		dX = X0[i]-X[i];
        		dY = Y0[i]-Y[i];
        		if (dX == 0 && dY < 0) O = 180.;
				else if (dX == 0 && dY > 0) O = 0.;
				else if (dX > 0) O = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
				else if (dX < 0) O = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
				else cout << "es ist unmÃ¶glich" << endl;
			
        		ShorefaceLength = sqrt(pow(X0[i]-X[i],2.) + pow(Y0[i]-Y[i],2.))*cos((M_PI/180.)*(Orientation[i]-90.-O));
				ShorefaceDepth[i] = ShorefaceLength*ShorefaceSlope;
        		        
        		Tracker[i] = 0;
        		Tracker[i+1] = 2;
        		Tracker[i-1] = 3;
        
        		MeshOrientation[i+1] = (MeshOrientation[i-1]+MeshOrientation[i]+MeshOrientation[i+1])/3.;
        		MeshOrientation[i-1] = (MeshOrientation[i-1]+MeshOrientation[i]+MeshOrientation[i+1])/3.;
    
    			CurrentNode.ShorefaceDepth = ShorefaceDepth[i];
				CurrentNode.i = i+1;
				CoastlineQueue.push(CurrentNode);
				CurrentNode.i = i-1;
				CoastlineQueue.push(CurrentNode);
			}
			else ShorefaceDepth[i] = ClosureDepth;
		}

		else if (Tracker[i] == 2)
		{
			//If tracker is 2 then i-1 is/was a triangle
			//Find intersection
			O1 = MeshOrientation[i]-90.;
			O2 = Orientation[i]-90.+e2[i];
			M1 = 1./tan((M_PI/180.)*O1);
        	M2 = 1./tan((M_PI/180.)*O2);
			X0[i] = (YMidPoints[i] - M2*XMidPoints[i] + M1*X0[i-1] - Y0[i-1])/(M1 - M2);
			Y0[i] = M1*X0[i] + Y0[i-1] - M1*X0[i-1];
			
			//get equivalent shoreface depth
			//may need to multiply by cos angle diff
        	dX = X0[i]-X[i];
        	dY = Y0[i]-Y[i];
        	if (dX == 0 && dY < 0) O = 180.;
			else if (dX == 0 && dY > 0) O = 0.;
			else if (dX > 0) O = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
			else if (dX < 0) O = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
			else cout << "es ist unmÃ¶glich" << endl;
			
        	ShorefaceLength = sqrt(pow(X0[i]-X[i],2.) + pow(Y0[i]-Y[i],2.))*cos((M_PI/180.)*(Orientation[i]-90.-O));
			ShorefaceDepth[i] = ShorefaceLength*ShorefaceSlope;
			
			if (ShorefaceDepth[i] < ClosureDepth)
			{
            	if (Tracker[i+1] == 1) Tracker[i+1] = 2;
	            else if (Tracker[i+1] == 3) cout << "This might cause trouble" << endl;
    	        CurrentNode.ShorefaceDepth = ShorefaceDepth[i];
				CurrentNode.i = i+1;
				CoastlineQueue.push(CurrentNode);
				MeshOrientation[i+1] = MeshOrientation[i];
        	}
        	else
        	{
            	if (Tracker[i+1] != 0) Tracker[i+1] = 4;
	            ShorefaceDepth[i] = ClosureDepth;
	            ShorefaceLength = (ShorefaceDepth[i]/ShorefaceSlope)/cos((M_PI/180.)*e2[i]);
    	        X0[i] = XMidPoints[i]-ShorefaceLength*cos((M_PI/180.)*Orientation[i]+e2[i]);
    	        Y0[i] = YMidPoints[i]+ShorefaceLength*sin((M_PI/180.)*Orientation[i]+e2[i]);
    		}
        	Tracker[i] = 0;
		}
		else if (Tracker[i] == 3)
		{
			//get gradients of intersecting lines
			O1 = MeshOrientation[i]-90.;
			O2 = Orientation[i]-90.-e1[i];
			M1 = 1./tan((M_PI/180.)*O1);
        	M2 = 1./tan((M_PI/180.)*O2);
			X0[i] = (YMidPoints[i-1] - M2*XMidPoints[i-1] + M1*X0[i+1] - Y0[i+1])/(M1 - M2);
			Y0[i] = M1*X0[i] + Y0[i+1] - M1*X0[i+1];
        
			//get equivalent shoreface depth
			//may need to multiply by cos angle diff
    		dX = X0[i]-X[i];
    		dY = Y0[i]-Y[i];
    		if (dX == 0 && dY < 0) O = 180.;
			else if (dX == 0 && dY > 0) O = 0.;
			else if (dX > 0) O = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
			else if (dX < 0) O = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
			else cout << "es ist unmÃ¶glich" << endl;
		
    		ShorefaceLength = sqrt(pow(X0[i]-X[i],2.) + pow(Y0[i]-Y[i],2.))*cos((M_PI/180.)*(Orientation[i]-90.-O));
			ShorefaceDepth[i] = ShorefaceLength*ShorefaceSlope;

        	if (ShorefaceDepth[i] < ClosureDepth)
        	{ 
         	   if (Tracker[i-1] == 1) Tracker[i-1] = 3;
            	else if (Tracker[i-1] == 2) cout << "This might cause trouble" << endl;
            	CurrentNode.ShorefaceDepth = ShorefaceDepth[i];
				CurrentNode.i = i-1;
				CoastlineQueue.push(CurrentNode);
				MeshOrientation[i-1] = MeshOrientation[i];
			}
			else
			{   
		        if (Tracker[i+1] != 0) Tracker[i+1] = 4;
		        ShorefaceDepth[i] = ClosureDepth;
		        ShorefaceLength = (ClosureDepth/ShorefaceSlope)/cos((M_PI/180.)*e1[i]);
		        X0[i] = XMidPoints[i-1]-ShorefaceLength*cos((M_PI/180.)*Orientation[i]-e1[i]);
		        Y0[i] = YMidPoints[i-1]+ShorefaceLength*sin((M_PI/180.)*Orientation[i]-e1[i]);
			}
			Tracker[i] = 0;
		}
	}
}


void Coastline::TransportSediment(double &TimeDelta, Wave TheWave, int FluxType, int RefDiffFlag, double FluxFraction)
{
	/* DESCRIPTION GOES HERE
	
	Function to calculate longshore sediment flux and move sediment along the coast
	
	*/

	//Declare temporary variables
	double G, dt;
	double PositionChange;
	int i; //iterators
	
	//get wave properties to individual variables
	OffshoreWavePeriod = TheWave.Get_WavePeriod();
	OffshoreWaveHeight = TheWave.Get_WaveHeight();
	OffshoreWaveDirection = TheWave.Get_WaveDirection();
	
	//handle offshore waves
	if ((OffshoreWaveDirection < (FluxOrientation[0]-170.)) || (OffshoreWaveDirection > (FluxOrientation[NoNodes-1]-10.))) return;

	//Find Shadow
	GetShadows();
	
	//refraction and Diffraction
	if (ShadowFlag == 1 && RefDiffFlag == 1) RefractDiffractShadowZone();
	
	//Transform waves
	TransformWaves();
	
	BuildCellGeometries();
	
	//Calculate Stable Timestep, Courant criterion following Dean and Dalrymple 2001
	for (i=0; i<NoNodes; ++i)
	{
		G = (0.41*pow(BreakingWaveHeight[i],2.5))/(ClosureDepth);
		if (G == G) dt = (0.66*MeanNodeSpacing)*(0.66*MeanNodeSpacing)/(2*G);
		if (dt < TimeDelta) 
		{
			TimeDelta = dt;
			cout << "TimeDelta = " << TimeDelta << endl;
		}
	}
	
	//Compute Sediment Fluxes and move coast
	//should now be able to do this in a single loop
	for (i=0; i<NoNodes; ++i)
	{
		//CALCULATE LONGSHORE FLUX
		CalculateFlux(i, FluxType);
		
		//CALCULATE VOLUME CHANGE
		//Handle Boundary Conditions
		if (i <= 1 || i>=NoNodes-3)
		{
			//PERIODIC BOUNDARY
			if (StartBoundary == 1 && EndBoundary == 1)
			{
				if (i == 0)
				{
					//Sediment Flux for last node
					CalculateFlux(NoNodes-2, FluxType);
					VolumeChange[i] = (LongshoreFlux[NoNodes-2]-LongshoreFlux[i])*(1.-LostFluxFraction)*TimeDelta;
				}
				else if (i == NoNodes-1) VolumeChange[i] = VolumeChange[0];
			}
			
			//FIXED BOUNDARY (with/without flux)
			else if (StartBoundary == 2 && EndBoundary == 2)
			{
				CalculateFlux(NoNodes-1, FluxType);
				if (i < 2)
				{
					if (LongshoreFlux[i] > 0) LongshoreFlux[i] = LongshoreFlux[0]*FluxFraction;
					else if (Shadows[i] == 2) LongshoreFlux[i] = 0;
					VolumeChange[i] = 0;
				}
				else if (i>NoNodes-3)
				{
					if (LongshoreFlux[i] < 0) LongshoreFlux[i] = LongshoreFlux[NoNodes-1]*FluxFraction;
					VolumeChange[i] = 0;
				}
				else if (i==NoNodes-3)
				{
					if (LongshoreFlux[i] < 0) LongshoreFlux[i] = LongshoreFlux[NoNodes-1]*FluxFraction;
					else if (Shadows[i] == 4) LongshoreFlux[i] = 0;
					VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction)*TimeDelta;
				}
			}
		}
		//all other cells
		else if (Fixed[i] == 1)
		{
			LongshoreFlux[i] = 0;
			VolumeChange[i] = 0;
		}
		else 
		{
			//volume change
			//this may not be needed after triangles are handled
			VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction)*TimeDelta;
		}
		
		//MOVE THE COAST
		//check if a triangle
		if (ShorefaceDepth[i] < ClosureDepth)
		{
			//Modify volume change for cell depth
			VolumeChange[i] = VolumeChange[i]*(ShorefaceDepth[i]/ClosureDepth);
			
			//Call cubic solver
			if (VolumeChange[i] != 0) PositionChange = SolveCubic(i, ShorefaceDepth[i]);
			else PositionChange = 0;
			
			if (fabs(PositionChange) > 10.)
			{
				cout << "Break" << endl;
			}
			
			//update coastal position
			X[i] -= PositionChange*cos((M_PI/180.)*Orientation[i]);
			Y[i] += PositionChange*sin((M_PI/180.)*Orientation[i]);
		}
		else
		{
			//solve quadratic equation to get position change
			//if cell is rectilinear so dont need to solve the quadratic
			if (VolumeChange[i] != 0) PositionChange = SolveQuadratic(i, ShorefaceDepth[i]);
			else PositionChange = 0;
			
			if (fabs(PositionChange) > 10.)
			{
				cout << "Break" << endl;
			}

			if (PositionChange != PositionChange)
			{
				//Check for NaN in PositionChange
				cout << "I am not equal to myself!!" << endl;
				cout << "Found NaN in Position Change: Node is " << i << endl;
			}

			//Decompose Position change into X and Y components
			X[i] -= PositionChange*cos((M_PI/180.)*Orientation[i]);
			Y[i] += PositionChange*sin((M_PI/180.)*Orientation[i]);
		}
	}
	//Check for intersections and check node spacing
	IntersectionAnalysis();
}

//	//MOVE THE COASTLINE
//	//replace this loop with a priority queue starting with most acute triangles?
//	//cout << "New LOOP!" << endl;
//	while (!CoastlineQueue.empty())
//	{
//		//pop a node out of the priority queue
//		CoastNode CurrentNode = CoastlineQueue.top();
//		i = CurrentNode.i;
//		D_sf = CurrentNode.ShorefaceDepth;
//		CoastlineQueue.pop();
//	
//		//check the node hasn't already been dealt with
//		if (EvolveFlag[i] != 1)
//		{
//			//check if a triangle
//			if (D_sf < ClosureDepth)
//			{
//				//first evolve the triangle cell
//				//Reduce D_sf slightly to allow quadratic solution
//				//TempVolumeChange = 0;
//				PositionChange = 0;
//				
//				//Calculate Volume Change
//				VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction)*TimeDelta;
//				
//				//Call cubic solver
//				if (VolumeChange[i] != 0) PositionChange = SolveCubic(i, ShorefaceDepth[i]);
//				
//				//update coastal position
//				X[i] -= PositionChange*cos((M_PI/180.)*Orientation[i]);
//				Y[i] += PositionChange*sin((M_PI/180.)*Orientation[i]);
//				
//				//mark cell as dealt with
//				EvolveFlag[i] = 1;
//			}
//		}
//		//else if not a triangle evolve in the normal way			
//		else
//		{
//			//solve quadratic equation to get position change
//			//if cell is rectilinear so dont need to solve the quadratic
//			if (VolumeChange[i] != 0) PositionChange = SolveQuadratic(i, ShorefaceDepth[i]);
//			else PositionChange = 0;

//			if (PositionChange != PositionChange)
//			{
//				//Check for NaN in PositionChange
//				cout << "I am not equal to myself!!" << endl;
//				cout << "Found NaN in Position Change: Node is " << i << endl;
//			}

//			//Decompose Position change into X and Y components
//			//Going to need a bit more logic here to handle different orientations?
//			//current setup x change is negative, y change is positive
//			X[i] -= PositionChange*cos((M_PI/180.)*Orientation[i]);
//			Y[i] += PositionChange*sin((M_PI/180.)*Orientation[i]);
//		
//			//Keep track of nodes dealt with
//			EvolveFlag[i] = 1;
//		}
//	}
//	
//	//Track Volume
//	//TotalVolumeOld = TotalVolume;
//	//VolumeOut = LongshoreFlux[NoNodes-3]*(1.-LostFluxFraction)*TimeDelta;;
//	
//	//Check for intersections and check node spacing
//	IntersectionAnalysis();
//}

////Try using 3cell orientation to get flux @i-1 and i?
//					if ((i > 2) && (i <NoNodes-2))
//					{
//						for (int j=i-1; j<(i+1); ++j)
//						{
//							e1[j] = Orientation[j]-FluxOrientation[j-1];
//							e2[j] = FluxOrientation[j]-Orientation[j];
//							FluxOrientation[j] = Orientation[j];
//							TransformWave(j);
//							CalculateFlux(j, FluxType);
//							VolumeChange[j] = (LongshoreFlux[j-1]-LongshoreFlux[j])*(1.-LostFluxFraction)*TimeDelta;
//							c = -VolumeChange[j]/ClosureDepth;
//							if (VolumeChange[j] != 0)
//							{
//								if (a==0) PositionChange = VolumeChange[j]/(CellWidth[j]*ClosureDepth);
//								//quadratic solver
//								else PositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//							}
//							if (PositionChange != PositionChange)
//							{
//								cout << "I am not equal to myself!!" << endl;
//								cout << "Found NaN in Position Change: Node is " << j << endl;
//								VolumeChange[j] = 0;
//								LongshoreFlux[j] = LongshoreFlux[j-1];
//								PositionChange = 0;
//							}
//						}
//					}		
//				m = i;
//				n = i;
//				
//				//spread outward to deal with adjacent cells
//				//m and n are indices for the spreading (m is left, n is right)
//				TriangleFlag = 1;
//				while (TriangleFlag == 1)
//				{
//					if (m > 2) --m;
//					if (n < NoNodes-3) ++n;
//					
//					//deal with m
//					if (EvolveFlag[m] == 1 && EvolveFlag[n] == 1) break;
//					if (EvolveFlag[m] != 1)
//					{
//						//reset vol and position change holders
//						TempVolumeChange = 0;
//						TempPositionChange = 0;
//						
//						//Calculate Volume Change
//						TempVolumeChange = (D_sf/ClosureDepth)*(LongshoreFlux[m-1]-LongshoreFlux[m])*(1.-LostFluxFraction)*TimeDelta;
//					
//						//quadratic coefficients
//						//note use of effective shoreface depth D_sf rather than "closure depth"
//						a = tan((M_PI/180.)*e1[m]) + tan((M_PI/180.)*e2[m]);
//						b = 2.*CellWidth[m] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[m]) + tan((M_PI/180.)*e2[m]));
//						c = -2.*TempVolumeChange/D_sf;
//					
//						//solve quadratic equation to get position change
//						//if cell is rectilinear so dont need to solve the quadratic
//						if (TempVolumeChange != 0)
//						{
//							if (fabs(a)<0.0001) TempPositionChange = TempVolumeChange/(CellWidth[m]*D_sf);
//							else TempPositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//						}
//					
//						X[m] -= TempPositionChange*cos((M_PI/180.)*Orientation[m]);
//						Y[m] += TempPositionChange*sin((M_PI/180.)*Orientation[m]);
//				
//						//mark cell as dealt with
//						EvolveFlag[m] = 1;
//					}
//					if (EvolveFlag[n] != 1)
//					{
//						//deal with n
//						//reset vol and position change holders
//						TempVolumeChange = 0;
//						TempPositionChange = 0;
//						
//						//Calculate Volume Change
//						TempVolumeChange = (D_sf/ClosureDepth)*(LongshoreFlux[n-1]-LongshoreFlux[n])*(1.-LostFluxFraction)*TimeDelta;
//					
//						//quadratic coefficients
//						//note use of effective shoreface depth D_sf rather than "closure depth"
//						a = tan((M_PI/180.)*e1[n]) + tan((M_PI/180.)*e2[n]);
//						b = 2.*CellWidth[n] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[n]) + tan((M_PI/180.)*e2[n]));
//						c = -2.*TempVolumeChange/D_sf;
//					
//						//solve quadratic equation to get position change
//						//if cell is rectilinear so dont need to solve the quadratic
//						if (TempVolumeChange != 0)
//						{
//							if (fabs(a)<0.0001) TempPositionChange = TempVolumeChange/(CellWidth[n]*D_sf);
//							else TempPositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//						}
//					
//						X[n] -= TempPositionChange*cos((M_PI/180.)*Orientation[n]);
//						Y[n] += TempPositionChange*sin((M_PI/180.)*Orientation[n]);
//				
//						//mark cell as dealt with
//						EvolveFlag[n] = 1;
//					}
//					
//					//Get Orientation simple version for now might change this
//					//Next get curvature and orientation (3-node)
//					dX = (X[n+1]-X[m-1]);
//					dY = (Y[n+1]-Y[m-1]);
//		
//					//convert to azimuths
//					if (dX == 0 && dY < 0) TempOrientation = 180.;
//					else if (dX == 0 && dY > 0) TempOrientation = 0.;
//					else if (dX > 0) TempOrientation = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
//					else if (dX < 0) TempOrientation = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
//				
//					//Cell Width parallel to orientation AT the node position // this needs changing
//					TempCellWidth = (0.5*sqrt(pow(X[i]-X[m-1],2.0) + pow(Y[i]-Y[m-1],2.0)))/cos(M_PI/180.*(TempOrientation-FluxOrientation[m-1]));
//					TempCellWidth += (0.5*sqrt(pow(X[n+1]-X[i],2.0) + pow(Y[n+1]-Y[i],2.0)))/cos(M_PI/180.*(TempOrientation-FluxOrientation[n]));
//		
//					//Get new shoreface depth based on edges of the two cells m and n
//					//Get cell width at bottom of the shoreface
//					Temp_e1 = TempOrientation-FluxOrientation[m-1];
//					Temp_e2 = FluxOrientation[n] - TempOrientation;
//				
//					//Check Width at shoreface depth
//					Wb = TempCellWidth + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2));	
//					if (Wb > 0) 
//					{
//						TriangleFlag = 0;
//						D_sf = ClosureDepth;
//					}
//					else
//					{
//						D_sf_new = -(TempCellWidth*ShorefaceSlope)/(tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2)) - 0.2;
//						if (D_sf_new > D_sf) 
//						{
//							TriangleFlag = 0;
//							D_sf = D_sf_new;
//						}
//					}
//				}
//			}
			
			
		
//		//Deal with triangular cells first, treat as larger cells so as to still use quadratic solver
//		
//		//check the node hasn't already been dealt with
//		if (EvolveFlag[i] != 1)
//		{
//			if (D_sf < ClosureDepth)
//			{
//				//Deal with triangle cells and their surroundings here
//				//a and b are indexes to the edges of the larger trapezium
//				m = i;
//				n = i;
//				
//				TriangleFlag = 1;
//				while (TriangleFlag == 1)
//				{

//					///Expand to find smallest larger cell that can be a trapezium
//					//stop trying to make trapezium bigger at boundaries
//					if ((EvolveFlag[m] == 1) || (EvolveFlag[n] == 1)) break;
//					
//					//Get Orientation simple version for now might change this
//					//Next get curvature and orientation (3-node)
//					dX = (X[n+1]-X[m-1]);
//					dY = (Y[n+1]-Y[m-1]);
//		
//					//convert to azimuths
//					if (dX == 0 && dY < 0) TempOrientation = 180.;
//					else if (dX == 0 && dY > 0) TempOrientation = 0.;
//					else if (dX > 0) TempOrientation = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
//					else if (dX < 0) TempOrientation = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
//				
//					//Cell Width parallel to orientation AT the node position // this needs changing
//					if (m==2) TempCellWidth = (0.5*sqrt(pow(X[i]-X[m-1],2.0) + pow(Y[i]-Y[m-1],2.0)))/cos(M_PI/180.*(TempOrientation-90.));
//					else TempCellWidth = (0.5*sqrt(pow(X[i]-X[m-1],2.0) + pow(Y[i]-Y[m-1],2.0)))/cos(M_PI/180.*(TempOrientation-FluxOrientation[m-1]));
//					TempCellWidth += (0.5*sqrt(pow(X[n+1]-X[i],2.0) + pow(Y[n+1]-Y[i],2.0)))/cos(M_PI/180.*(TempOrientation-FluxOrientation[n]));
//		
//					//Get cell width at bottom of the shoreface
//					Temp_e1 = TempOrientation-FluxOrientation[m-1];
//					Temp_e2 = FluxOrientation[n] - TempOrientation;
//				
//					//Check Width at shoreface depth
//					Wb = TempCellWidth + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2));	
//				
//					//Flag triangles and calculate mass of sand assuming a fixed beach width
//					if (Wb > 0) 
//					{
//						TriangleFlag = 0;
//						D_sf = ClosureDepth;
//					}
//					else if (fabs(Orientation[n]-Orientation[m]) > 30.) 
//					{
//						TriangleFlag = 0;
//						D_sf = -(TempCellWidth*ShorefaceSlope)/(tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2)) - 0.5;
//						if (D_sf < 0) cout << "FUCK!" << endl;
//					}
//					else
//					{
//						if (EvolveFlag[m-1] == 1 && EvolveFlag[n+1] == 1) break;
//						if ((m > 2) && (EvolveFlag[m-1] != 1)) m-=1;
//						if ((n < NoNodes-3) && (EvolveFlag[n+1] != 1)) n+=1;
//					}
//				}
//				
//				//VolumeChange of large cell
//				TempVolumeChange = (LongshoreFlux[m-1]-LongshoreFlux[n])*(1.-LostFluxFraction)*TimeDelta;
//					
//				//move the coast
//				//quadratic coefficients
//				a = tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2);
//				b = 2.*TempCellWidth + (tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2));
//				c = -2.*TempVolumeChange/D_sf;

//				//solve quadratic equation to get position change
//				//if cell is rectilinear dont need to solve the quadratic
//				if (TempVolumeChange != 0)
//				{
//					if (a==0) TempPositionChange = TempVolumeChange/(TempCellWidth*ClosureDepth);
//					else TempPositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//				}
//				
//				//update position of all nodes in large cell
//				for (int k=m; k<n+1; ++k)
//				{
//					X[k] -= TempPositionChange*cos((M_PI/180.)*TempOrientation);
//					Y[k] += TempPositionChange*sin((M_PI/180.)*TempOrientation);
//				}
//				
//				D_sf = CurrentNode.ShorefaceDepth-0.5;
//	
//				//LOOP THROUGH NODES IN LARGE CELL AND ADJUST SHORE
//				for (int k=m; k<n+1; ++k)
//				{
//					//left boundary becomes no flux since we've dealt with it already
//					if (k==m) TempLongshoreFluxLeft = 0;
//					else TempLongshoreFluxLeft = LongshoreFlux[k-1]*(D_sf/ClosureDepth);
//					
//					//right boundary becomes no flux since we've already dealt with it
//					if (k==n) TempLongshoreFluxRight = 0;
//					else TempLongshoreFluxRight = LongshoreFlux[k]*(D_sf/ClosureDepth);
//					
//					//Calculate Volume Change
//					TempVolumeChange = (TempLongshoreFluxLeft-TempLongshoreFluxRight)*(1.-LostFluxFraction)*TimeDelta;
//					
//					//quadratic coefficients
//					//note use of effective shoreface depth D_sf rather than "closure depth"
//					a = tan((M_PI/180.)*e1[k]) + tan((M_PI/180.)*e2[k]);
//					b = 2.*CellWidth[k] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[k]) + tan((M_PI/180.)*e2[k]));
//					c = -2.*TempVolumeChange/D_sf;

//					//solve quadratic equation to get position change
//					//if cell is rectilinear so dont need to solve the quadratic
//					if (VolumeChange[k] != 0)
//					{
//						if (a==0) TempPositionChange = TempVolumeChange/(CellWidth[k]*D_sf);
//						else TempPositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//					}
//					
//					X[k] -= TempPositionChange*cos((M_PI/180.)*Orientation[k]);
//					Y[k] += TempPositionChange*sin((M_PI/180.)*Orientation[k]);
//					
//					EvolveFlag[k] = 1;
//				}	
//			}
//			
//			//else evolve normally
//			else
//			{
//				//quadratic coefficients
//				a = tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]);
//				b = 2.*CellWidth[i] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
//				c = -2.*VolumeChange[i]/D_sf;

//				//solve quadratic equation to get position change
//				//if cell is rectilinear so dont need to solve the quadratic
//				if (VolumeChange[i] != 0)
//				{
//					if (a==0) PositionChange = VolumeChange[i]/(CellWidth[i]*D_sf);
//					else PositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//				}
//				else PositionChange = 0;

//				if (PositionChange != PositionChange)
//				{
//					//Check for NaN in PositionChange
//					cout << "I am not equal to myself!!" << endl;
//					cout << "Found NaN in Position Change: Node is " << i << endl;
//				
//					//Try using 3cell orientation to get flux @i-1 and i?
//					if ((i > 2) && (i <NoNodes-2))
//					{
//						for (int j=i-1; j<(i+1); ++j)
//						{
//							e1[j] = Orientation[j]-FluxOrientation[j-1];
//							e2[j] = FluxOrientation[j]-Orientation[j];
//							FluxOrientation[j] = Orientation[j];
//							TransformWave(j);
//							CalculateFlux(j, FluxType);
//							VolumeChange[j] = (LongshoreFlux[j-1]-LongshoreFlux[j])*(1.-LostFluxFraction)*TimeDelta;
//							c = -VolumeChange[j]/ClosureDepth;
//							if (VolumeChange[j] != 0)
//							{
//								if (a==0) PositionChange = VolumeChange[j]/(CellWidth[j]*ClosureDepth);
//								//quadratic solver
//								else PositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//							}
//							if (PositionChange != PositionChange)
//							{
//								cout << "I am not equal to myself!!" << endl;
//								cout << "Found NaN in Position Change: Node is " << j << endl;
//								VolumeChange[j] = 0;
//								LongshoreFlux[j] = LongshoreFlux[j-1];
//								PositionChange = 0;
//							}
//						}
//					}
//				}

//				//Decompose Position change into X and Y components
//				//Going to need a bit more logic here to handle different orientations?
//				//current setup x change is negative, y change is positive
//				X[i] -= PositionChange*cos((M_PI/180.)*Orientation[i]);
//				Y[i] += PositionChange*sin((M_PI/180.)*Orientation[i]);
//			
//				//Keep track of nodes dealt with
//				EvolveFlag[i] = 1;
//			}
//		}
//	}
//	//Check for intersections and check node spacing
//	IntersectionAnalysis();
//}

	
//	for (i=0; i<NoNodes; ++i)
//	{				
//		if (FlagTriangles[i] == 1)
//		{
//			if (i == 2)
//			{
//				//Include more nodes and see if still a triangle
//				j = i+1;
//				TriangleFlag = 1;
//				while (TriangleFlag == 1)
//				{
//					//Get Orientation
//					TempOrientation = (Orientation[i] + Orientation[j])/2.;
//				
//					//Cell Width parallel to orientation AT the node position // this needs changing
//					TempCellWidth = (0.5*sqrt(pow(X[i]-X[i-1],2.0) + pow(Y[i]-Y[i-1],2.0)))/cos(M_PI/180.*(TempOrientation-FluxOrientation[i-1]));
//					TempCellWidth += (0.5*sqrt(pow(X[j+1]-X[j],2.0) + pow(Y[j+1]-Y[j],2.0)))/cos(M_PI/180.*(TempOrientation-FluxOrientation[j]));
//		
//					//Get cell width at bottom of the shoreface
//					Temp_e1 = TempOrientation-FluxOrientation[i-1];
//					Temp_e2 = FluxOrientation[j] - TempOrientation;
//				
//					//Check Width at shoreface depth
//					Wb = TempCellWidth + (ClosureDepth/ShorefaceSlope)*(tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2));	
//				
//					//Flag triangles and calculate mass of sand assuming a fixed beach width
//					if (Wb > 0) TriangleFlag = 0;
//					else ++j;
//				}
//				
//				//VolumeChange of largest cell
//				TempVolumeChange = (LongshoreFlux[i-1]-LongshoreFlux[j])*(1.-LostFluxFraction)*TimeDelta;
//					
//				//move the coast
//				//quadratic coefficients
//				a = tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2);
//				b = 2.*TempCellWidth + (tan((M_PI/180.)*Temp_e1) + tan((M_PI/180.)*Temp_e2));
//				c = -2.*TempVolumeChange/ClosureDepth;

//				//solve quadratic equation to get position change
//				//if cell is rectilinear so dont need to solve the quadratic
//				if (TempVolumeChange != 0)
//				{
//					if (a==0) TempPositionChange = TempVolumeChange/(TempCellWidth*ClosureDepth);
//					//quadratic solver
//					else TempPositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//				}
//				
//				//update position of all nodes in large cell
//				for (int k=i; k<j+1; ++k)
//				{
//					X[k] -= TempPositionChange*cos((M_PI/180.)*TempOrientation);
//					Y[k] += TempPositionChange*sin((M_PI/180.)*TempOrientation);
//				}
//					
//				//now deal with the shoreface
//				D_sf = 10.;
//				for (int k=i; k<j+1; ++k)
//				{
//					Temp_D_sf = -(CellWidth[i]*ShorefaceSlope)/(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i])) - 1.;
//					if (Temp_D_sf < D_sf) D_sf = Temp_D_sf;
//				}
//				
//				//LOOP THROUGH NODES IN LARGE CELL AND ADJUST SHORE
//				for (int k=i; k<j+1; ++k)
//				{
//					if (k==i) TempLongshoreFluxLeft = 0;
//					else TempLongshoreFluxLeft = LongshoreFlux[k-1]*(D_sf/ClosureDepth);
//					
//					if (k==j) TempLongshoreFluxRight = 0;
//					else TempLongshoreFluxRight = LongshoreFlux[k]*(D_sf/ClosureDepth);
//					
//					TempVolumeChange = (TempLongshoreFluxLeft-TempLongshoreFluxRight)*(1.-LostFluxFraction)*TimeDelta;
//					
//					//quadratic coefficients
//					a = tan((M_PI/180.)*e1[k]) + tan((M_PI/180.)*e2[k]);
//					b = 2.*CellWidth[i] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[k]) + tan((M_PI/180.)*e2[k]));
//					c = -2.*TempVolumeChange/D_sf;

//					//solve quadratic equation to get position change
//					//if cell is rectilinear so dont need to solve the quadratic
//					if (VolumeChange[k] != 0)
//					{
//						if (a==0) TempPositionChange = TempVolumeChange/(CellWidth[k]*D_sf);
//						//quadratic solver
//						else TempPositionChange = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//					}
//					
//					X[k] -= TempPositionChange*cos((M_PI/180.)*Orientation[k]);
//					Y[k] += TempPositionChange*sin((M_PI/180.)*Orientation[k]);
//					
//					//And deal with next cell too
//					TempVolumeChange = TempLongshoreFluxRight*(1.-LostFluxFraction)*TimeDelta;
//					
//					//update i to skip
//					i = k;
//				}
//				
//				
//			}									
//			else
//			{
//				//accumulation in concave wedge
//				cout << "Found Triangle" << endl;
//				if (VolumeChange[i] > 0) D_sf = -(CellWidth[i]*ShorefaceSlope)/(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i])) - 1.;
//				else D_sf = -(CellWidth[i]*ShorefaceSlope)/(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
//			}
//		}
//		else 
//		{
//			D_sf = ClosureDepth;
//				
//			//quadratic coefficients
//			a = tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]);
//			b = 2.*CellWidth[i] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
//			c = -2.*VolumeChange[i]/D_sf;

//			//solve quadratic equation to get position change
//			//if cell is rectilinear so dont need to solve the quadratic
//			if (VolumeChange[i] != 0)
//			{
//				if (a==0) PositionChange[i] = VolumeChange[i]/(CellWidth[i]*D_sf);
//				//quadratic solver
//				else PositionChange[i] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//			}
//			else PositionChange[i] = 0;

//			if (PositionChange[i] != PositionChange[i])
//			{
//				//Check for NaN in PositionChange
//				cout << "I am not equal to myself!!" << endl;
//				cout << "Found NaN in Position Change: Node is " << i << endl;
//				
//				//Try using 3cell orientation to get flux @i-1 and i?
//				if ((i > 2) && (i <NoNodes-2))
//				{
//					for (int j=i-1; j<(i+1); ++j)
//					{
//						e1[j] = Orientation[j]-FluxOrientation[j-1];
//						e2[j] = FluxOrientation[j]-Orientation[j];
//						FluxOrientation[j] = Orientation[j];
//						TransformWave(j);
//						CalculateFlux(j, FluxType);
//						VolumeChange[j] = (LongshoreFlux[j-1]-LongshoreFlux[j])*(1.-LostFluxFraction)*TimeDelta;
//						c = -VolumeChange[j]/ClosureDepth;
//						if (VolumeChange[j] != 0)
//						{
//							if (a==0) PositionChange[j] = VolumeChange[j]/(CellWidth[j]*ClosureDepth);
//							//quadratic solver
//							else PositionChange[j] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
//						}
//						if (PositionChange[j] != PositionChange[j])
//						{
//							cout << "I am not equal to myself!!" << endl;
//							cout << "Found NaN in Position Change: Node is " << j << endl;
//							VolumeChange[j] = 0;
//							LongshoreFlux[j] = LongshoreFlux[j-1];
//							PositionChange[j] = 0;
//						}
//					}
//				}
//			}

//			//Decompose Position change into X and Y components
//			//Going to need a bit more logic here to handle different orientations?
//			//current setup x change is negative, y change is positive
//			X[i] -= PositionChange[i]*cos((M_PI/180.)*Orientation[i]);
//			Y[i] += PositionChange[i]*sin((M_PI/180.)*Orientation[i]);
//		}
//	}
//	//Check for intersections and check node spacing
//	IntersectionAnalysis();
//}

void Coastline::CalculateFlux(int i, int FluxType)
{
	/*
	Seidment transport functions goes in here
	
	Type of transport dependent on specified flux type:
		0 = Simple Diffusion
		1 = CERC Diffusion
		2 = Kamphuis Diffusion
		3 = Bailard?
		4 = Deigaard?
		5 = Damgaard-Soulsby (Bedload/Shingle only)?
		6 = Bayram?
		7 = Morfett?

	86400 multiplier converts from units of m3/s to m3/day since timesteps in 
	model are in days.
	
	*/
	
	if (FluxType == 0)
	{
		//Simple Diffusion this will need fixing
		double DiffCoef = 1.0;
		LongshoreFlux[i] = DiffCoef*FluxOrientation[i];
		//cout << "FluxType doesn't exist yet" << endl;
		//exit(EXIT_FAILURE);
	}
	else if (FluxType == 1)
	{
		//CERC Equation
		LongshoreFlux[i] = 86400.*0.41
						*pow(BreakingWaveHeight[i],2.5)
						*sin((M_PI/180.0)*BreakingWaveAngle[i])
						*cos((M_PI/180.0)*BreakingWaveAngle[i]);
		
		if (LongshoreFlux[i] != LongshoreFlux[i])
		{
			//Check for NaN in PositionChange
			cout << "I am not equal to myself!!" << endl;
			cout << "Found NaN in Position Change: Node is " << i << endl;
		}
	}
	else if (FluxType == 11)
	{
		//CERC Equation, Ozasa & Brampton (1980) style! (see also van den Berg 2012).
		//Justification still unclear as noted by Ashton et al 2006.
		double dX = X[i+1]-X[i];
		double dY = Y[i+1]-Y[i];
		double Distance = sqrt(dX*dX + dY*dY);
		double K1 = 0.41;
		double K2 = 2.*K1/ShorefaceSlope;
		LongshoreFlux[i] = 86400.*(BreakingWaveHeight[i],2.5)
							*(K1*sin(2.*BreakingWaveAngle[i]) 
							- K2*cos(BreakingWaveAngle[i])*((BreakingWaveAngle[i+1]-BreakingWaveAngle[i])/Distance));
						
		if (LongshoreFlux[i] != LongshoreFlux[i])
		{
			//Check for NaN in PositionChange
			cout << "I am not equal to myself!!" << endl;
			cout << "Found NaN in Position Change: Node is " << i << endl;
		}
	}
	else if (FluxType == 2)
	{
		//Kamphius Equation
//		LongshoreFlux[i] = 86400.*2.27*pow(ShorefaceSlope,0.75)*pow(Period,1.5)
//								*pow(BreakingWaveHeight[i],2.0)
//								*pow(cos((M_PI/180.0)*BreakingWaveAngle[i]),0.6)
//								*pow(sin((M_PI/180.0)*BreakingWaveAngle[i]),0.6);
		cout << "FluxType doesn't exist yet" << endl;
		exit(EXIT_FAILURE);
	}
	else if (FluxType == 3)
	{
		//Bailard Equation
		//This aint ready yet, what are u_mb and W?
		//u_mb is the maximum oscillatory velocity of the breaking wave = (gamma/2)*sqrt(gh)
		//W is sediment fall velocity
		//gamma is 0.78
//		LongshoreFlux[i] = 86400.*(0.05 + 2.6*pow(sin((M_PI/180.0)*2.*BreakingWaveAngle[i]),2.) + 0.007*u_mb/FallVelocity)
//								*pow(BreakingWaveHeight[i],2.5)
//								*sin((M_PI/180.0)*BreakingWaveAngle[i])
//								*cos((M_PI/180.0)*BreakingWaveAngle[i]);
		cout << "FluxType doesn't exist yet" << endl;
		exit(EXIT_FAILURE);
	}
	else if (FluxType == 4)
	{
//		//Deigaard Equation
//		LongshoreFlux[i] = 86400.*MaxLongshoreFlux
//								*pow(sin((2.*(M_PI/180.0)*BreakingWaveAngle[i]))
//										*(1.-0.4*(BreakingWaveAngle[i]/90.)*(1.-BreakingWaveAngle[i]/90.)) ,2.5);
		cout << "FluxType doesn't exist yet" << endl;
		exit(EXIT_FAILURE);
	}
	else if (FluxType == 5)
	{
		cout << "FluxType doesn't exist yet" << endl;
		exit(EXIT_FAILURE);
	}
	else if (FluxType == 6)
	{
		cout << "FluxType doesn't exist yet" << endl;
		exit(EXIT_FAILURE);
	}
	else if (FluxType == 7)
	{
		cout << "FluxType doesn't exist yet" << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		cout << "FluxType not recognised" << endl;
		exit(EXIT_FAILURE);
	}
}

void Coastline::SimpleDiffusion(double DiffCoeff, double TimeDelta)
{
	//do diffusion
	//Diffuse coastline as a function of curvature
	//Coastline adjusts perpendicular to orientation
	//Calls CalculateMorphology() at the end to update

	//temp variables
	double PositionChange;

	for (int i=0; i<NoNodes; ++i)
	{
		if (i > 1 || i < NoNodes-2)
		{
			//Calculate Position Change
			PositionChange = -DiffCoeff*Curvature[i]*TimeDelta;
		
			//Decompose Position change into X and Y components
			//Going to need a bit more logic here to handle different orientations
			//current setup x change is negative, y change is positive
			X[i] -= PositionChange*cos((M_PI/180.)*Orientation[i]);
			Y[i] += PositionChange*sin((M_PI/180.)*Orientation[i]);
		}
	}

	//Update Coastal Morphology
}

double Coastline::SolveQuadratic(int i, double D_sf)
{
	//quadratic coefficients
	double a, b, c;
	//note use of effective shoreface depth D_sf rather than "closure depth"
	//D_sf may need to become a vector now
	a = tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]);
	b = 2.*CellWidth[i] + (D_sf/ShorefaceSlope)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
	c = -2.*VolumeChange[i]/D_sf;
	
	//return solution
	if (a == 0) return VolumeChange[i]/(CellWidth[i]*ClosureDepth);
	return (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
}

double Coastline::SolveCubic(int i, double D_sf)
{
	//Solve Cubic equation of the form ax^3 + bx^2 + cx + d = 0
	//following Press at al. 1992 for solution to positionchange
	//for cells which do not reach closure depth
	double a, b, c, d, a1, b1, c1, Q, Q3, R, R2, A, B, Theta, min2rootQ, aover3, x1, x2, x3;
	double D2 = D_sf;
	double D0 = D_sf;
	double epsilon = 0.00001;
	
	//setup coefficients within the function
	//if cell is triangular use these coefficients
	if (1 == 1)
	{
		a = -(1./3.)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
		b = -CellWidth[i]/2.;
		c = (CellWidth[i]*D_sf/ShorefaceSlope) - ((D_sf*D_sf)/(2.*ShorefaceSlope*ShorefaceSlope))*(tan((M_PI/180.)*-e1[i]) + tan((M_PI/180.)*-e2[i]));
		d = -VolumeChange[i]/ShorefaceSlope;
	}
	else
	{
		a = -(1./3.)*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
		b = -CellWidth[i]/2.;
		c = 0.5*CellWidth[i]*(D2-D0+D0/ShorefaceSlope)+(((D2-D0)*ShorefaceSlope-2.*D0*D0)/(4.*ShorefaceSlope*ShorefaceSlope))*(tan((M_PI/180.)*e1[i]) + tan((M_PI/180.)*e2[i]));
		d = -VolumeChange[i]/ShorefaceSlope;
	}
	
	if (fabs(a) < epsilon) return VolumeChange[i]/(CellWidth[i]*ClosureDepth);
	
	//convert to x^3 + ax^2 + bx + c = 0
	a1 = b/a;
	b1 = c/a;
	c1 = d/a;
	
	//Compute Q and R
	Q = (a1*a1 - 3.*b1) / 9.;
	R = (2*a1*a1*a1 - 9.*a1*b1 + 27.*c1) / 54.;
	
	//Check if there are 3 real roots
	R2 = R*R;
	Q3 = Q*Q*Q;
	aover3 = a1/3.;
	
	if (R2 < Q3)
	{
		//find the 3 roots
		Theta = acos(R/sqrt(Q3));
		min2rootQ = -2.*sqrt(Q);
		x1 = min2rootQ*cos(Theta/3.) - aover3;
		x2 = min2rootQ*cos((Theta+2.*M_PI)/3.) - aover3;
		x3 = min2rootQ*cos((Theta-2.*M_PI)/3.) - aover3;
		
		//need some logic here for what is real solution...?
		if (e1[i] < 0) 
		{
			if (x1 > 10. || x1 < -10.)
			{
				cout << "problem" << endl;
			}
			return x1;
		}
		else 
		{
			if (x3 > 10. || x3 < -10.)
			{
				cout << "problem" << endl;
			}
			return x3;
		}
	}
	else
	{
		//find the only real root
		A = -pow(R + sqrt(R2-Q3),1./3.);
		if (A == 0) B = 0;
		else B = Q/A;
		x1 = (A + B) - aover3;
		if (x1 > 10. || x1 < -10.)
		{
			cout << "problem" << endl;
		}
		return x1;
	}
}



#endif

