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
#include "coastline_v0.25.hpp"
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
	//vector<double> Zeros(NoNodes, 0);
	//Fixed = Zeros;
	
	ClosureDepth = 10.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	ShorefaceSlope = 0.01;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	
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
	//vector<double> Zeros(NoNodes, 0);
	//Fixed = Zeros;
	
	ClosureDepth = 10.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	ShorefaceSlope = 0.01;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
		
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
	//Trend = Trend;

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

	X.reserve(128);
	Y.reserve(128);
	Distance.reserve(128);
	Orientation.reserve(128);
	FluxOrientation.reserve(128);
	RegionalOrientation.reserve(128);
	Curvature.reserve(128);
	BreakingWaveHeight.reserve(128);
	BreakingWaveAngle.reserve(128);
	CellWidth.reserve(128);
	LongshoreFlux.reserve(128);
	VolumeChange.reserve(128);
	PositionChange.reserve(128);
	Shadows.reserve(128);
	ShadowZoneWaveDirection.reserve(128);
	ShadowZoneWaveHeight.reserve(128);
	
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
	
	ClosureDepth = 10.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.0;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	ShorefaceSlope = 0.01;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	
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

		//

		//Write Boundary Conditions to screen
		if 		(StartBoundary == 1 && EndBoundary == 1) cout << "Coastline.ReadCoast: Boundaries are periodic" << endl;
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

				while (iss >> TempY)
				{
					Y.push_back(TempY);
				}

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
	if (oftest)
	{
		FileExists = 1;
	}
	oftest.close();

	//open the output filestream and write headers
	ofstream WriteCoastFile;
	if (FileExists == 0)
	{
		WriteCoastFile.open(OutputFileName.c_str());
		if (WriteCoastFile.is_open())
		{
			WriteCoastFile << StartBoundary << " " << EndBoundary << endl;
		}
	}
	WriteCoastFile.close();

	//open output filestream again to append coastline data
	WriteCoastFile.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteCoastFile.is_open())
	{
		//write X
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i<NoNodes; ++i)
		{
			WriteCoastFile << setprecision(10) << " " << X[i];
		}
		WriteCoastFile << endl;

		//write Y
		WriteCoastFile << setprecision(4) << Time;
		for (int i=0; i< NoNodes; ++i)
		{
			WriteCoastFile << setprecision(10) << " " << Y[i];
		}
		WriteCoastFile << endl;
	}
	
	else
	{
		//report errors
		cout << "Coastline.WriteCoast: Error, the file " << OutputFileName << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
	
	//cout << "File Written" << endl;
}

//void Coastline::GetArea()
//{
	//or incorporate into get morphology
//}

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
	double SpatialGradient, SpatialGradientLast;
	
	//check for nodes to remove
//	if (NodesToRemove.size() > 0)
//	{
//		//loop through flags and delete required nodes
//		for (int j=0; j<(int)NodesToRemove.size(); ++j)
//		{
//			if (j>=NoNodes-2 || j <= 1)
//			{
//				cout << "Removing Fixed Node!" << j << "/" << NoNodes << endl;
//			}
//			X.erase(X.begin()+j);
//			Y.erase(Y.begin()+j);
//			EmptyVector.erase(EmptyVector.end()-1);
//			
//			NoNodes -= 1;
//		}
//		//empty the vector with delete flags
//		NodesToRemove.erase(NodesToRemove.begin(),NodesToRemove.end());
//	}
		
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
			else
			{
				cout << "This is impossible!" << endl;
			}

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
			//CellWidth[i] = cos(M_PI/180.*(Trend-Orientation[i]))*sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));
			CellWidth[i] = (sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i]));
		}

		//Deal with EndBoundary
		else if (i==NoNodes-1)
		{
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
			
			//Cell Width parallel to orientation AT the node position
			CellWidth[i] = (0.5*sqrt(pow(X[i]-X[i-1],2.0) + pow(Y[i]-Y[i-1],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i-1]));
			CellWidth[i] += (0.5*sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i]));
		}
	}
}

void Coastline::CheckNodeSpacing()
{
	//declare variables
	int SpacingFlag = 0;
	double Distance;
	//need to recalculate morphology adn start check again after each new node is added
	while (SpacingFlag == 0)
	{
		for (int i=0; i<NoNodes-1; ++i)
		{
			if (i<1) {}
			else if (i>NoNodes-3) {}
			else 
			{	
				Distance = sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));
				
				//if distance between two cells is too small, 
				//replace with a linearly interpolated point between
				if (Distance < 0.66*DesiredNodeSpacing)
				{
					//if cell is too narrow remove its adjacent nodes
					if (i == 1 || i == NoNodes-3) 
					{
//						X.erase(X.begin()+i+1); 
//						Y.erase(Y.begin()+i+1); 
//						EmptyVector.erase(EmptyVector.begin()+i+1);
//						NoNodes -= 1;
//						SpacingFlag = 0;
					}
					else
					{
						X[i] = ((X[i+1]+X[i])/2.);
						Y[i] = ((Y[i+1]+Y[i])/2.);
						X.erase(X.begin()+i+1); 
						Y.erase(Y.begin()+i+1); 
						//Fixed.erase(Fixed.begin()+i+1);
						NoNodes -= 1;
						SpacingFlag = 0;
					}
					break;
				}
				else if 
				(Distance > 1.5*DesiredNodeSpacing) 
				{
					//if cell too far apart add new node
					//just do with a straightline for now
					X.insert(X.begin()+i+1,(X[i]+X[i+1])/2.);
					Y.insert(Y.begin()+i+1,(Y[i]+Y[i+1])/2.);
					//Fixed.insert(Fixed.begin()+i+1,0.0);
					NoNodes += 1;
					SpacingFlag = 0;
					break;
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

//	//handle special case just upcoast of shadow for negative (upcoast) flux
//	//use "headland" orientation
//	if(Shadows[i] == 4)
//	{
//		//determine alpha_0 angle between coast and wave crest approach
//		if (Theta_0 <= Orientation[i]) Alpha_0 = FluxOrientation[i]-Theta_0-90.;
//    	else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-Theta_0;
//	   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i]);
//	}
//	else if (Shadows[i] == 2)
//	{
//		//determine alpha_0 angle between coast and wave crest approach for this and next cell!
//		//cout << "Haven't cheked this" << endl;
//		if (Theta_0 <= Orientation[i]) Alpha_0 = FluxOrientation[i]-Theta_0-90.;
//    	else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-Theta_0;
//	   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i]);
//	}
//	
//	else
//	{
		//determine alpha_0 angle between coast and wave crest approach
		if (Shadows[i] != 0)
	   	{
	   		if (ShadowZoneWaveDirection[i] <= FluxOrientation[i]) Alpha_0 = FluxOrientation[i]-ShadowZoneWaveDirection[i]-90.;
    		else if (ShadowZoneWaveDirection[i] > FluxOrientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-ShadowZoneWaveDirection[i];
	    	else Alpha_0= 270.-(ShadowZoneWaveDirection[i]-FluxOrientation[i]);
	   	}
	   	else
	   	{
	   		if (Theta_0 <= FluxOrientation[i]) Alpha_0 = FluxOrientation[i]-Theta_0-90.;
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
	   	if (Alpha_0 <= -90. || Alpha_0 >= 90.) Alpha_0 = 0;
	   	//next if alpha_0 and alpha_0_last have the same sign
	   	else if ((Alpha_0_Last > 0) && (Alpha_0 > 0))
	   	{
	    	if ((Alpha_0_Last < 45.) && (Alpha_0 > 45.)) Alpha_0 = 45.;
	    	else if ((Alpha_0_Last > 45.) && (Alpha_0 < 45.)) Alpha_0 = 45.;
	    }
	    //next if alpha_0 and alpha_0_next have the same sign
	    else if ((Alpha_0_Next < 0) && (Alpha_0 < 0))
	    {
	    	if ((Alpha_0_Next > -45.) && (Alpha_0 < -45.)) Alpha_0 = -45.;
	    	else if ((Alpha_0_Next < -45.) && (Alpha_0 > -45.)) Alpha_0 = -45.;
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
	for (i=1; i<NoNodes-2; ++i)
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
        
        	for (int j=0; j<NoNodes-2; ++j)
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

	        //set point for casting shadow    
	        X1 = X[ShadowStart];
    	    Y1 = Y[ShadowStart];
        
        	//create a hypothetical point for vector in direction of ShadowAngle
	        X2 = X[ShadowStart] + 100000.*sin((M_PI/180.)*ShadowAngle);
    	    Y2 = Y[ShadowStart] + 100000.*cos((M_PI/180.)*ShadowAngle);
        	dX12 = X2-X1;
	        dY12 = Y2-Y1;
        
        	CasterFlag = 0;
        
        	for (int j=0; j<NoNodes-2; ++j)
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
    	            else CasterFlag = 1; break;
				}
        	}
                
        	//if shadowstart casts the shadow set Shadows to 2
        	//otherwise shadowend must be the caster
        	if (CasterFlag == 1) Shadows[ShadowStart] = 2;
        	else
        	{
        	    Shadows[ShadowEnd] = 3;
        	    Shadows[ShadowEnd-1] = 4;
        	    Shadows[ShadowEnd-2] = 1;
        	}    
        
        	//continue from after shadowed section
        	i = ShadowEnd;
        }
	}
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
	int i, j, k;
	
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
			j=i+1;
			while ((Shadows[j] != 0) && (j<NoNodes-1)) ++j;
			k = j;
			--j;

			if (Shadows[j] == 3) Shadows[j] = 1;
			if (Shadows[k] == 3) Shadows[k] = 1;
			
			//get distance the length of the shore across the shadow zone
			//get gradients of intersecting lines
			M1 = 1./tan((M_PI/180.)*ShadowAngle);
			M2 = 1./tan((M_PI/180.)*FluxOrientation[j]);
			//find point of intersection
			X0 = (Y[j] - M2*X[j] + M1*X[i] - Y[i])/(M1 - M2);
			Y0 = M1*X0 + Y[i] - M1*X[i];
			//Find distance difference
			dX = X0-X[j];
			dY = Y0-Y[j];
			Dist = sqrt(dX*dX + dY*dY);
			ShadowZoneDistance = Distance[j]-Distance[i]+Dist;
			
			//loop back through and set BreakingWaveHeight and Direction in the shadow
			while (j > i)
			{
				//calculate shadow zone angle
				dX = X[j]-X[i];
				dY = Y[j]-Y[i];
				SpatialGradient = dY/dX;
							
				//convert to azimuths and determine refracted diffracted wave properties
				if (dX == 0 && dY < 0) AngleDiff = ShadowAngle-180.;
				else if (dX == 0 && dY > 0) AngleDiff = ShadowAngle-360;
				else if (dX > 0) AngleDiff = (180./M_PI)*(M_PI*0.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				else if (dX < 0) AngleDiff = (180./M_PI)*(M_PI*1.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*AngleDiff;
				if (1.5*AngleDiff > 90.) ShadowZoneWaveHeight[j] = 0;
				else ShadowZoneWaveHeight[j] = OffshoreWaveHeight*0.5*(1-sin((M_PI/180.)*(AngleDiff)));
				--j;
			}
			
			//Reduce wave heights outside the shadow zone over the length that they were redued in the shadow zone
			//This should preserve energy.
			
			if (Shadows[i] == 1) cout << "What the fuck!?" << endl << endl << endl;
			
			Dist = (Distance[k]-Distance[k-1])-Dist;
			while(Dist <= ShadowZoneDistance)
			{
				ShadowZoneWaveHeight[k] = OffshoreWaveHeight*(1.-0.5*(1.-sin((M_PI/180.)*(90.*(Dist/ShadowZoneDistance)))));
				ShadowZoneWaveDirection[k] = OffshoreWaveDirection;
				Shadows[k] = 1;
				k+=1;
				Dist += (Distance[k]-Distance[k-1]);
				if (k == NoNodes-1 || k == 0) 
				{
					//cout << endl << "got to end of array" << endl;
					Dist = ShadowZoneDistance+1;
				}
			}
		}
		
		else if (Shadows[i] == 3)
		{
			//find end of shadow zone
			j=i-1;
			while ((Shadows[j] != 0) && (j>0)) --j;
			k = j-1;
			
			//get distance the length of the shore across the shadow zone
			//get gradients of intersecting lines1
			M1 = 1./tan((M_PI/180.)*ShadowAngle);
			M2 = 1./tan((M_PI/180.)*FluxOrientation[j]);
			//find point of intersection
			X0 = (Y[j] - M2*X[j] + M1*X[i] - Y[i])/(M1 - M2);
			Y0 = M1*X0 + Y[i] - M1*X[i];
			//Find distance difference
			dX = X0-X[j];
			dY = Y0-Y[j];
			Dist = sqrt(dX*dX + dY*dY);
			ShadowZoneDistance = Distance[i]-Distance[j]-Dist;
			
			//loop back through and set BreakingWaveHeight and Direction in the shadow
			while (j < i-1)
			{
				//calculate shadow zone angle to node SUPPLYING sediment
				dX = X[j+1]-X[i];
				dY = Y[j+1]-Y[i];
				SpatialGradient = dY/dX;
							
				//convert to azimuths and determine refracted diffracted wave properties
				if (dX == 0 && dY < 0) AngleDiff = ShadowAngle-180.;
				else if (dX == 0 && dY > 0) AngleDiff = ShadowAngle-360;
				else if (dX > 0) AngleDiff = (180./M_PI)*(M_PI*0.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				else if (dX < 0) AngleDiff = (180./M_PI)*(M_PI*1.5-atan(SpatialGradient)-(M_PI/180.)*(ShadowAngle));
				if (1.5*AngleDiff < -90.) TempWaveHeight = 0;
				else TempWaveHeight = OffshoreWaveHeight*0.5*(1-sin((M_PI/180.)*(-AngleDiff)));
				if (TempWaveHeight < ShadowZoneWaveHeight[j] || ShadowZoneWaveHeight[j] == -9999)
				{
					ShadowZoneWaveHeight[j] = TempWaveHeight;
					ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*AngleDiff;
				}
				++j;
			}
			
			//Reduce wave heights outside the shadow zone over the length that they were redued in the shadow zone
			//This should preserve energy.
			//Dist = (Distance[k+1]-Distance[k]-Dist);
			while(Dist < ShadowZoneDistance)
			{
				ShadowZoneWaveHeight[k] = OffshoreWaveHeight*(1.-0.5*(1.-sin((M_PI/180.)*(90.*(Dist/ShadowZoneDistance)))));
				ShadowZoneWaveDirection[k] = OffshoreWaveDirection;
				if (Shadows[k+1] != 1) 
				{
					Shadows[k+1] = 1;
					k-=1;
					Dist += (Distance[k+1]-Distance[k]);
				}
				else Dist = ShadowZoneDistance+1;
				if (k<2 || k >NoNodes-3) 
				{
					cout << "Fell off the front!" << endl;
					Dist = ShadowZoneDistance+1;
				}
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
					//if ((Fixed[i] == 0) && (Fixed[j] == 0))
					//{	
						cout << "Intersection in coastline detected!!" << endl;
						cout << "i = " << i << "; j = " << j << "; NoNodes = " << NoNodes << endl;
					//}

					//if (i==0) Fixed[j] = 1;
					//else if (j==NoNodes-2) Fixed[i] = 1;
					//else
					//{
						X.erase(X.begin()+i+1, X.begin()+j+1);
						Y.erase(Y.begin()+i+1, Y.begin()+j+1);
						//Fixed.erase(Fixed.begin()+i+1, Fixed.begin()+j+1);
						NoNodes -= (j-i);
					//}
				}
			}
		}
	}
	CheckNodeSpacing();
	CalculateMorphology();
}

void Coastline::TransportSediment(double &TimeDelta, double MaxTimeStep, Wave TheWave, int FluxType, int RefDiffFlag, double FluxFraction)
{
	/* DESCRIPTION GOES HERE
	
	Function to calculate longshore sediment flux and move sediment along the coast
	
	*/

	//Declare temporary variables
	double e1, e2, a, b, c;
	//TimeDelta = 0;

	//get wave properties to individual variables
	OffshoreWavePeriod = TheWave.Get_WavePeriod();
	OffshoreWaveHeight = TheWave.Get_WaveHeight();
	OffshoreWaveDirection = TheWave.Get_WaveDirection();
	
	//handle offshore waves
	if ((OffshoreWaveDirection < (FluxOrientation[0]-170.)) || (OffshoreWaveDirection > (FluxOrientation[NoNodes-1]-10.))) return;

	//Check Node Spacing is ok
	CheckNodeSpacing();

	//Find Shadows
	GetShadows();
	if (ShadowFlag == 1 && RefDiffFlag == 1) RefractDiffractShadowZone();
	
	//Transform waves
	TransformWaves();
	
	//Compute Sediment Fluxes and move coast
	for (int i=0; i<NoNodes; ++i)
	{
		//Calculate Longshore Flux
		CalculateFlux(i, FluxType);
		
		//Calculate volume change
		
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
		//else if (Fixed[i] == 1)
		//{
		//	LongshoreFlux[i] = 0;
		//	VolumeChange[i] = 0;
		//}
		else VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction)*TimeDelta;		

		//move the coastline
				
		//new position change for trapezium-shaped cells
		//a, b, c are parameters for quadratic equation solver
		if (i==0) e1 = Orientation[i]-FluxOrientation[i];
		else e1 = Orientation[i]-FluxOrientation[i-1];
		e2 = FluxOrientation[i]-Orientation[i];
		
		//one of these coefficients needs some sort of closure depth dependency?
		a = 0.5*(tan((M_PI/180.)*e1) + tan((M_PI/180.)*e2));
		b = CellWidth[i];
		c = -VolumeChange[i]/ClosureDepth;

		//solve quadratic equation to get position change
		//if cell is rectilinear so dont need to solve the quadratic
		if (VolumeChange[i] != 0)
		{
			if (a==0) PositionChange[i] = VolumeChange[i]/(CellWidth[i]*ClosureDepth);
			//quadratic solver
			else PositionChange[i] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
		}
		else PositionChange[i] = 0;

		if (PositionChange[i] != PositionChange[i])
		{
			//Check for NaN in PositionChange
			cout << "I am not equal to myself!!" << endl;
			cout << "Found NaN in Position Change: Node is " << i << endl;
					
			//Try using 3cell orientation to get flux @i-1 and i?
			if ((i > 2) && (i <NoNodes-2))
			{
				for (int j=i-1; j<(i+1); ++j)
				{
					e1 = Orientation[j]-FluxOrientation[j-1];
					e2 = FluxOrientation[j]-Orientation[j];
					FluxOrientation[j] = Orientation[j];
					TransformWave(j);
					CalculateFlux(j, FluxType);
					VolumeChange[j] = (LongshoreFlux[j-1]-LongshoreFlux[j])*(1.-LostFluxFraction)*TimeDelta;
					c = -VolumeChange[j]/ClosureDepth;
					if (VolumeChange[j] != 0)
					{
						if (a==0) PositionChange[j] = VolumeChange[j]/(CellWidth[j]*ClosureDepth);
						//quadratic solver
						else PositionChange[j] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
					}
					if (PositionChange[j] != PositionChange[j])
					{
						cout << "I am not equal to myself!!" << endl;
						cout << "Found NaN in Position Change: Node is " << j << endl;
						VolumeChange[j] = 0;
						LongshoreFlux[j] = LongshoreFlux[j-1];
						PositionChange[j] = 0;
					}
				}
			}
			
			//This means cell destroys itself!?
			//VolumeChange[i] = 0;
			//LongshoreFlux[i] = LongshoreFlux[i-1];
			//PositionChange[i] = 0;
			
			//Flag to destroy
//			if ((i > 1) && (i <NoNodes-2)) 
//			{
//				NodesToRemove.push_back(i);
//				cout << "Deleting Node" << i << "/" << NoNodes << endl;
//			}
		}

		//effective angle for trend //more logic needed here
		// htis isnt needed for oreientation -| migration
		//Angle = 180.-Trend;

		//Decompose Position change into X and Y components
		//Going to need a bit more logic here to handle different orientations?
		//current setup x change is negative, y change is positive
		X[i] -= PositionChange[i]*cos((M_PI/180.)*Orientation[i]);
		Y[i] += PositionChange[i]*sin((M_PI/180.)*Orientation[i]);
				
	}
	if ((LongshoreFlux[1] != 0) && (Shadows[1] == 2))
	{
		cout << "FUCK!" << endl;
		TransformWave(NoNodes-4);
	}

	//FOR DEBUGGING, WRITE SOME STUFF TO A FILE
	//test if output file already exists
	string OutputFileName = "outputdebug.txt";
	int FileExists = 0;
	ifstream oftest(OutputFileName.c_str());
	if (oftest) FileExists = 1;
	oftest.close();

	//open the output filestream and write headers
	ofstream WriteDebugOutput;
	if (FileExists == 0)
	{
		WriteDebugOutput.open(OutputFileName.c_str());
		if (WriteDebugOutput.is_open())
		{
			WriteDebugOutput << "Header line" << endl;
		}
	}
	WriteDebugOutput.close();

	//open output filestream again to append coastline data
	WriteDebugOutput.open(OutputFileName.c_str(), fstream::app|fstream::out);

	//Check if file exists if not open a new one and write headers
	if (WriteDebugOutput.is_open())
	{
		//write X
		WriteDebugOutput << setprecision(4) << OffshoreWaveDirection << " ";
		for (int i=0; i<NoNodes; ++i) WriteDebugOutput << setprecision(10) << " " << Shadows[i];
		WriteDebugOutput << endl;

		//write Y
		WriteDebugOutput << setprecision(4) << OffshoreWaveHeight;
		for (int i=0; i< NoNodes; ++i) WriteDebugOutput << setprecision(10) << " " << LongshoreFlux[i];
		WriteDebugOutput << endl;
	}
	else
	{
		//report errors
		cout << "Error, the file " << WriteDebugOutput << " is not open or cannot be read." << endl;
		exit(EXIT_FAILURE);
	}
	
	//Check for intersection
	IntersectionAnalysis();
	
	//Update Coastal Morphology
	CalculateMorphology();
}

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
	CalculateMorphology();
}
#endif

