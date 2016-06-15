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
#include "coastline_v0.12.hpp"
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
	vector<double> Empty(NoNodes, -9999);
	EmptyVector = Empty;
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
	CalculateMeanNodeSpacing();
	DesiredNodeSpacing = MeanNodeSpacing;

}

void Coastline::Initialise(string xyfilename, int StartTime)
{

	/*	Read coastline.XY text file */

	cout << "\nCoastline.Initialise: Initialising Coastline from XY file: " << xyfilename << " at Time: " << StartTime << endl;

	ReadCoast(xyfilename, StartTime);

	//Populate empty vectors
	vector<double> Empty(NoNodes, -9999);
	EmptyVector = Empty;
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
	double Distance = 0;
	MeanNodeSpacing = (double)NodeSpacing;
	DesiredNodeSpacing = MeanNodeSpacing;
	
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
			tempX = tempX + MeanNodeSpacing*sin(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
			tempY = tempY + MeanNodeSpacing*cos(Trend*M_PI/180) + 2*(-0.5+((double)rand()/RAND_MAX));
			X.push_back(tempX);
			Y.push_back(tempY);
			NoNodes += 1;
		}
		Distance += MeanNodeSpacing;
	}
	
	//Populate empty vectors
	vector<double> Empty(NoNodes, -9999);
	EmptyVector = Empty;
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
	cout << "Coastline: Time is " << setprecision(3) << fixed << Time << " years\r";

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
	if (NodesToRemove.size() > 0)
	{
		//loop through flags and delete required nodes
		for (unsigned int j=0; j<NodesToRemove.size(); ++j)
		{
			X.erase(X.begin()+j);
			Y.erase(Y.begin()+j);
			EmptyVector.erase(EmptyVector.end()-1);
			NoNodes -= 1;
		}
		//empty the vector with delete flags
		NodesToRemove.erase(NodesToRemove.begin(),NodesToRemove.end());
	}
		
	//Reset all vectors (could builda function for this)
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
		else if (i>NoNodes-3)
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
			else cout << "es ist unmöglich" << endl;

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
			if (i==0) {}
			else if (i==NoNodes-1) {}
			else 
			{	
				Distance = sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));
				
				//if distance between two cells is too small, 
				//replace with a linearly interpolated point between
				if (Distance < 0.66*DesiredNodeSpacing)
				{
					//if cell is too narrow remove its adjacent nodes
					X[i] = ((X[i+1]+X[i])/2.);
					Y[i] = ((Y[i+1]+Y[i])/2.);
					X.erase(X.begin()+i+1); 
					Y.erase(Y.begin()+i+1); 
					EmptyVector.erase(EmptyVector.begin()+i+1);
					NoNodes -= 1;
					SpacingFlag = 0;
					break;
				}
				else if 
				(Distance > 1.5*DesiredNodeSpacing) 
				{
					//if cell too wide add two new nodes
					//just do with a straightline for now
					X.insert(X.begin()+i+1,(X[i]+X[i+1])/2.);
					Y.insert(Y.begin()+i+1,(Y[i]+Y[i+1])/2.);
					EmptyVector.insert(EmptyVector.begin()+i+1,-9999);
					NoNodes += 1;
					SpacingFlag = 0;
					break;
				}
				else
				{
					SpacingFlag = 1;
				}
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

	double C_0, L_0, Alpha_0, H_0, Theta_0;
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

	//handle special case just upcoast of shadow for negative (upcoast) flux
	//use "headland" orientation
	if(Shadows[i] == 4)
	{
		//determine alpha_0 angle between coast and wave crest approach
		if (Theta_0 <= Orientation[i+1]) Alpha_0 = FluxOrientation[i+1]-Theta_0-90.;
    	else if (Theta_0 > Orientation[i+1]+270.) Alpha_0 = (FluxOrientation[i+1]+270.)-Theta_0;
	   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i+1]);
	}
	else
	{
		//determine alpha_0 angle between coast and wave crest approach
		if (Theta_0 <= FluxOrientation[i]) Alpha_0 = FluxOrientation[i]-Theta_0-90.;
    	else if (Theta_0 > FluxOrientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-Theta_0;
	    else Alpha_0= 270.-(Theta_0-FluxOrientation[i]);
	    //if high angle waves then use flux orientation of the previous node (updrift)
	    if (Alpha_0 > 40.0)
	    {
	    	if (Theta_0 <= Orientation[i]) Alpha_0 = FluxOrientation[i-1]-Theta_0-90.;
	    	else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (FluxOrientation[i-1]+270.)-Theta_0;
		   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i-1]);
	    }
		else if (Alpha_0 < -40.0)
		{
			if (Theta_0 <= Orientation[i]) Alpha_0 = FluxOrientation[i+1]-Theta_0-90.;
	    	else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (FluxOrientation[i+1]+270.)-Theta_0;
		   	else Alpha_0= 270.-(Theta_0-FluxOrientation[i+1]);
		}
	}

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

	These Codes are required for handling sediment transport and
	for the refraction/diffraction code
	
	*/
	
	//declare temporary variables
	double ShadowAngle, X2, Y2, ShadowTest, TempOrientation;
	int i, I;
	
	//reset shadows
	vector<double> Zeros(NoNodes, 0);
	Shadows = Zeros;
	ShadowZoneWaveDirection = Zeros;
	ShadowZoneWaveHeight = Zeros;
	
	//get incloming wave direction and project to shadow angle
	ShadowAngle = OffshoreWaveDirection+180.;
	if (ShadowAngle >= 360.) ShadowAngle -= 360.;

	ShadowFlag = 0;
	
	//loop through
	for (i=1; i<NoNodes-1; ++i)
	{
		//for backward loop simultaneously
		//I = NoNodes-i-1;
		
		//test if already in shadow for forward along coast
		if (Shadows[i] < 1)
		{
			//create a hypothetical point for shadow vector
			X2 = X[i] + 10.*sin((M_PI/180.)*ShadowAngle);
            Y2 = Y[i] + 10.*cos((M_PI/180.)*ShadowAngle);
            
            //test if next point along coast is in shadow
            ShadowTest = ((X2 - X[i])*(Y[i+1] - Y[i]) - (Y2 - Y[i])*(X[i+1] - X[i]));
            if ((ShadowTest <= 0) && (abs(ShadowAngle-FluxOrientation[i]) < 90.))
            {
            	//label the shadow-casting node
            	Shadows[i] = 2;
            	//label the shadowed node
            	Shadows[i+1] = 1;
            	ShadowFlag = 1;
            	
            	//Continue through nodes to test if still in shadow
            	for (int j=i+2; j<NoNodes-1; ++j)
            	{
                    ShadowTest = ((X2 - X[i])*(Y[j] - Y[i]) - (Y2 - Y[i])*(X[j] - X[i]));
                    //Test if still in Shadow
                    if (ShadowTest <= 0) Shadows[j] = 1;
                    else break;
                }
            }
        }
     }
     for (I=NoNodes-2; I>1; --I)
     {
        if (Shadows[I] < 1)
        {
        	//create a hypothetical point for shadow vector
			X2 = X[I] + 10.*sin((M_PI/180.)*ShadowAngle);
            Y2 = Y[I] + 10.*cos((M_PI/180.)*ShadowAngle);
            
            //test if next point along coast is in shadow
            ShadowTest = ((X2 - X[I])*(Y[I-1] - Y[I]) - (Y2 - Y[I])*(X[I-1] - X[I]));
            
            //Invert FluxOrientaion
            TempOrientation = FluxOrientation[I-1]+180.;
            if (TempOrientation >= 360.) TempOrientation = TempOrientation-360.;
                        
            //if ((ShadowTest >= 0) && (abs(ShadowAngle - TempOrientation) < 45.))
            if ((ShadowTest >= 0 && (abs(ShadowAngle-TempOrientation) < 90.)))
            {
               //First cell should not be in shadow when working backward through vector
               //since flux is always calculated to/from the right boundary
               //Give it a special value (4)
               //Label shadow-casting node (3)
               Shadows[I] = 3;            
               Shadows[I-1] = 4;
               ShadowFlag = 1;
               ;
               //Continue through nodes to test if still in shadow
            	for (int J=I-2; J>0; --J)             
                {
                	ShadowTest = ((X2 - X[I])*(Y[J] - Y[I]) - (Y2 - Y[I])*(X[J] - X[I]));
                    //Test if still in Shadow
                    if (ShadowTest >= 0) Shadows[J] = 1;
                    else 
                    {
                    	Shadows[J] = 1;
                    	break;
                    }
				}
			}
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
	double ShadowAngle, SpatialGradient, TempWaveAngle, TempWaveHeight, dX, dY, AngleDiff; //D1, D2, 
	int i, j;
	
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
			while (Shadows[j] != 0) ++j;
			
			//Get breaking wave height and angle at shadow zone edge
			TransformWave(j);
			TempWaveAngle = FluxOrientation[j]-90.-BreakingWaveAngle[j];
			TempWaveHeight = BreakingWaveHeight[j];
			
			//get distance across shadow zone
			dX = X[j]-X[i];
			dY = Y[j]-Y[i];
			//D1 = sqrt(dX*dX + dY*dY);
				
			--j;
			
			//loop back through and set BreakingWaveHeight and Direction in the shadow
			while ((j > i) && (j < NoNodes-1))
			{
				if (j < 1) cout << "PROBLEM!" << endl;
				if (j > NoNodes-2) cout << "PROBLEM!" << endl;
				if (i < 1) cout << "PROBLEM!" << endl;
				if (i > NoNodes-2) cout << "PROBLEM!" << endl;
				
				//calculate shadow zone angle
				dX = X[j]-X[i];
				dY = Y[j]-Y[i];
				SpatialGradient = dY/dX;
							
				//convert to azimuths and determine refracted diffracted wave properties
				if (dX == 0 && dY < 0) 
				{
					ShadowZoneWaveDirection[j] = (OffshoreWaveDirection + 1.5*(180.-TempWaveAngle))+90.;
					ShadowZoneWaveHeight[j] = TempWaveHeight*cos((M_PI/180.)*(OffshoreWaveDirection-180.));
				}
				else if (dX == 0 && dY > 0) 
				{
					ShadowZoneWaveDirection[j] = (TempWaveAngle + 1.5*(360.-TempWaveAngle))+90.;
					ShadowZoneWaveHeight[j] = TempWaveHeight*cos((M_PI/180.)*(OffshoreWaveDirection-360.));
				}
				else if (dX > 0)
				{
					AngleDiff = (M_PI*0.5-atan(SpatialGradient))-((M_PI/180.)*(ShadowAngle));
					ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*90.*(180./M_PI)*sin(AngleDiff)*sin(AngleDiff);
					ShadowZoneWaveHeight[j] = OffshoreWaveHeight*cos(AngleDiff)*cos(AngleDiff)*cos(AngleDiff);
				}
				else if (dX < 0) 
				{
					AngleDiff = (M_PI*1.5-atan(SpatialGradient))-((M_PI/180.)*(ShadowAngle));
					ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*(180./M_PI)*AngleDiff;
					ShadowZoneWaveHeight[j] = OffshoreWaveHeight*cos(AngleDiff)*cos(AngleDiff)*cos(AngleDiff);
					}
				else 
				{
					cout << "This is impossible!" << endl;
				}
				// some checks
				if (ShadowZoneWaveHeight[j] < 0.) ShadowZoneWaveHeight[j] = 0.;
				if (ShadowZoneWaveDirection[j] == -9999 || ShadowZoneWaveHeight[j] == -9999) cout << "Diffraction FAILS!" << endl;
				--j;
			}
		}
		else if (Shadows[i] == 3)
		{
			//find end of shadow zone
			j=i-1;
			while (Shadows[j] != 0) --j;
			
			//flux is on the right of each cell so go back 1 more nodes
						
			//Get breaking wave height and angle at shadow zone edge
			TransformWave(j);
			TempWaveAngle = BreakingWaveAngle[j];
			TempWaveHeight = BreakingWaveHeight[j];
			
			//loop back through and set BreakingWaveHeight and Direction in the shadow
			while ((j < i-1) && (j > 1))
			{
				if (j < 1) cout << "PROBLEM!" << endl;
				if (j > NoNodes-2) cout << "PROBLEM!" << endl;
				if (i < 1) cout << "PROBLEM!" << endl;
				if (i > NoNodes-2) cout << "PROBLEM!" << endl;
				
				//calculate shadow zone angle to node SUPPLYING sediment
				dX = X[j+1]-X[i];
				dY = Y[j+1]-Y[i];
				SpatialGradient = dY/dX;
			
				//convert to azimuths and determine refracted diffracted wave properties
				if (dX == 0 && dY < 0) 
				{
					ShadowZoneWaveDirection[j] = (OffshoreWaveDirection + 1.5*(180.-TempWaveAngle))+90.;
					ShadowZoneWaveHeight[j] = TempWaveHeight*cos((M_PI/180.)*(OffshoreWaveDirection-180.));
				}
				else if (dX == 0 && dY > 0) 
				{
					ShadowZoneWaveDirection[j] = (TempWaveAngle + 1.5*(360.-TempWaveAngle))+90.;
					ShadowZoneWaveHeight[j] = TempWaveHeight*cos((M_PI/180.)*(OffshoreWaveDirection-360.));
				}
				else if (dX > 0)
				{
					AngleDiff = (M_PI*0.5-atan(SpatialGradient))-((M_PI/180.)*(ShadowAngle));
					ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*90.*(180./M_PI)*sin(AngleDiff)*sin(AngleDiff);
					ShadowZoneWaveHeight[j] = OffshoreWaveHeight*cos(AngleDiff)*cos(AngleDiff)*cos(AngleDiff);
				}
				else if (dX < 0) 
				{
					AngleDiff = (M_PI*1.5-atan(SpatialGradient))-((M_PI/180.)*(ShadowAngle));
					ShadowZoneWaveDirection[j] = OffshoreWaveDirection + 1.5*(180./M_PI)*AngleDiff;
					ShadowZoneWaveHeight[j] = OffshoreWaveHeight*cos(AngleDiff)*cos(AngleDiff)*cos(AngleDiff);
				}
				else 
				{
					cout << "This is impossible!" << endl;
				}
				// some checks
				if (ShadowZoneWaveHeight[j] < 0.) ShadowZoneWaveHeight[j] = 0.;
				if (ShadowZoneWaveDirection[j] == -9999 || ShadowZoneWaveHeight[j] == -9999) cout << "Diffraction FAILS!" << endl;
				++j;
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
					cout << "Intersection in coastline detected!!" << endl;
					cout << "i = " << i << "; j = " << j << endl;
				}
			}
		}
	}
}

void Coastline::TransportSediment(double &TimeDelta, double MaxTimeStep, Wave TheWave, int FluxType)
{
	/* DESCRIPTION GOES HERE
	
	Function to calculate longshore sediment flux and move sediment along the coast
	
	*/

	//Declare temporary variables
	double e1, e2, a, b, c;
	TimeDelta = 0;

	//get wave properties to individual variables
	OffshoreWavePeriod = TheWave.Get_WavePeriod();
	OffshoreWaveHeight = TheWave.Get_WaveHeight();
	OffshoreWaveDirection = TheWave.Get_WaveDirection();
	
	//Check Node Spacing is ok
	CheckNodeSpacing();

	//Find Shadows
	GetShadows();
	if (ShadowFlag == 1) 
	{
		RefractDiffractShadowZone();
	}
	
	//Transform waves
	TransformWaves();
	
	//Compute Sediment Fluxes
	for (int i=0; i<NoNodes; ++i)
	{
		//Calculate Longshore Flux
		CalculateFlux(i, FluxType);
		
		//Handle Boundary Conditions
		//PERIODIC BOUNDARY
		if (i==0 && StartBoundary == 1)
		{
			//Sediment Flux for last node
			CalculateFlux(NoNodes-2, FluxType);
										
			//Get Volume Change Vout - Vin - Vlost
			VolumeChange[i] = (LongshoreFlux[NoNodes-2]-LongshoreFlux[i])*(1.-LostFluxFraction);
		}
		else if (i==NoNodes-1 && EndBoundary == 1) VolumeChange[i] = VolumeChange[0];
		
		//FIXED BOUNDARY //only allow transport towards end node 
		else if (i <= 1 && StartBoundary == 2)
		{
			//nothing gets washed in, only washed out
			VolumeChange[i] = 0;
			if (LongshoreFlux[i] > 0) LongshoreFlux[i] = 0;
		}
		else if (i>=NoNodes-2 && EndBoundary == 2) 
		{
			VolumeChange[i] = 0;
			if (LongshoreFlux[i] < 0) LongshoreFlux[i] = 0;
		}
		else if (i==NoNodes-3 && EndBoundary == 2)
		{
			if (LongshoreFlux[i] < 0) LongshoreFlux[i] = 0;
			VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction);
		}
		//all other cells
		else VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction);		
	}
	
	//see previous versions for dynamic timestep code if need to put it back in
	TimeDelta = MaxTimeStep;

	//Move the coast
	for (int i=0; i<NoNodes; ++i)
	{
		//Convert to Shoreline position change
		VolumeChange[i] = VolumeChange[i]*TimeDelta;
		
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
						
			//This means cell destroys itself!?
			PositionChange[i] = 0;
			//Flag to destroy
			NodesToRemove.push_back(i);
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
