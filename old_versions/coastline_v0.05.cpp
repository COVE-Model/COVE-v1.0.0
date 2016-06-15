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
// 27/11/13 Non-rectilinear cell shapes added and node migration perpendicular
//				to the orientation rather than the trend
//
// 22/11/13 CERC diffusion added, rectilinear style
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
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "coastline_v0.05.hpp"
#include "waveclimate_v0.01.hpp"
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
	Curvature = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	CellWidth = EmptyVector;
	LongshoreFlux = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	Shadows = EmptyVector;
	ClosureDepth = 10.0;			//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.1;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE

	//Get general orientation
	//double dX = (X[NoNodes-1]-X[0]);
	//double dY = (Y[NoNodes-1]-Y[0]);
	//double TempSlope = dY/dX;

	//if (dX == 0 && dY < 0) Trend = 180;
	//else if (dX == 0 && dY > 0) Trend = 0;
	//else if (dX > 0) Trend = (180.0/M_PI)*((0.5*M_PI)-atan(TempSlope));
	//else if (dX < 0) Trend = (180.0/M_PI)*((1.5*M_PI)-atan(TempSlope));
	//else cout << "Why do you hate me?" << endl;

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
	Curvature = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	CellWidth = EmptyVector;
	LongshoreFlux = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	Shadows = EmptyVector;
	ClosureDepth = 10.0;			//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.1;		//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE

	//Get general orientation
	//double dX = (X[NoNodes-1]-X[0]);
	//double dY = (Y[NoNodes-1]-Y[0]);
	//double TempSlope = dY/dX;

	//if (dX == 0 && dY < 0) Trend = 180;
	//else if (dX == 0 && dY > 0) Trend = 0;
	//else if (dX > 0) Trend = (180.0/M_PI)*((0.5*M_PI)-atan(TempSlope));
	//else if (dX < 0) Trend = (180.0/M_PI)*((1.5*M_PI)-atan(TempSlope));
	//else cout << "Why do you hate me?" << endl;

	//Calculate Orientation, Slope and Curvature
	CalculateMorphology();
	CalculateMeanNodeSpacing();
	DesiredNodeSpacing = MeanNodeSpacing;

}

void Coastline::Initialise(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundaryInput, int EndBoundaryInput)
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
			<< "\n\t MeanNodeSpacing = " << MeanNodeSpacing
			<< "\n\t CoastLength = " << CoastLength
			//<< "\n\t Trend = " << Trend
			<< "\n\t Boundary Conditions = " << StartBoundary << ", " << EndBoundary << endl << endl;

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
	Curvature = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	CellWidth = EmptyVector;
	LongshoreFlux = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	Shadows = EmptyVector;
	
	ClosureDepth = 10.0;			//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	LostFluxFraction = 0.1;	//NEED TO FIND SOMEWHERE BETTER TO DEFINE THESE
	//Calculate Orientation, Slope and Curvature
	CalculateMorphology();

}

void Coastline::ReadCoast(string InputFileName, int Time)
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

void Coastline::WriteCoast(string OutputFileName, int Time)
{
	/*	Writes a coastline object X and Y coordinates to file for a given time
		If the file already exists the data will be appended else a new file is
		created.

		File format is 	Time | X[0] | X[1] | X[2] =====> X[NoNodes]
								Time | Y[0] | Y[1] | Y[2] =====> Y[NoNodes]			*/

	//Print to screen
	cout.flush();
	cout << "Coastline: Time is " << Time << " years\r";

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
		WriteCoastFile << Time;
		for (int i=0; i<NoNodes; ++i)
		{
			WriteCoastFile << " " << X[i];
		}
		WriteCoastFile << endl;

		//write Y
		WriteCoastFile << Time;
		for (int i=0; i< NoNodes; ++i)
		{
			WriteCoastFile << " " << Y[i];
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
	
	//Reset all vectors (could builda function for this)
	Orientation = EmptyVector;
	FluxOrientation = EmptyVector;
	Curvature = EmptyVector;
	CellWidth = EmptyVector;
	BreakingWaveHeight = EmptyVector;
	BreakingWaveAngle = EmptyVector;
	VolumeChange = EmptyVector;
	PositionChange = EmptyVector;
	LongshoreFlux = EmptyVector;
	Shadows = EmptyVector;
	
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
			else cout << "This is impossible!" << endl;

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
				Curvature[i] = Curvature[0];
				CellWidth[i] = CellWidth[0];
			}
			//FIXED BOUNDARY
			else if (EndBoundary == 2)
			{
				FluxOrientation[i] = FluxOrientation[i-1];
				Orientation[i] = Orientation[i-1];
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
			SpatialGradientLast = SpatialGradient;
			SpatialGradient = dY/dX;
			//convert to azimuths
			if (dX == 0 && dY < 0) FluxOrientation[i] = 180.;
			else if (dX == 0 && dY > 0) FluxOrientation[i] = 0.;
			else if (dX > 0) FluxOrientation[i] = (180./M_PI)*(M_PI*0.5 - atan(SpatialGradient));
			else if (dX < 0) FluxOrientation[i] = (180./M_PI)*(M_PI*1.5 - atan(SpatialGradient));
			else 
			{
				cout << "es imposible!" << endl;
			}

			//Next get curvature and orientation (3-node)
			dX = (X[i+1]-X[i-1]);
			dY = (Y[i+1]-Y[i-1]);
			Curvature[i] = (SpatialGradient-SpatialGradientLast)/sqrt(dX*dX + dY*dY);
			//convert to azimuths
			if (dX == 0 && dY < 0) Orientation[i] = 180.;
			else if (dX == 0 && dY > 0) Orientation[i] = 0.;
			else if (dX > 0) Orientation[i] = (180./M_PI)*(M_PI*0.5 - atan(dY/dX));
			else if (dX < 0) Orientation[i] = (180./M_PI)*(M_PI*1.5 - atan(dY/dX));
			else 
			{	
				cout << "es ist unmöglich" << endl;
			}

			//Get cell widths
			//comment out perpendicular to trend calcs
			//CellWidth[i] += 0.5*cos(M_PI/180.*(Trend-Orientation[i]))*sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));
			//CellWidth[i+1] = 0.5*cos(M_PI/180.*(Trend-Orientation[i]))*sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));

			//Cell Width parallel to orientation AT the node position
			CellWidth[i] = (0.5*sqrt(pow(X[i]-X[i-1],2.0) + pow(Y[i]-Y[i-1],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i-1]));
			CellWidth[i] += (0.5*sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0)))/cos(M_PI/180.*(Orientation[i]-FluxOrientation[i]));
		}
	}
	//CalculateMeanNodeSpacing(); //This may be needed once we allow cells to move in different directions
}

void Coastline::CheckNodeSpacing()
{
	//vector<int> SpacingFlag(NoNodes,0);
	int SpacingFlag = 0;
	double Distance;
	while (SpacingFlag == 0)
	{
		for (int i=0; i<NoNodes-1; ++i)
		{
			if (i==0) {}
			else if (i==NoNodes-1) {}
			else 
			{	
				Distance = sqrt(pow(X[i+1]-X[i],2.0) + pow(Y[i+1]-Y[i],2.0));
				if (CellWidth[i] < 0.66*DesiredNodeSpacing)
				{
					//if cell is too narrow remove its adjacent nodes
					X.erase(X.begin()+i+1); 
					X.erase(X.begin()+i-1);
					Y.erase(Y.begin()+i+1); 
					Y.erase(Y.begin()+i-1);
					EmptyVector.erase(EmptyVector.begin()+i-1);
					NoNodes -= 1;
					SpacingFlag = 0;
					//cout << "Removing Cells" << endl;
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
					Orientation.insert(Orientation.begin()+i+1,-9999);
					FluxOrientation.insert(FluxOrientation.begin()+i+1,-9999);
					Curvature.insert(Curvature.begin()+i+1,-9999);
					BreakingWaveHeight.insert(BreakingWaveHeight.begin()+i+1,-9999);
					BreakingWaveAngle.insert(BreakingWaveAngle.begin()+i+1,-9999);
					CellWidth.insert(CellWidth.begin()+i+1,-9999);
					VolumeChange.insert(VolumeChange.begin()+i+1,-9999);
					PositionChange.insert(PositionChange.begin()+i+1,-9999);
					LongshoreFlux.insert(LongshoreFlux.begin()+i+1,-9999);
					Shadows.insert(Shadows.begin()+i+1,-9999);
					
					NoNodes += 1;
					SpacingFlag = 0;
					//cout << "Adding Cells" << endl;
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
	double dx, dy, MeanDistanceDouble;
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
	MeanDistanceDouble = TotalDistance/NoNodes;
	//Convert to int
	MeanNodeSpacing = (int) MeanDistanceDouble;
}

void Coastline::SimpleDiffusion(double DiffCoeff, double TimeDelta)
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

void Coastline::TransformWaves(double T, double Theta_0, double H_0)
{
	//Transforms waves assume shore-parallel contours and fixed shoreface gradient
	//Wave base depth calculated based on wave conditions?

	// C_0 is offshore wave celerity, L_0 is offshore wave length,
	// C is wave speed, H is wave height, L is wavelength, alpha is wave angle to the coast
	// h is depth, n is a shoaling factor, k is wave number,
	//	Ks is shoaling coefficient, Kr is refraction coefficient

	//cout << "Inside TransformWaves" << endl;

	// See Komar (1998)

	double C_0, L_0, Alpha_0;
	double C, H, L, Alpha;
	double h, n, k, Ks, Kr;

	int BreakCondition = 0;

	C_0 = (g*T)/(2.0*M_PI);			//Deep water wave speed (m/s)
	L_0 = C_0*T;		//Deep water wave length (m)

	//loop through coastline and perform wave transformation
	for (int i=0; i<NoNodes; ++i)
	{
		//cout << i << " ";
		BreakCondition = 0;	//for testing wave breaking
		h = 3.*H_0;		// water depth at wave base (metres) calculate this later based on period and height

		//Catch offshore waves
		if (		((Theta_0 >= FluxOrientation[i]) && (Theta_0 <= FluxOrientation[i]+180.0))
				||	((Theta_0 >= FluxOrientation[i]) && (Theta_0 <= FluxOrientation[i]-180.0)))
		{
			BreakingWaveHeight[i] = 0;
			BreakingWaveAngle[i] = 0;
			BreakCondition = 1; // dont bother with refraction for this cell
		}
		//otherwise find incidence angle
		//Use FluxOrientation for low angle waves and Orientation for high angle waves
		else
		{
			if (Theta_0 <= FluxOrientation[i]) Alpha_0 = FluxOrientation[i]-Theta_0-90.;
    		else if (Theta_0 > FluxOrientation[i]+270.) Alpha_0 = (FluxOrientation[i]+270.)-Theta_0;
	    	else Alpha_0= 270.-(Theta_0-FluxOrientation[i]);
	    	if (Alpha_0 > 40.0 || Alpha_0 < -40.0)
	    	{
	    		if (Theta_0 <= Orientation[i]) Alpha_0 = Orientation[i]-Theta_0-90.;
	    		else if (Theta_0 > Orientation[i]+270.) Alpha_0 = (Orientation[i]+270.)-Theta_0;
		    	else Alpha_0= 270.-(Theta_0-Orientation[i]);
	    	}
		}

		while (BreakCondition == 0)
		{
			//Calculate new wave height
			L = L_0*sqrt(tanh((2.*M_PI*h)/L_0));		//Wavelength (m) in intermediate-shallow waters
			//L = L_0*pow(tanh(pow(pow(2.*M_PI/T,2.)*h/g,.75)),2.0/3.0); //Fenton & McKee (1990) formulation
         C = C_0*tanh((2*M_PI*h)/L);					//Wave speed (m/s) set by water depth, L and C_0
         k = 2*M_PI/L;										//Wave number (1/m)
         n = ((2*h*k)/(sinh(2*h*k)) + 1)/2;			//Shoaling factor

         //
         Ks = sqrt(C_0/(n*C*2));							//Shoaling coefficient
         Alpha = (180./M_PI)*asin((C/C_0)*sin((M_PI/180.)*Alpha_0));		//update theta
         Kr = sqrt(cos((M_PI/180.)*Alpha_0)/cos((M_PI/180.)*Alpha));		//refraction coefficient
         H = H_0*Ks*Kr;										//calculate new wave height

			//test if wave breaks
			if (H > h*0.78) //or something like that
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

void Coastline::GetShadows(double OffshoreWaveDirection)
{
	//declare temporary variables
	double ShadowAngle, dX, dY,SpatialGradient,TempOrientation;
	//reset shadows
	vector<double> Zeros(NoNodes, 0);
	Shadows = Zeros;

	//loop through but ignore start and end cells
	for (int i=1; i<NoNodes-1; ++i)
	{
		//Find equation for shadow from cell i based on BREAKING Wave Angle for now
		ShadowAngle = FluxOrientation[i-1]+OffshoreWaveDirection;
		//More Logic needed here

		//Loop through to find shadow
		for (int j=i+1; j<NoNodes-2; ++j)
		{
			//First get orientation to proposed shadow cell
			dX = X[j]-X[i];
			dY = Y[j]-Y[i];
			SpatialGradient = dY/dX;
			//convert to azimuths
			if (dX == 0 && dY < 0) TempOrientation = 180.;
			else if (dX == 0 && dY > 0) TempOrientation = 0.;
			else if (dX > 0) TempOrientation = (180./M_PI)*(M_PI*0.5 - atan(SpatialGradient));
			else if (dX < 0) TempOrientation = (180./M_PI)*(M_PI*1.5 - atan(SpatialGradient));

			//Check if shadowed
			//Need more logic here to handle the round the clock effect
			if (Shadows[j] == 1) break;
			else if (TempOrientation < ShadowAngle) break;
			else Shadows[j] = 1;
		}
		//Find equation for shadow from cell NoNodes-1-i based on BREAKING Wave Angle for now
		ShadowAngle = FluxOrientation[NoNodes-i-1]+OffshoreWaveDirection;
		//More Logic needed here

		//Loop through to find shadow
		for (int j=NoNodes-i-2; j>1; --j)
		{
			//First get orientation to proposed shadow cell
			dX = X[NoNodes-i-1]-X[j];
			dY = Y[NoNodes-i-1]-Y[j];
			SpatialGradient = dY/dX;
			//convert to azimuths
			if (dX == 0 && dY < 0) TempOrientation = 180.;
			else if (dX == 0 && dY > 0) TempOrientation = 0.;
			else if (dX > 0) TempOrientation = (180./M_PI)*(M_PI*0.5 - atan(SpatialGradient));
			else if (dX < 0) TempOrientation = (180./M_PI)*(M_PI*1.5 - atan(SpatialGradient));

			//Check if shadowed
			//Need more logic here to handle the round the clock effect
			if (Shadows[j] == 1) break;
			if (TempOrientation < ShadowAngle) break;
			else Shadows[j] = 1;
		}
	}
}

void Coastline::CERCDiffusion(double &TimeDelta, double MaxTimeStep, double OffshoreWavePeriod, double OffshoreWaveDirection, double OffshoreWaveHeight)
{
	/* DESCRIPTION GOES HERE
		CERC Equation Qls = 0.41*(Hb^2.5)*sin(alpha)*cos(alpha)
		predicts transport rate in m3/sec
	*/

	//cout << "inside CERCDiffusion" << endl;

	//Declare temporary variables
	double dt, e1, e2, a, b, c;
	TimeDelta = 0;

	//Check Node Spacing is ok
	CheckNodeSpacing();
	
	//Transform waves
	TransformWaves(OffshoreWavePeriod, OffshoreWaveDirection, OffshoreWaveHeight);

	//Find Shadows
	GetShadows(OffshoreWaveDirection);

	//cout << "Transformed Waves" << endl;

	//Compute Sediment Fluxes, find stable TimeDelta
	for (int i=0; i<NoNodes; ++i)
	{
		//CERC Equation
		if (Shadows[i] == 0) LongshoreFlux[i] = 86400.*0.41*pow(BreakingWaveHeight[i],2.5)*sin((M_PI/180.0)*BreakingWaveAngle[i])*cos((M_PI/180.0)*BreakingWaveAngle[i]);
		else LongshoreFlux[i] = 0;
		//Handle Boundary Conditions
		//PERIODIC BOUNDARY
		if (i==0 && StartBoundary == 1)
		{
			//CERC Equations for last node
			LongshoreFlux[NoNodes-2] = 86400.*0.41*pow(BreakingWaveHeight[NoNodes-2],2.5)*sin((M_PI/180.0)*BreakingWaveAngle[NoNodes-2])*cos((M_PI/180.0)*BreakingWaveAngle[NoNodes-2]);
			//Get Volume Change Vout - Vin - Vlost
			VolumeChange[i] = (LongshoreFlux[NoNodes-2]-LongshoreFlux[i])*(1.-LostFluxFraction);
		}
		else if (i==NoNodes-1 && StartBoundary == 1)
		{
			VolumeChange[i] = VolumeChange[0];
		}
		//FIXED BOUNDARY
		else if (i <= 1 && StartBoundary == 2)
		{
			//only allow transport towards start node
			if (LongshoreFlux[i] > 0) LongshoreFlux[i] = 0;
			//Get Volume Change
			VolumeChange[i] = 0;
		}
		//FIXED BOUNDARY //only allow transport towards end node 
		else if (i>=NoNodes-2 && EndBoundary == 2) VolumeChange[i] = 0;
		else VolumeChange[i] = (LongshoreFlux[i-1]-LongshoreFlux[i])*(1.-LostFluxFraction);
		
		//Find Stable Timestep (Courant-Friedrichs-Lewy Criterion)
		dt = (MeanNodeSpacing/abs(VolumeChange[i]));
		if (TimeDelta == 0) TimeDelta = dt;
		if (dt < TimeDelta) TimeDelta = dt;
	}
	
	//limit time step to MaxTimeStep (usually 1 day)
	if (TimeDelta > MaxTimeStep) TimeDelta = MaxTimeStep;
	cout << "dt = " << dt << endl;
	
	//Move the coast
	for (int i=0; i<NoNodes; ++i)
	{
		//Convert to Shoreline position change
		//check angles in here!

		//old position change for rectilinear grid
		//PositionChange[i] = VolumeChange[i]*cos((M_PI/180.0)*(Trend-Orientation[i]))/(CellWidth[i]*ClosureDepth);
		
		VolumeChange[i] = VolumeChange[i]*TimeDelta;
		
		//new position change for trapezium-shaped cells
		//a, b, c are parameters for quadratic equation solver
		if (i==0) e1 = Orientation[i]-FluxOrientation[i];
		else e1 = Orientation[i]-FluxOrientation[i-1];
		e2 = FluxOrientation[i]-Orientation[i];
		//one of these coefficients needs some sort of closure depth dependency
		a = 0.5*(tan((M_PI/180.)*e1) + tan((M_PI/180.)*e2));
		b = CellWidth[i];
		c = -VolumeChange[i]/ClosureDepth;

		//solve quadratic equation to get position change
		//if cell is rectilinear so dont need to solve the quadratic
		if (VolumeChange[i] != 0)
		{
			if (a==0) PositionChange[i] = VolumeChange[i]/(CellWidth[i]*ClosureDepth);
			else 
			{
				//whether you need the positive or the negative solution here will depend
				// on whether the shoreline is convex or concave?
				// pyhton script suggests this should always be negative
				//if (CONDITION HERE) PositionChange[i] = (-b+sqrt(b*b-4.*a*c))/(2.*a);
				//else PositionChange[i] = (-b-sqrt(b*b-4.*a*c))/(2.*a);
				if (Curvature[i] < 0) PositionChange[i] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
				else PositionChange[i] = (-b+sqrt((b*b)-(4.*a*c)))/(2.*a);
			}
		}
		else PositionChange[i] = 0;

		if (PositionChange[i] != PositionChange[i])
		{
			cout << "I am not equal to myself!!" << endl;
		}//effective angle for trend //more logic needed here
		// htis isnt needed for oreientation -| migration
		//Angle = 180.-Trend;

		//Decompose Position change into X and Y components
		//Going to need a bit more logic here to handle different orientations
		//current setup x change is negative, y change is positive
		X[i] -= PositionChange[i]*cos((M_PI/180.)*Orientation[i]);
		Y[i] += PositionChange[i]*sin((M_PI/180.)*Orientation[i]);
				
	}

	//Update Coastal Morphology
	CalculateMorphology();
}

void Coastline::KamphiusDiffusion()
{
	//Do sediment flux
	//Update morphology
}

#endif
