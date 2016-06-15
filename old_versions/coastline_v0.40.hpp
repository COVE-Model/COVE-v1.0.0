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
#include <queue>
#include "waveclimate_v0.04.hpp"

using namespace std;

#ifndef coastline_HPP
#define coastline_HPP

/**
* \class Coastline
* \brief Coastline class
*
* The coastline class object is the main object used in evolving a one-line coast. It consists of paired X, Y vectors which describe the location of the coastline
* and a series of functions with which to transform waves impinging on the coast, calculate alongshore sediment transport and update the coastline position during
* a model timestep. More to follow...
*
* MDH 6/1/14
*
*/
	
class Coastline
{
	private:
	//data members
	int NoNodes;										//Number of nodes along coastline
	vector<double> X;									//position in x (m) (private)
	vector<double> Y;									//position in y (m) (private)
	vector<double> Distance;						//Distance along the vector
	vector<double> Orientation;					//orientation/azimuth of shoreline across i-1, i+1 (private)
	vector<double> FluxOrientation;				//shoreline orientation bewteen i, i+1 (private)
	vector<double> LeftOrientation;			  //orientation for meshing left boundary
	vector<double> RightOrientation;			//orientation for meshing right boundary
	vector<double> CellWidth;						//Width of individual cells parallel to orientation
	vector<double> BreakingWaveHeight;			//Heights of breaking waves
	vector<double> BreakingWaveAngle;			//Angle of breaking waves
	vector<double> Shadows;							//Shadow Zone Vector
	vector<double> ShadowZoneWaveDirection;	//Wave direction modified for ref/diff in shadow zone
	vector<double> ShadowZoneWaveHeight;		//Wave height modified for ref/diff in shadow zone
	vector<double> LongshoreFlux;					//Volume of sediment transported alongshore in m3/day
	vector<double> VolumeChange;					//Volume change in each cell during a particular timestep (m3/day)
	vector<double> PositionChange;				//Magnitude of change in shoreline position
	vector<double> Volume;
	vector<double> Area;
	vector<double> e1;
	vector<double> e2;
	vector<double> XMidPoints;						//X position of downcoast midpoint to next node
	vector<double> YMidPoints;						//Y position of downcoast midpoint to next node
	vector<double> alpha_0;
	vector<double> X0;
	vector<double> Y0;
	vector<double> XL;
	vector<double> YL;
	vector<double> XR;
	vector<double> YR;
	vector<double> Dsf;
	
	vector<int> Fixed;						//Is coastline fixed?
	
	vector<int> NodesToRemove;				//For removing nodes that are getting squashed
	
	int StartBoundary;						//Boundary type periodic=1 /no flux=2 /continous flux=3 /sinks=4 (private)
	int EndBoundary;						//Boundary type periodic=1 /no flux=2 /continous flux=3 /sinks=4 (private)
	int MeanNodeSpacing;					//average spacing between nodes (m) (private)
	int DesiredNodeSpacing;					//for adding and removing nodes?
	int ShadowFlag;
	int TriangleFlag;
	int BuildCellsFlag;
	int NodeAddedFlag;
	
	double Trend;							//redundant?
	double TotalVolume;
	
	double OffshoreWavePeriod, OffshoreWaveHeight, OffshoreWaveDirection;
	double ClosureDepth;					// } This needs to 
	double ShorefaceSlope;					// }
	double BeachWidth;						// }
	
	//Use equals operator for assignment
	// CODE GOES HERE ? //
	void UpdateMorphology();
	void CalculateMorphology();
	void CheckNodeSpacing();
	void CalculateMeanNodeSpacing();
	void IntersectionAnalysis();
	void TransformWave(int i);
	void TransformWaves();
	void GetShadows();
	void RefractDiffractShadowZone();
	void BuildCellGeometries();
	
	//Use priority queue to deal with triangle nodes in order of increasing shoreface depth?
	struct CoastNode
	{
		double ShorefaceDepth;
		int i;
	};

	class CompareNodes
	{
		public:
		bool operator()(CoastNode& Node1, CoastNode& Node2)
		{
			if (Node1.ShorefaceDepth > Node2.ShorefaceDepth) return true;
			else return false;
		}	
	};

	priority_queue<CoastNode, vector<CoastNode>, CompareNodes> CoastlineQueue;
	vector< vector<CoastNode> > Recievers;
	
	//Structure to handle single nodes
	struct Node
	{
		double X, Y;
	};
	vector< vector<Node> > Vertices;
	
	void Initialise();
	void Initialise(string xyfilename);
	void Initialise(string xyfilename, float StartTime);
	void Initialise(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary);
	
	public:

	//temporarily public
	
	
	/*******************************************\
	| Coastline Object Initialisation Functions |
	\*******************************************/
	
	///Initialise
	Coastline()
	{
		Initialise();
	}
	
	///Initialise the coastline from xyfile.
	Coastline(string xyfilename)
	{
		Initialise(xyfilename);
	}
	
	///Initialise coastline from xyfile at a specific time.
	Coastline(string xyfilename, float StartTime)
	{
		Initialise(xyfilename, StartTime);
	}
	
	///Initialise coastline as a straight line.
	Coastline(int MeanNodeSpacing, double CoastLength, double Trend, int StartBoundary, int EndBoundary)
	{
		Initialise(MeanNodeSpacing, CoastLength, Trend, StartBoundary, EndBoundary);
	}
	
	///Read coastline from xyfile.
	void ReadCoast(string InputFileName, double Time);
	
	///Write coastline to xyfile.
	void WriteCoast(string OutputFileName, double Time);

	
	/*Flux Type 	0 = Simple Diffusion
					1 = CERC equation (default)
					2 = Kamphius equation
					3 = Bailard equation
					4 = Deigaard equation
	*/
	
	///Create priority queue for cell builder
	void SetupQueue();
	///Perform sediment transport and evolve coast
	void TransportSediment(double &TimeDelta, Wave TheWave, int FluxType=1, int RefDiffFlag=1, double FluxFraction=0, double LostFluxFraction=0);
	///Calculate volumetric sediment fluxes
	void CalculateFlux(int i, int FluxType);
	///Evolve case using simple diffusion model (no waves)
	void SimpleDiffusion(double DiffCoeff, double TimeDelta);
	///Perform quadratic solution for position change of a cell
	double SolveQuadratic(int i, double D_sf);
	///Perform cubic solution for position change of a cell
	double SolveCubic(int i, double D_sf);
	///Perform cubic solution for position change of a cell
	double SolveCubicNew(int i);
	
	/*****************************************\
	| Get Functions to return private members |
	\*****************************************/
	
	/**
	* \brief Return number of nodes in coastline.
	*
	* Function to return the private member NoNodes value, the number of nodes in the coastline object
	*
	* MDH 6/1/14
	*
	*/
	int get_NoNodes() const						{return NoNodes;}			//get Number of nodes along coastline
	/**
	* \brief Return coastline X-values.
	*
	* Function to return the private member X values
	*
	* MDH 6/1/14
	*
	*/
	vector<double> get_X() const				{return X;}					//get X vector position in x (m)
	/**
	* \brief Return coastline Y-values.
	*
	* Function to return the private member X values
	*
	* MDH 6/1/14
	*
	*/
	vector<double> get_Y() const				{return Y;}					//get Y vectir position in y (m)
	/**
	* \brief Return coastline orientation.
	*
	* Function to return the private member Orientation, the azimuth direction of the line connecting the two adjacent nodes to the node of interest (3-cell orientation). 
	* The evolving cosatal node advances or retreats perpendicular to this orientation.
	*
	* MDH 6/1/14
	*
	*/
	vector<double> get_Orientation() const 	 	{return Orientation;}		//get orientation of shoreline between i-1 and i+1
	/**
	* \brief Return coastline flux orientation.
	*
	* Function to return the private member FluxOrientation, the azimuth direction of the line connecting the node of interest to the next node in the vector (2-cell 
	* orientation). This FluxOrientation is used to transform waves and for determining the amount of alongshore sediment transport.
	*
	* MDH 6/1/14
	*
	*/
	vector<double> get_FluxOrientation() const 	{return FluxOrientation;}	//get shoreline orientation bewteen i and i+1
	/**
	* \brief Return coastal shadows.
	*
	* Function to return the private member Shadows, an integer vector which describes whether or not a node is in shadow, or casts a shadow, or is open coast.
	*
	* MDH 6/1/14
	*
	*/
	vector<double> get_Shadows() const			{return Shadows;}			//get Shadows	
	/**
	* \brief Return start boundary condition.
	*
	* Function to return the private member StartBoundary, the boundary condition at the start end of the vector. Currently, StartBoundary and EndBoundary must have the 
	* same value. 1 = Periodic boundary, 2 = Fixed Boundary (with or without sediment flux).
	*
	* MDH 6/1/14
	*
	*/
	int get_StartBoundary() const				{return StartBoundary;}		//get Boundary type start
	/**
	* \brief Return end boundary condition.
	*
	* Function to return the private member EndBoundary, the boundary condition at the start end of the vector. Currently, StartBoundary and EndBoundary must have the 
	* same value. 1 = Periodic boundary, 2 = Fixed Boundary (with or without sediment flux).
	*
	* MDH 6/1/14
	*
	*/
	int get_EndBoundary() const					{return EndBoundary;}		//get Boundary type end
	/**
	* \brief Return desired mean node spacing.
	*
	* Function to return the private member MeanNodeSpacing, the desired mean node spacing for the model.
	*
	* MDH 6/1/14
	*
	*/
	int get_MeanNodeSpacing() const				{return MeanNodeSpacing;}	//get average spacing between nodes (m)

};

#endif
