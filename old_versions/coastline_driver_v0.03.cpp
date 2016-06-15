#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "coastline_v0.03.hpp"
#include "waveclimate_v0.01.hpp"

using namespace std;


/*
coast always on the left of vector as you progress from i=0 to i=NoNodes
Node i=0 and i=NoNodes-1 are the same Node copied
*/

int main()
{
	//declare time shit
	int EndTime = 1000;					// End time (years)
	double Time = 0.;						// Start Time (years)
	double TimeDelta;						// Time step (years)
	double MaxTimeStep = 0.5;			// Maximum timestep (days)
	int PrintTimeDelta = 1;				// how often to print coastline (years)
	int PrintTime = PrintTimeDelta;	// Print time (years)
	
	//setup coastline
	//int MeanNodeSpacing = 100;			//mean spacing between nodes (m)
	//double CoastLength = 10000;		//length of coastline (m)
	//double Trend = 120.0;				//general orientation of coastline to initialise (degrees)
	//int StartBoundary = 1;				//boundary condition at first node (both = 1 for periodic)
	//int EndBoundary = 1;					//boundary condition at last node (both = 1 for periodic)
	
	//sediment transport parameters
	//double DiffCoeff = 0.01; 			//coastline diffusion coefficient m2/s (Limber & Murray 2011)
	
	//initialise coast
	//Coastline CoastVector1 = Coastline(); //test the empty
	//Coastline CoastVector = Coastline("coastlineXY_periodic.txt"); //load file with periodic boundary
	//Coastline CoastVector = Coastline("coastlineXY_fixed.txt"); //load file with fixed boundary
	//Coastline CoastVector = Coastline(MeanNodeSpacing, CoastLength, Trend, StartBoundary, EndBoundary);
	Coastline CoastVector = Coastline("coastline180.xy"); //load file with fixed boundary
	
	//declare wave stuff
		double OffshoreMeanWavePeriod = 6.0; //secs
		double OffshoreStDWavePeriod = 1.0;
		double OffshoreMeanWaveDirection = 40.0;
		double OffshoreStDWaveDirection = 20.0; //approaches from in degrees
		double OffshoreMeanWaveHeight = 1.0;
		double OffshoreStDWaveHeight = 0.2;	//metres
		double OffshoreWavePeriod = 6.0;
		double OffshoreWaveDirection = 100.0;
		double OffshoreWaveHeight = 0.8;
	
	//initialise wave climate
	//SingleWaveClimate WaveClimate = SingleWaveClimate(OffshoreWavePeriod, OffshoreWaveDirection, OffshoreWaveHeight);
	//test the empty
	GaussianWaveClimate WaveClimate = GaussianWaveClimate(OffshoreMeanWavePeriod, OffshoreStDWavePeriod,
						OffshoreMeanWaveDirection, OffshoreStDWaveDirection,
						OffshoreMeanWaveHeight, OffshoreStDWaveHeight);
		
	//test read and write functions here
	string WriteCoastFile = "CoastlineEvolution.xy";
	CoastVector.WriteCoast(WriteCoastFile, Time);

	//Simple diffusion here to start with
	//while (Time < EndTime)
	//{
//		Time += TimeDelta;
//		CoastVector.SimpleDiffusion(DiffCoeff, TimeDelta);
//		if (Time >= PrintTime)
//		{
//			CoastVector.WriteCoast(WriteCoastFile, Time);
//			PrintTime += PrintTimeDelta;
//		}
//	}
//	cout << endl << endl;

	//cout << "Starting Model Run" << endl;
	//cout << MaxTimeStep << endl;
	//Simple CERCDiffusion here to start with
	
	//ofstream WriteWaves;
	//WriteWaves.open("Waves.dat","w");
	
	while (Time < EndTime)
	{
		//get wave conditions
		OffshoreWavePeriod = WaveClimate.Get_WavePeriod();
		OffshoreWaveDirection = WaveClimate.Get_WaveDirection();
		OffshoreWaveHeight = WaveClimate.Get_WaveHeight();
		
		
		//Do sediment transport
		CoastVector.CERCDiffusion(TimeDelta, MaxTimeStep, OffshoreWavePeriod, OffshoreWaveDirection, OffshoreWaveHeight);
		Time += TimeDelta/365.;
		
		//cout << Time << endl;
		//CoastVector.WriteCoast(WriteCoastFile, PrintTime);
		
		//if (Time >= PrintTime)
		//{
			CoastVector.WriteCoast(WriteCoastFile, PrintTime);
			//PrintTime += PrintTimeDelta;
		//}
		
	}
	CoastVector.WriteCoast(WriteCoastFile, EndTime);
	cout 	<< endl << "Coastline: Results written to " << WriteCoastFile 
			<< endl << "Coastline: Model run complete." << endl << endl;
	//WriteWaves.close();
	exit(EXIT_FAILURE);
	
}


	//declare wave stuff
	//assume normal distribution pdfs of wave direction, height and period for now
	//sample from these distributions at random
	//double MeanOffshoreWaveHeight = 1.0;
	//double StDOffshoreWaveHeight = 0.2;
	//double MeanWavePeriod = 12.0; //secs
	//double StDWavePeriod = 2.0; //secs
	//double MeanOffshoreWaveDirection = 10.0; //approaches from in degrees
	//double StDOffshoreWaveDirection = 5.0; //degrees
	//double OffshoreWaveHeight, OffshoreWaveDirection;
			
	//declare coast morphology
	//double MaxWaveBase = 10.0;			//NB this could be dynamic based on the actual wave data
												//calculate explicitly the depth at which to start refraction
	//double ShoreFaceGradient = 0.02;	//Gradient of shoreface
	//double ShoreFaceDepth = 10.0;		//Depth of beach over which it is assumed sediment transport occurs
												//however there shouldn't be any sedi transport below MaxWaveBase
	
	//sediment transport parameters
	//double DiffCoeff = 0.01; 			//coastline diffusion coefficient m2/s (Limber & Murray 2011)
	
	
//	int pNoNodes = CoastVector.get_NoNodes();
//	vector<double> pX = CoastVector.get_X();
//	vector<double> pY = CoastVector.get_Y();
//	vector<double> pOrientation = CoastVector.get_Orientation();
//	vector<double> pSlope = CoastVector.get_Slope();
//	vector<double> pCurvature = CoastVector.get_Curvature();
//	int pStartBoundary = CoastVector.get_StartBoundary();
//	int pEndBoundary = CoastVector.get_EndBoundary();
//	int pMeanNodeSpacing = CoastVector.get_MeanNodeSpacing();
//	
//	int n = 2;
//	cout << pNoNodes << endl;		
//	cout << pX[n] << endl;				
//	cout << pY[n] << endl;				
//	cout << pOrientation[n] << endl;	
//	cout << pSlope[n] << endl;		
//	cout << pCurvature[n] << endl;	
//	cout << pStartBoundary << endl;	
//	cout << pEndBoundary << endl;	
//	cout << pMeanNodeSpacing << endl;
	
//	for (int i=0; i<NoNodes; ++i)
//	{
//		if (i==0)
//		{
//			CoastX[i] = 0;
//			CoastY[i] = 0;
//		}
//		else
//		{
//			//sort out random addition
//			CoastX[i] = CoastX[i-1] + MeanNodeSpacing*cos(TrendRads) + 2*(-0.5+((double)rand()/RAND_MAX));
//			CoastY[i] = CoastY[i-1] - MeanNodeSpacing*sin(TrendRads) + 2*(-0.5+((double)rand()/RAND_MAX));
//		}
//	}
//	
//	//output the initial vector
//	ofstream coastline_out;
//	coastline_out.open("coastlineXY.txt");
//	coastline_out << "EndTime " << EndTime << " TimeDelta " << TimeDelta << " NoNodes " << NoNodes << endl;
//	coastline_out << Time;
//	for (int i=0; i<NoNodes; ++i)	coastline_out << " " << CoastX[i];
//	coastline_out << endl << Time;
//	for (int i=0; i<NoNodes; ++i) coastline_out << " " << CoastY[i];
//	coastline_out << endl;
//	
//	//write waves to test
//	//ofstream waves_out;
//	//waves_out.open("waves.txt");
//	//waves_out << "dir hgt" << endl;
//	
//	//Loop through time and evolve the coast w/ periodic boundary
//	while (Time < EndTime)
//	{
//		//Update Time
//		Time += TimeDelta;
//		
//		//Print time to screeen
//		cout << flush << Time << "\r";
//		
//		//Get Waves
//		//Sample from PDFs
//		//Random normal distribution from Box-Muller transform sqrt(-2*log((double)rand()/RAND_MAX))*sin(2*PI*((double)rand()/RAND_MAX))
//		//double rand1 = (double)rand()/RAND_MAX;
//		//double rand2 = (double)rand()/RAND_MAX;
//		//OffshoreWaveHeight = MeanOffshoreWaveHeight + StDOffshoreWaveHeight*sqrt(-2*log(rand1))*sin(2*PI*(rand2));
//		//OffshoreWaveDirection = MeanOffshoreWaveDirection + StDOffshoreWaveDirection*sqrt(-2*log(rand1))*cos(2*PI*(rand2));
//					
//		//Transform Waves
//		//CODE HERE WHEN NEEDED
//		
//		// Calculate slope and curvature for periodic boundary coastline
//		// in this definition of periodic NoNode=0 and NoNode=NoNodes-1
//		// are essentially the same node but only used seperately for
//		// calculating slope and curvature? 
//			
////		//Loop through coastline nodes 
////		for (int i=0; i<NoNodes; ++i) 
////		{
////			//Calculate Slope
////			//Slope (orientation between i and i+1) is dy/dx
////			//first handle periodic boundary
////			//will need to add logic for different boundaries
////			if (i==NoNodes-1) Slope[i] = Slope[0];
////			else Slope[i] = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
////		}
////		//Loop through coastline nodes
////		double dX, dY, Slope0, Slope1; 
////		for (int i=0; i<NoNodes; ++i) 
////		{
////			//Calculate Curvature
////			//Slope (orientation between i and i+1) is dy/dx
////			//Curv is change in Slope over total distance between i-1 and i+1
////			//first handle periodic boundary
////			//will need to add logic for different boundaries
////			if (i==0) continue;
////			else if (i==NoNodes-1)
////			{
////				dX = CoastX[i]-CoastX[i-1]+CoastX[1]-CoastX[0];
////				dY = CoastY[i]-CoastY[i-1]+CoastY[1]-CoastY[0];
////				Slope0 = (CoastY[i]-CoastY[i-1])/(CoastX[i]-CoastX[i-1]);
////				Slope1 = (CoastY[1]-CoastY[0])/(CoastX[1]-CoastX[0]);
////				Curv[i] = (Slope1-Slope0)/sqrt(dX*dX + dY*dY);
////				Curv[0] = Curv[i];
////			}
////			else 
////			{
////				dX = (CoastX[i+1]-CoastX[i-1]);
////				dY = (CoastX[i+1]-CoastX[i-1]);
////				Slope0 = (CoastY[i]-CoastY[i-1])/(CoastX[i]-CoastX[i-1]);
////				Slope1 = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
////				Curv[i] = (Slope1-Slope0)/sqrt(dX*dX + dY*dY);
////			}
////		}
////		//calculate orientation
////		for (int i=0; i<NoNodes; ++i) 
////		{
////			if (i==0) 
////			{
////				dX = CoastX[NoNodes-1]-CoastX[NoNodes-2]+CoastX[1]-CoastX[0];
////				dY = CoastY[NoNodes-1]-CoastY[NoNodes-2]+CoastY[1]-CoastY[0];
////			}
////			else if (i==NoNodes-1)
////			{
////				dX = CoastX[i]-CoastX[i-1]+CoastX[1]-CoastX[0];
////				dY = CoastY[i]-CoastY[i-1]+CoastY[1]-CoastY[0];
////			}
////			else 
////			{
////				dX = (CoastX[i+1]-CoastX[i-1]);
////				dY = (CoastX[i+1]-CoastX[i-1]);
////			}
////			if (dX > 0) Orientation[i] = (PI*0.5 - atan((dY/dX)));
////			else if (dX < 0) Orientation[i] = (PI*1.5 - atan((dY/dX)));
////			else Orientation[i] = 0;
////		}
//		
//		for (int i=0; i<NoNodes; ++i) 
//		{
//			//Calculate Slope AND Curvature AND Orientation
//			//Slope (orientation between i and i+1) is dy/dx
//			//Curv is change in Slope over total distance between i-1 and i+1
//			//Orientation is the slope over total distance between i-1 and i+1
//			//first handle periodic boundary
//			//will need to add logic for different boundaries
//			if (i==0) Slope[i] = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
//			else if (i==NoNodes-1)
//			{
//				dX = CoastX[i]-CoastX[i-1]+CoastX[1]-CoastX[0];
//				dY = CoastY[i]-CoastY[i-1]+CoastY[1]-CoastY[0];
//				Slope[i] = Slope[0];
//				Curv[i] = (Slope[i]-Slope[i-1])/sqrt(dX*dX + dY*dY);
//				Curv[0] = Curv[i];
//				if (dX > 0) Orientation[i] = (PI*0.5 - atan((dY/dX)));
//				else if (dX < 0) Orientation[i] = (PI*1.5 - atan((dY/dX)));
//				Orientation[0] = Orientation[i];
//			}
//			else 
//			{
//				dX = (CoastX[i+1]-CoastX[i-1]);
//				dY = (CoastY[i+1]-CoastY[i-1]);
//				Slope[i] = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
//				Curv[i] = (Slope[i]-Slope[i-1])/sqrt(dX*dX + dY*dY);
//				if (dX > 0) Orientation[i] = (PI*0.5 - atan((dY/dX)));
//				else if (dX < 0) Orientation[i] = (PI*1.5 - atan((dY/dX)));
//			}
//		}
//		
//		//do diffusion
//		for (int i=0; i<NoNodes; ++i) 
//		{
//			//Diffuse coastline as a function of curvature
//			//Coastline adjusts perpendicular to orientation
//			Change = DiffCoeff*Curv[i]*TimeDelta;
//			//update X and Y
//			if (Orientation[i] < PI/2) 
//			{
//				alpha = Orientation[i];
//				CoastX[i] += Change*cos(alpha);
//				CoastY[i] -= Change*sin(alpha);
//			}
//			else if (Orientation[i] < PI)
//			{
//				alpha = PI-Orientation[i];
//				CoastX[i] -= Change*cos(alpha);
//				CoastY[i] -= Change*sin(alpha);
//			}
//			else if (Orientation[i] < 3*PI/2)
//			{
//				alpha = Orientation[i]-PI;
//				CoastX[i] -= Change*cos(alpha);
//				CoastY[i] += Change*sin(alpha);
//			}
//			else
//			{
//				alpha = Orientation[i]-PI;
//				CoastX[i] += Change*cos(alpha);
//				CoastY[i] += Change*sin(alpha);
//			}
//		}
//		
//		//Loop from wave base to calculate breaking waveheight?
//		WaterDepth = WaveBase;
//		
//		while (WaveBreak == 0)
//		{
//			//Calculate offshore wave properties
//			DeepWaterWaveVelocity = g*WavePeriod/(2*M_PI);
//			DeepWaterWaveLength = DeepWaterWaveVelocity*WavePeriod;
//			
//			//Wavelength in intermediate water depths from Fenton and Mackee
//			WaveLength = DeepWaterWaveLength * pow(tanh(pow(pow(2*M_PI/Period,2)*WaterDepth/g,.75)),(2/3));
//			
//			//Wave Velocity from WaveLength and Period
//			WaveVelocity = WaveLength/Period;
//			
//			//Calculate Shoaled Wave Angle
//			Angle = asin(WaveVelocity/DeepWaterWaveVelocity * sin(DeepWaterWaveAngle));
//			
//			// Determine n = 1/2(1+2kh/tanh(kh)) Komar 5.21 //
//			/* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L		*/
//			k = 2*M_PI/WaveLength;
//			n =0.5*(1+(2*k*WaterDepth/sinh(2*k*WaterDepth)));
//			
//			// Determine Wave height from refract calcs - Komar 5.49
//			//	Calculate Ks and Kr!!!???
//			WaveHeight = OffShoreWavHeight * pow(DeepWaterWaveVelocity*cos(DeepWaterWaveAngle)/(WaveVelocity*2.0*n*cos(Angle)),0.5);
//			
//			//Check wave breaking condition and update water depth
//			if (WvHeight > Depth*KBreak || Depth <= RefractStep) break;
//			else Depth -= DepthIncrement;
//		}
//			
//			//
//			// write out derivation of diffusion equation
//			//
//			
//			//Or simply assume H = 0.78*Ho
//			//Is there a simplified refraction we can do? 0.78*alphao
//			
//		//Compute Sediment Flux
//			//Need shoreline orientation
//			//CERC Eq.
//			
//		//Update Coast	
//			//change x/y
//			//resample? at even spaces?
//			//recalc orientation?
//			
//		//Output Vectors
//		coastline_out << Time;
//		for (int i=0; i<NoNodes; ++i)	coastline_out << " " << CoastX[i];
//		coastline_out << endl << Time;
//		for (int i=0; i<NoNodes; ++i) coastline_out << " " << CoastY[i];
//		coastline_out << endl;
//		//output x,y,qs?
//	}
//	cout << endl;
//	//waves_out.close();
//	coastline_out.close();
//}
