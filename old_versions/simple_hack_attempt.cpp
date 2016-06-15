#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


/*
coast always on the left of vector as you progress from i=0 to i=NoNodes
Node i=0 and i=NoNodes-1 are the same Node copied
*/

int main()
{
	//declare shit
	static const double PI = 3.14159;
	static const double EU = 2.71828;
	
	int NoNodes = 101.0;
	double MeanNodeSpacing = 10.0;
	double Trend = 120.0;

	double TrendRads = (180-Trend)*(PI/180);
	double alpha, dX, dY, Change;
	vector<double> CoastX(NoNodes);
	vector<double> CoastY(NoNodes);
	vector<double> Orientation(NoNodes);
	vector<double> Slope(NoNodes);
	vector<double> Curv(NoNodes);

	//declare time shit
	int Time = 0;
	int EndTime = 1000;
	int TimeDelta = 1;
	
	//declare wave stuff
	//assume normal distribution pdfs of wave direction, height and period for now
	//sample from these distributions at random
	double MeanOffshoreWaveHeight = 1.0;
	double StDOffshoreWaveHeight = 0.2;
	double MeanWavePeriod = 12.0; //secs
	double StDWavePeriod = 2.0; //secs
	double MeanOffshoreWaveDirection = 10.0; //approaches from in degrees
	double StDOffshoreWaveDirection = 5.0; //degrees
	double OffshoreWaveHeight, OffshoreWaveDirection;
			
	//declare coast morphology
	double MaxWaveBase = 10.0;			//NB this could be dynamic based on the actual wave data
												//calculate explicitly the depth at which to start refraction
	double ShoreFaceGradient = 0.02;	//Gradient of shoreface
	double ShoreFaceDepth = 10.0;		//Depth of beach over which it is assumed sediment transport occurs
												//however there shouldn't be any sedi transport below MaxWaveBase
	
	//sediment transport parameters
	double DiffCoeff = 0.01; 			//coastline diffusion coefficient m2/s (Limber & Murray 2011)
													
	//initialise coast
	for (int i=0; i<NoNodes; ++i)
	{
		if (i==0)
		{
			CoastX[i] = 0;
			CoastY[i] = 0;
		}
		else
		{
			//sort out random addition
			CoastX[i] = CoastX[i-1] + MeanNodeSpacing*cos(TrendRads) + 2*(-0.5+((double)rand()/RAND_MAX));
			CoastY[i] = CoastY[i-1] - MeanNodeSpacing*sin(TrendRads) + 2*(-0.5+((double)rand()/RAND_MAX));
		}
	}
	
	//output the initial vector
	ofstream coastline_out;
	coastline_out.open("coastlineXY.txt");
	coastline_out << "EndTime " << EndTime << " TimeDelta " << TimeDelta << " NoNodes " << NoNodes << endl;
	coastline_out << Time;
	for (int i=0; i<NoNodes; ++i)	coastline_out << " " << CoastX[i];
	coastline_out << endl << Time;
	for (int i=0; i<NoNodes; ++i) coastline_out << " " << CoastY[i];
	coastline_out << endl;
	
	//write waves to test
	//ofstream waves_out;
	//waves_out.open("waves.txt");
	//waves_out << "dir hgt" << endl;
	
	//Loop through time and evolve the coast w/ periodic boundary
	while (Time < EndTime)
	{
		//Update Time
		Time += TimeDelta;
		
		//Print time to screeen
		cout << flush << Time << "\r";
		
		//Get Waves
		//Sample from PDFs
		//Random normal distribution from Box-Muller transform sqrt(-2*log((double)rand()/RAND_MAX))*sin(2*PI*((double)rand()/RAND_MAX))
		//double rand1 = (double)rand()/RAND_MAX;
		//double rand2 = (double)rand()/RAND_MAX;
		//OffshoreWaveHeight = MeanOffshoreWaveHeight + StDOffshoreWaveHeight*sqrt(-2*log(rand1))*sin(2*PI*(rand2));
		//OffshoreWaveDirection = MeanOffshoreWaveDirection + StDOffshoreWaveDirection*sqrt(-2*log(rand1))*cos(2*PI*(rand2));
					
		//Transform Waves
		//CODE HERE WHEN NEEDED
		
		// Calculate slope and curvature for periodic boundary coastline
		// in this definition of periodic NoNode=0 and NoNode=NoNodes-1
		// are essentially the same node but only used seperately for
		// calculating slope and curvature? 
			
//		//Loop through coastline nodes 
//		for (int i=0; i<NoNodes; ++i) 
//		{
//			//Calculate Slope
//			//Slope (orientation between i and i+1) is dy/dx
//			//first handle periodic boundary
//			//will need to add logic for different boundaries
//			if (i==NoNodes-1) Slope[i] = Slope[0];
//			else Slope[i] = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
//		}
//		//Loop through coastline nodes
//		double dX, dY, Slope0, Slope1; 
//		for (int i=0; i<NoNodes; ++i) 
//		{
//			//Calculate Curvature
//			//Slope (orientation between i and i+1) is dy/dx
//			//Curv is change in Slope over total distance between i-1 and i+1
//			//first handle periodic boundary
//			//will need to add logic for different boundaries
//			if (i==0) continue;
//			else if (i==NoNodes-1)
//			{
//				dX = CoastX[i]-CoastX[i-1]+CoastX[1]-CoastX[0];
//				dY = CoastY[i]-CoastY[i-1]+CoastY[1]-CoastY[0];
//				Slope0 = (CoastY[i]-CoastY[i-1])/(CoastX[i]-CoastX[i-1]);
//				Slope1 = (CoastY[1]-CoastY[0])/(CoastX[1]-CoastX[0]);
//				Curv[i] = (Slope1-Slope0)/sqrt(dX*dX + dY*dY);
//				Curv[0] = Curv[i];
//			}
//			else 
//			{
//				dX = (CoastX[i+1]-CoastX[i-1]);
//				dY = (CoastX[i+1]-CoastX[i-1]);
//				Slope0 = (CoastY[i]-CoastY[i-1])/(CoastX[i]-CoastX[i-1]);
//				Slope1 = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
//				Curv[i] = (Slope1-Slope0)/sqrt(dX*dX + dY*dY);
//			}
//		}
//		//calculate orientation
//		for (int i=0; i<NoNodes; ++i) 
//		{
//			if (i==0) 
//			{
//				dX = CoastX[NoNodes-1]-CoastX[NoNodes-2]+CoastX[1]-CoastX[0];
//				dY = CoastY[NoNodes-1]-CoastY[NoNodes-2]+CoastY[1]-CoastY[0];
//			}
//			else if (i==NoNodes-1)
//			{
//				dX = CoastX[i]-CoastX[i-1]+CoastX[1]-CoastX[0];
//				dY = CoastY[i]-CoastY[i-1]+CoastY[1]-CoastY[0];
//			}
//			else 
//			{
//				dX = (CoastX[i+1]-CoastX[i-1]);
//				dY = (CoastX[i+1]-CoastX[i-1]);
//			}
//			if (dX > 0) Orientation[i] = (PI*0.5 - atan((dY/dX)));
//			else if (dX < 0) Orientation[i] = (PI*1.5 - atan((dY/dX)));
//			else Orientation[i] = 0;
//		}
		
		for (int i=0; i<NoNodes; ++i) 
		{
			//Calculate Slope AND Curvature AND Orientation
			//Slope (orientation between i and i+1) is dy/dx
			//Curv is change in Slope over total distance between i-1 and i+1
			//Orientation is the slope over total distance between i-1 and i+1
			//first handle periodic boundary
			//will need to add logic for different boundaries
			if (i==0) Slope[i] = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
			else if (i==NoNodes-1)
			{
				dX = CoastX[i]-CoastX[i-1]+CoastX[1]-CoastX[0];
				dY = CoastY[i]-CoastY[i-1]+CoastY[1]-CoastY[0];
				Slope[i] = Slope[0];
				Curv[i] = (Slope[i]-Slope[i-1])/sqrt(dX*dX + dY*dY);
				Curv[0] = Curv[i];
				if (dX > 0) Orientation[i] = (PI*0.5 - atan((dY/dX)));
				else if (dX < 0) Orientation[i] = (PI*1.5 - atan((dY/dX)));
				Orientation[0] = Orientation[i];
			}
			else 
			{
				dX = (CoastX[i+1]-CoastX[i-1]);
				dY = (CoastY[i+1]-CoastY[i-1]);
				Slope[i] = (CoastY[i+1]-CoastY[i])/(CoastX[i+1]-CoastX[i]);
				Curv[i] = (Slope[i]-Slope[i-1])/sqrt(dX*dX + dY*dY);
				if (dX > 0) Orientation[i] = (PI*0.5 - atan((dY/dX)));
				else if (dX < 0) Orientation[i] = (PI*1.5 - atan((dY/dX)));
			}
		}
		
		//do diffusion
		for (int i=0; i<NoNodes; ++i) 
		{
			//Diffuse coastline as a function of curvature
			//Coastline adjusts perpendicular to orientation
			Change = DiffCoeff*Curv[i]*TimeDelta;
			//update X and Y
			if (Orientation[i] < PI/2) 
			{
				alpha = Orientation[i];
				CoastX[i] += Change*cos(alpha);
				CoastY[i] -= Change*sin(alpha);
			}
			else if (Orientation[i] < PI)
			{
				alpha = PI-Orientation[i];
				CoastX[i] -= Change*cos(alpha);
				CoastY[i] -= Change*sin(alpha);
			}
			else if (Orientation[i] < 3*PI/2)
			{
				alpha = Orientation[i]-PI;
				CoastX[i] -= Change*cos(alpha);
				CoastY[i] += Change*sin(alpha);
			}
			else
			{
				alpha = Orientation[i]-PI;
				CoastX[i] += Change*cos(alpha);
				CoastY[i] += Change*sin(alpha);
			}
		}
		
		//Loop from wave base to calculate breaking waveheight?
		WaterDepth = WaveBase;
		
		while (WaveBreak == 0)
		{
			//Calculate offshore wave properties
			DeepWaterWaveVelocity = g*WavePeriod/(2*M_PI);
			DeepWaterWaveLength = DeepWaterWaveVelocity*WavePeriod;
			
			//Wavelength in intermediate water depths from Fenton and Mackee
			WaveLength = DeepWaterWaveLength * pow(tanh(pow(pow(2*M_PI/Period,2)*WaterDepth/g,.75)),(2/3));
			
			//Wave Velocity from WaveLength and Period
			WaveVelocity = WaveLength/Period;
			
			//Calculate Shoaled Wave Angle
			Angle = asin(WaveVelocity/DeepWaterWaveVelocity * sin(DeepWaterWaveAngle));
			
			// Determine n = 1/2(1+2kh/tanh(kh)) Komar 5.21 //
			/* First Calculate kh = 2 pi Depth/L  from k = 2 pi/L		*/
			k = 2*M_PI/WaveLength;
			n =0.5*(1+(2*k*WaterDepth/sinh(2*k*WaterDepth)));
			
			// Determine Wave height from refract calcs - Komar 5.49
			//	Calculate Ks and Kr!!!???
			WaveHeight = OffShoreWavHeight * pow(DeepWaterWaveVelocity*cos(DeepWaterWaveAngle)/(WaveVelocity*2.0*n*cos(Angle)),0.5);
			
			//Check wave breaking condition and update water depth
			if (WvHeight > Depth*KBreak || Depth <= RefractStep) break;
			else Depth -= DepthIncrement;
		}
			
			//
			// write out derivation of diffusion equation
			//
			
			//Or simply assume H = 0.78*Ho
			//Is there a simplified refraction we can do? 0.78*alphao
			
		//Compute Sediment Flux
			//Need shoreline orientation
			//CERC Eq.
			
		//Update Coast	
			//change x/y
			//resample? at even spaces?
			//recalc orientation?
			
		//Output Vectors
		coastline_out << Time;
		for (int i=0; i<NoNodes; ++i)	coastline_out << " " << CoastX[i];
		coastline_out << endl << Time;
		for (int i=0; i<NoNodes; ++i) coastline_out << " " << CoastY[i];
		coastline_out << endl;
		//output x,y,qs?
	}
	cout << endl;
	//waves_out.close();
	coastline_out.close();
}
