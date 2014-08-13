/** \file TimingInformation.cpp
 *\brief Structures and methods for high res timing
 *
 *Contains common data structures and methods needed for high 
 *resolution timing analysis
 *
 *\author S. V. Paulauskas 
 *\date 09 May 2011
 */
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include <cmath>
#define PI 3.14159265358979323846

#include "PathHolder.hpp"
#include "RawEvent.hpp"
#include "TimingInformation.hpp"
#include "Trace.hpp"

using namespace std;

map<string, double> TimingInformation::constantsMap;
TimingInformation::TimingCalMap TimingInformation::calibrationMap;

//********** Data (Default)**********
TimingInformation::TimingData::TimingData(void) : trace(emptyTrace)
{
    aveBaseline    = numeric_limits<double>::quiet_NaN();
    discrimination = numeric_limits<double>::quiet_NaN();
    highResTime    = numeric_limits<double>::quiet_NaN();
    maxpos         = numeric_limits<double>::quiet_NaN();
    maxval         = numeric_limits<double>::quiet_NaN();
    phase          = numeric_limits<double>::quiet_NaN();
    snr            = numeric_limits<double>::quiet_NaN();
    stdDevBaseline = numeric_limits<double>::quiet_NaN();
    tqdc           = numeric_limits<double>::quiet_NaN();
    walk           = numeric_limits<double>::quiet_NaN();
    walkCorTime    = numeric_limits<double>::quiet_NaN();
    
    numAboveThresh = -1;
}

//********** Data **********
TimingInformation::TimingData::TimingData(ChanEvent *chan) : trace(chan->GetTrace())
{
    //put all the times as ns
    // 11/18/2012 KM: after removal of Trace member, trace is
    // refered here directly. Should not change anything but allows to 
    // remove global variable emptyTrace
    const Trace& trace = chan->GetTrace();
    aveBaseline    = trace.GetValue("baseline");
    discrimination = trace.GetValue("discrim");
    highResTime    = chan->GetHighResTime()*1e+9;  
    maxpos         = trace.GetValue("maxpos");
    maxval         = trace.GetValue("maxval");
    numAboveThresh = trace.GetValue("numAboveThresh");
    phase          = trace.GetValue("phase")*(pixie::adcClockInSeconds*1e+9);
    stdDevBaseline = trace.GetValue("sigmaBaseline");
    tqdc           = trace.GetValue("tqdc")/qdcCompression;
    walk           = trace.GetValue("walk");
        
    //Calculate some useful quantities.
    //snr = pow(maxval/stdDevBaseline,2); 
    snr = 20*log10(maxval/stdDevBaseline); 
    walkCorTime   = highResTime - walk;

    //validate data and set a flag saying it's ok
    // clean up condition at some point
    if((maxval == maxval) && (phase == phase) && 
       (tqdc == tqdc) && (highResTime == highResTime) &&
       (stdDevBaseline == stdDevBaseline))
		dataValid = true;
    else
		dataValid = false;
}

#ifdef useroot
//********** DataRoot **********
TimingInformation::DataRoot::DataRoot(void) //--- initiallizes DataRoot variable arrays
{
  
 /*   for (size_t i = 0; i < maxMultiplicity; i++) {
	tof[i]    	= numeric_limits<double>::quiet_NaN();
	lqdc[i] 	= numeric_limits<double>::quiet_NaN();
	rqdc[i]    	= numeric_limits<double>::quiet_NaN();
	tsLow[i]        = numeric_limits<double>::quiet_NaN();
	tsHigh[i]       = numeric_limits<double>::quiet_NaN();
	lMaxVal[i]	= numeric_limits<double>::quiet_NaN();
	rMaxVal[i]	= numeric_limits<double>::quiet_NaN();
	qdc[i]		= numeric_limits<double>::quiet_NaN();
	energy[i]	= numeric_limits<double>::quiet_NaN();
//	location[i]     = numeric_limits<int>::quiet_NaN();
	location[i]     = -9;
    }
*/
	tof		= numeric_limits<double>::quiet_NaN();
	lqdc		= numeric_limits<double>::quiet_NaN();
	rqdc		= numeric_limits<double>::quiet_NaN();
	tsLow		= numeric_limits<double>::quiet_NaN();
	tsHigh		= numeric_limits<double>::quiet_NaN();
	lMaxVal		= numeric_limits<double>::quiet_NaN();
	rMaxVal		= numeric_limits<double>::quiet_NaN();
	qdc		= numeric_limits<double>::quiet_NaN();
	energy   	= numeric_limits<double>::quiet_NaN();
	ejectAngle  	= numeric_limits<double>::quiet_NaN();
	recoilAngle 	= numeric_limits<double>::quiet_NaN();
	recoilEnergy 	= numeric_limits<double>::quiet_NaN();
	exciteEnergy 	= numeric_limits<double>::quiet_NaN(); 
	flightPath	= numeric_limits<double>::quiet_NaN();
	xflightPath	= numeric_limits<double>::quiet_NaN();
	yflightPath 	= numeric_limits<double>::quiet_NaN();
	zflightPath 	= numeric_limits<double>::quiet_NaN();
    
    	multiplicity = 0;
   	dummy = numeric_limits<int>::quiet_NaN();
    	location = numeric_limits<int>::quiet_NaN();

}
#endif

//********** BarData **********
TimingInformation::BarData::BarData(const TimingData &Right, const TimingData &Left, const TimingCal &cal, const string &type) 
{
    //Clear the maps that hold this information
    timeOfFlight.clear();
    energy.clear();
    corTimeOfFlight.clear();
      
    //Set the values for the useful bar stuff. 
    lMaxVal   = Left.maxval;
    rMaxVal   = Right.maxval;
    lqdc      = Left.tqdc;
    rqdc      = Right.tqdc;
    lTime     = Left.highResTime;
    rTime     = Right.highResTime;
    qdc       = sqrt(Right.tqdc*Left.tqdc);
    qdcPos    = (lqdc - rqdc) / qdc;
    timeAve   = (Right.highResTime + Left.highResTime)*0.5; //in ns
    timeDiff  = (Left.highResTime-Right.highResTime) + cal.lrtOffset; //in ns
    walkCorTimeDiff = (Left.walkCorTime-Right.walkCorTime) + cal.lrtOffset; 
    walkCorTimeAve  = (Left.walkCorTime+Right.walkCorTime)*0.5; //in ns

    //Calculate some useful quantities for the bar analysis.
    event = BarEventCheck(timeDiff, type);
    flightPath = CalcFlightPath(timeDiff, cal, type, xflightPath, yflightPath, zflightPath);

};

//********** VMLData **********
TimingInformation::vmlData::vmlData(const BarData bar, double corTOF, double enrgy, double tLow, double tHigh, double recoilE)
{
    //Set the values for the useful bar stuff. 
    tof	     	 	= corTOF;
    lqdc      		= bar.lqdc;
    rqdc      		= bar.rqdc;
    tsLow     		= tLow;
    tsHigh    		= tHigh;
    lMaxVal   		= bar.lMaxVal;
    rMaxVal   		= bar.rMaxVal;
    qdc       		= bar.qdc;
    energy   		= enrgy;			//already in keV
    ejectAngle 		= bar.ejectAngle*180/PI;	//Conversion from radians->Degrees
    recoilEnergy 	= recoilE*1000;			//Energies conversion to keV
    recoilAngle 	= bar.recoilAngle*180/PI;
    exciteEnergy 	= bar.exciteEnergy*1000;
    flightPath		= bar.flightPath;
    xflightPath		= bar.xflightPath;
    yflightPath		= bar.yflightPath;
    zflightPath		= bar.zflightPath;
    	
    dummy 		= 1;
};

//********** BarEventCheck **********
bool TimingInformation::BarData::BarEventCheck(const double &timeDiff, const string &type)
{
    if(type == "small") {
		double lengthSmallTime = TimingInformation::GetConstant("lengthSmallTime");
		return(fabs(timeDiff) < lengthSmallTime+20);
    } 
	else if(type == "big") {
		double lengthBigTime = TimingInformation::GetConstant("lengthBigTime");
		return(fabs(timeDiff) < lengthBigTime+20);
    } 
	else
		return(false);
}


//********** CalcFlightPath **********
double TimingInformation::BarData::CalcFlightPath(double &timeDiff, const TimingCal &cal, const string &type, 
								double &xflightPath, double &yflightPath, double &zflightPath)
{
	if ( fabs(timeDiff) < 16 ){//maximum amount of time for time difference between signal, otherwise bad event
		
		double rad = PI/180;
		double speedOfLightInBar = 0.0;
		if(type == "small"){
			speedOfLightInBar = TimingInformation::GetConstant("speedOfLightSmall");
		}
		else if(type == "big"){
			speedOfLightInBar = TimingInformation::GetConstant("speedOfLightBig");
		}
		else{
			return(numeric_limits<double>::quiet_NaN());
		}
    		
		//Calculate components of flight path from center of bar
		xflightPath = cal.x - 0.5*speedOfLightInBar*timeDiff*sin(cal.orientTheta*rad)*cos(cal.orientPhi*rad);
		yflightPath = cal.y - 0.5*speedOfLightInBar*timeDiff*sin(cal.orientTheta*rad)*sin(cal.orientPhi*rad);
		zflightPath = cal.z - 0.5*speedOfLightInBar*timeDiff*cos(cal.orientTheta*rad);	
		return(sqrt( pow(xflightPath,2) + pow(yflightPath,2) + pow(zflightPath,2) ));
	
	}
	return(-1);//not a correct event that we want to calculate, so set to -1
}

//********** CalculateEnergy **********
double TimingInformation::CalcEnergy(const double &corTOF, const double &r)
{	
	//calculates energy of the ejected particle (neutron)
	double speedOfLight = TimingInformation::GetConstant("speedOfLight");
   	double neutronMass  = TimingInformation::GetConstant("neutronMass");

   	if (corTOF > 5) //do not calculate neutron energies for gammas
    		return((0.5*neutronMass*pow((r/corTOF)/speedOfLight, 2))*1000);

	return(0);
}

//********** CalculateRecoilEnergy *********
double TimingInformation::CalcRecoilEnergy(const double &energy, const double &flightPath, double &zflightPath, double &ejectAngle, double &recoilAngle, double &exciteEnergy)
{
    /* There are four particles involved in calculating the energy of the recoiling particle:
	m1 is the beam particle, m2 is the target, m3 is the ejectile that is measured
	and m4 is the particle whose energy and angle we would like to calculate. 
    */
	
    //read in constants
	double m1 = TimingInformation::GetConstant("alphaMass");	//masses in MeV/c^2
	double m2 = TimingInformation::GetConstant("carbon13Mass");
	double m3 = TimingInformation::GetConstant("neutronMass");
	double m4 = TimingInformation::GetConstant("oxygen16Mass");
	double Ebeam = TimingInformation::GetConstant("beamEnergy");	//beam energy in MeV

    //calculate needed quantities
    	double ejectEnergy = energy/1000; 		//energy of eject in MeV
	double Etotal = (Ebeam + m1) + m2;		//if m2 is stationary target
	double E1 = Ebeam + m1;
	double E3 = ejectEnergy + m3; 			//ejectile energy (neutron)
	double P1 = sqrt(Ebeam*Ebeam + 2*m1*Ebeam);	//relativistic momentum for beam and eject
	double P3 = sqrt(ejectEnergy*ejectEnergy + 2*m3*ejectEnergy); 
	ejectAngle = acos(zflightPath/flightPath); 	//angle from z-axis to ejected neutron path
	recoilAngle = atan( (-P3*sin(ejectAngle))/ (P1-P3*cos(ejectAngle)) ) + PI/2; // angle from z-axis to recoil path

    //calculate Qs
    	double Qconst = m1+m2-m3-m4;
	double Qreact = m1+m2-m3-sqrt(m1*m1 + m2*m2 + m3*m3 + 2*m2*E1 - 2*E3*(E1+m2) + 2*P1*P3*cos(ejectAngle));

	exciteEnergy = Qconst-Qreact;	//excitation energy of recoil in MeV

    //calulate Recoil energy
	double T4 = Etotal -E3 - (m1+m2-m3-Qreact); //Kinetic Energy of Recoil

	return(T4);
	//return(T4+m4); //total energy of Recoil
}
 
//********** GetConstant **********
double TimingInformation::GetConstant(const string &name)
{
    map<string, double>::iterator itTemp = 
	constantsMap.find(name);
    if(itTemp == constantsMap.end()) {
	cout << endl << endl 
	     << "Cannot Locate " << name << " in the Timing Constants Map!!" 
	     << endl << "Please check timingConstants.txt" << endl << endl;
	exit(EXIT_FAILURE); 
    } else {
	double value = (*constantsMap.find(name)).second;
	return(value);
    }
    return (numeric_limits<double>::quiet_NaN());
}


//********** GetTimingCalParameter **********
TimingInformation::TimingCal TimingInformation::GetTimingCal(const IdentKey &identity)
{
    map<IdentKey, TimingCal>::iterator itTemp = 
	calibrationMap.find(identity);
    if(itTemp == calibrationMap.end()) {
	cout << endl << endl 
	     << "Cannot locate detector named " << identity.second 
	     << " at location " << identity.first 
	     << "in the Timing Calibration!!" 
	     << endl << "Please check timingCal.txt" << endl << endl;
	exit(EXIT_FAILURE); 
    } else {
	TimingCal value = (*calibrationMap.find(identity)).second;
	return(value);
    }
}


//********** ReadTimingConstants **********
void TimingInformation::ReadTimingConstants(void)
{
    PathHolder* conf_path = new PathHolder();
    string constantsFileName = conf_path->GetFullPath("timingConstants_5300_13C.txt");
    delete conf_path;

    ifstream readConstants(constantsFileName.c_str());
    
    if (!readConstants) {
        cout << endl << "Cannot open file 'timingConstants.txt'" 
	     << "-- This is Fatal! Exiting..." << endl << endl;
	exit(EXIT_FAILURE);
    } else {
	while(readConstants) {
	    double value = 0.0;
	    string name = "";
	    
	    if (isdigit(readConstants.peek())) {
		readConstants >> value  >> name;		    
		constantsMap.insert(make_pair(name, value));
	    } else{
		readConstants.ignore(1000, '\n');
	    }
	} // end while (!readConstants) loop 
    }
    readConstants.close();

    double lengthSmallTime = 
	(*constantsMap.find("lengthSmallPhysical")).second /
	(*constantsMap.find("speedOfLightSmall")).second;
    double lengthBigTime   = 
	(*constantsMap.find("lengthBigPhysical")).second /
	(*constantsMap.find("speedOfLightBig")).second;

    constantsMap.insert(make_pair("lengthBigTime", lengthBigTime));
    constantsMap.insert(make_pair("lengthSmallTime", lengthSmallTime));
} //void TimingInformation::ReadTimingConstants


//********** ReadTimingCalibration **********
void TimingInformation::ReadTimingCalibration(void)
{
    TimingCal timingcal;
    PathHolder* conf_path = new PathHolder();
    string timeCalFileName = conf_path->GetFullPath("timingCal_5300_13C_oct14_1.txt"); 
    delete conf_path;

    ifstream timingCalFile(timeCalFileName.c_str());
   
    if (!timingCalFile) {
        cout << endl << "Cannot open file 'timingCal.txt'" 
	     << "-- This is Fatal! Exiting..." << endl << endl;
		exit(EXIT_FAILURE);
    } 
    else {
	double rad = PI/180;

	while(timingCalFile) {
	    if (isdigit(timingCalFile.peek())) {
		unsigned int location = -1;
		string type = "";

		timingCalFile >> location >> type >> timingcal.x
			      >> timingcal.y >> timingcal.z >> timingcal.r 
			      >> timingcal.barPosTheta >> timingcal.barPosPhi
			      >> timingcal.orientTheta >> timingcal.orientPhi 
			      >> timingcal.lrtOffset
			      >> timingcal.tofOffset0 >> timingcal.tofOffset1;

		//Coordinate Conversions
			if ( (timingcal.x!=0 || timingcal.y!=0 || timingcal.z!=0) && timingcal.r!=0){ //error -- only specify one coor. system
				cout << endl << "ERROR--Specify only one position coordinate system for bar #" << location << endl;
			}
			else if (timingcal.x==0 && timingcal.y==0 && timingcal.z==0 && timingcal.r!=0 ){//from spherical to cartesian
				
				timingcal.x = timingcal.r*sin(timingcal.barPosTheta*rad)*cos(timingcal.barPosPhi*rad);
				timingcal.y = timingcal.r*sin(timingcal.barPosTheta*rad)*sin(timingcal.barPosPhi*rad);
				timingcal.z = timingcal.r*cos(timingcal.barPosTheta*rad);
				cout << location << " : " << "x = " << timingcal.x << ", y = " << timingcal.y << ", z = " << timingcal.z << endl;
			   	
			}
			else if ( (timingcal.r==0 && timingcal.barPosTheta==0 && timingcal.barPosPhi==0) && 
					(timingcal.x!=0 || timingcal.y!=0 || timingcal.z!=0) ){
				//from cartesian to spherical
				timingcal.r = sqrt(pow(timingcal.x,2)+pow(timingcal.y,2)+pow(timingcal.z,2));
				timingcal.barPosTheta = acos(timingcal.z/timingcal.r)*180/PI;
				timingcal.barPosPhi = atan2(timingcal.y,timingcal.x)*180/PI; //specifies correct quadrant
				cout << location <<" : radius : " << timingcal.r
					<< ",  Theta : " << timingcal.barPosTheta
					<< ",  Phi : " << timingcal.barPosPhi << endl;
			}
			else{};			

		IdentKey calKey(location, type);
		calibrationMap.insert(make_pair(calKey, timingcal));

	    } else{
			timingCalFile.ignore(1000, '\n');
	    }
	} // end while (!timingCalFile) loop 
    }
    timingCalFile.close();
}
