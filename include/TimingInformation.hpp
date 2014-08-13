/** \file TimingInformation.hpp
 * \brief File holding structures and info for Timing Analysis
 */

#ifndef __TIMINGINFORMATION_HPP_
#define __TIMINGINFORMATION_HPP_

#include <map>
//#include "ChanEvent.hpp"

#ifdef useroot
#include "Rtypes.h"
#endif

// Forward declarations for hell of circular dependencies
// see Trace.hpp
class Trace;
// see ChanEvent.hpp
class ChanEvent;

class TimingInformation
{
 public:
    struct TimingCal {
		//constants from reading in timingCal.txt file
		double x;
		double y;
		double z;
		double r;
		double barPosTheta;
		double barPosPhi;
		double orientTheta;
		double orientPhi;
		double lrtOffset;
		double tofOffset0;
		double tofOffset1;
    };
    
    struct TimingData {
		TimingData(void);
		TimingData(ChanEvent *chan);
    	
		const Trace &trace;
	
		bool dataValid;
	
		double aveBaseline;
		double discrimination;
		double highResTime;
		double maxpos;
		double maxval;
		double phase;
		double snr;
		double stdDevBaseline;
		double tqdc;
		double walk;
		double walkCorTime;

		int numAboveThresh;
    };
    
#ifdef useroot
    struct DataRoot { //this struct is going into root tree
		
		static const size_t maxMultiplicity = 9;
		DataRoot(void);


		/*Double_t tof[maxMultiplicity];
		Double_t lqdc[maxMultiplicity];
		Double_t rqdc[maxMultiplicity];
		Double_t tsLow[maxMultiplicity];
		Double_t tsHigh[maxMultiplicity];
		Double_t lMaxVal[maxMultiplicity];
		Double_t rMaxVal[maxMultiplicity];
		Double_t qdc[maxMultiplicity];
		Double_t energy[maxMultiplicity];
		UInt_t location[maxMultiplicity]; //--- from key*/
		Double_t tof;
		Double_t lqdc;
		Double_t rqdc;
		Double_t tsLow;
		Double_t tsHigh;
		Double_t lMaxVal;
		Double_t rMaxVal;
		Double_t qdc;
		Double_t energy;
	
		Double_t ejectAngle;
		Double_t recoilEnergy;
		Double_t recoilAngle;
		Double_t exciteEnergy;
		Double_t flightPath;
		Double_t xflightPath;
		Double_t yflightPath;
		Double_t zflightPath;
			
		UInt_t   multiplicity;
		UInt_t   dummy;
		UInt_t location; //--- from key
				
    };
#endif

    struct BarData{
	
		BarData(const TimingData& Right, const TimingData& Left, const TimingCal &cal, const std::string &type);

		bool BarEventCheck(const double &timeDiff, const std::string &type);
		bool event;
		
		double CalcFlightPath(double &timeDiff, const TimingCal &cal, const std::string &type, 
					double &xflightPath, double &yflightPath, double &zflightPath);
					  
		
		double flightPath;
		double zflightPath;
		double yflightPath;
		double xflightPath;
		double ejectAngle;
		double recoilAngle;
		double recoilEnergy;
		double exciteEnergy;
		
		double lMaxVal;
		double lqdc;
		double lTime;
		double qdc;
		double qdcPos;
		double rMaxVal;
		double rqdc;
		double rTime;
		double theta;
		double timeAve;
		double timeDiff;
		double walkCorTimeDiff;
		double walkCorTimeAve;
		
		std::map<unsigned int, double> timeOfFlight;
		std::map<unsigned int, double> corTimeOfFlight;
		std::map<unsigned int, double> energy;
		
    };

    struct vmlData { //needs to have same quantities as DataRoot, fill this from already processed data
    
		vmlData(BarData bar, double corTOF, double enrgy, double tLow, double tHigh, double recoilE);	

		//bool dataValid;
			
		double tof;
		double lqdc;
		double rqdc;
		double tsLow;
		double tsHigh;	
		double lMaxVal;
		double rMaxVal;
		double qdc;
		double energy;
		double ejectAngle;
		double recoilEnergy;
		double recoilAngle;
		double exciteEnergy;
		double flightPath;
		double xflightPath;
		double yflightPath;
		double zflightPath;

		int multiplicity;
		int dummy;
		int location;
    };

    //define types for the keys and maps
    typedef std::pair<unsigned int, std::string> IdentKey;
    typedef std::map<IdentKey, struct TimingData> TimingDataMap;
    typedef std::map<IdentKey, struct BarData> BarMap;
    typedef std::map<unsigned int, struct vmlData> VMLMap; //--- key could be just an int, the bar location
    typedef std::map<IdentKey, struct TimingCal> TimingCalMap;
    typedef std::map<unsigned int, double> TimeOfFlightMap;

    double CalcEnergy(const double &timeOfFlight, const double &r);
    double CalcRecoilEnergy(const double &energy, const double &flightPath, double &zflightPath, double &ejectAngle, double &recoilAngle, double &exciteEnergy);
    
    static double GetConstant(const std::string &value);
    static TimingCal GetTimingCal(const IdentKey &identity);
    static void ReadTimingCalibration(void);
    static void ReadTimingConstants(void);
    
 private:
    static const double qdcCompression = 4.0;
    static std::map<std::string, double> constantsMap;
    static TimingCalMap calibrationMap;

}; // class TimingInformation
#endif //__TIMINGINFORMATION_HPP_
