/** \file VandleProcessor.hpp
 * \brief Class for handling Vandle Bars
 */

#ifndef __VANDLEPROCESSOR_HPP_
#define __VANDLEPROCESSOR_HPP_

#include "EventProcessor.hpp"
#include "TimingInformation.hpp"

class TFile;

class VandleProcessor : public EventProcessor, public TimingInformation
{
 public:
    VandleProcessor(TFile*); // no virtual c'tors
    VandleProcessor(const int VML_OFFSET, const int RANGE);
    VandleProcessor(const int RP_OFFSET, const int RANGE, int i);
    virtual void DeclarePlots(void);
    virtual bool Process(RawEvent &event);
    VMLMap vmlMap;
    DataRoot vandle;
	
    /// All processors with AddBranch() information
    vector<EventProcessor *> vecProcess;
    virtual bool Init(RawEvent& rawev);
    bool AddBranch(TTree *tree);
    UInt_t   multiplicity;
    UInt_t   vmllocation;

 protected:
    //define the maps
    BarMap barMap;
    TimingDataMap bigMap;
    TimingDataMap smallMap;
    TimingDataMap startMap;
    TimingDataMap tvandleMap;
 
 private:
    virtual bool RetrieveData(RawEvent &event);
    virtual double CorrectTOF(const double &TOF, const double &corRadius, const double &r) 
	{	
		return((r/corRadius)*TOF); //corRadius is the flight path of neutron from target to bar
	};

    virtual void AnalyzeData(RawEvent& rawev);
    virtual void BuildBars(const TimingDataMap &endMap, const std::string &type, BarMap &barMap);
    virtual void ClearMaps(void);
    virtual void CrossTalk(void);
    virtual void FillMap(const vector<ChanEvent*> &eventList, const std::string type, TimingDataMap &eventMap);
    virtual void Tvandle(void);
    virtual void WalkBetaVandle(const TimingInformation::TimingDataMap &beta, const TimingInformation::BarData &bar);
    virtual void FillRoot(const vmlData &tempData, UInt_t location);

    bool hasDecay;
    double decayTime;
    int counter;
	
    typedef std::pair<unsigned int, unsigned int> CrossTalkKey; 
    typedef std::map<CrossTalkKey, double> CrossTalkMap;
    std::map<CrossTalkKey, double> crossTalk;
	
}; //Class VandleProcessor
#endif // __VANDLEPROCESSOR_HPP_
