/** \file EventProcessor.hpp
 * \brief Generic event processor
 */

#ifndef __EVENTPROCESSOR_HPP_
#define __EVENTPROCESSOR_HPP_

#include <map>
#include <set>
#include <string>

#include <sys/times.h>
#include "Plots.hpp"
#include "TreeCorrelator.hpp"

#include "TFile.h"
#include "TTree.h"

// forward declarations
class DetectorSummary;
class RawEvent;

class EventProcessor {
 private:
    // things associated with timing
    tms tmsBegin;
    double userTime;
    double systemTime;
    double clocksPerSecond;

 protected:
    // define the associated detector types and only initialize if present
    std::string name;
    std::set<std::string> associatedTypes; //--- set of type string
    bool initDone;
    bool didProcess;
    TFile *topFile; //< Pointer to master root file (DetectorDriver)
    TTree *tree; //< ROOT tree where event branches are filled

    // map of associated detector summary
    std::map<std::string, const DetectorSummary *> sumMap;

    // Plots class for given Processor, takes care of declaration
    // and plotting within boundries allowed by PlotsRegistry
    Plots histo;

    virtual void plot(int dammId, double val1, double val2 = -1, double val3 = -1, const char* name="h") {
        histo.Plot(dammId, val1, val2, val3, name);
    }
    virtual void DeclareHistogram1D(int dammId, int xSize, const char* title) {
        histo.DeclareHistogram1D(dammId, xSize, title);
    }
    virtual void DeclareHistogram2D(int dammId, int xSize, int ySize, const char* title) {
        histo.DeclareHistogram2D(dammId, xSize, ySize, title);
    }

 public:
    EventProcessor();
    EventProcessor(int offset, int range);
    EventProcessor(TFile*, std::string);
    EventProcessor(int offset, int range, TFile*, std::string);
    virtual ~EventProcessor();
    virtual void Close();
    void SetTopFile(TFile* topFile_){ topFile = topFile_; }

    // declare associated damm plots (called by drrsub_)
    virtual void DeclarePlots(void);
    virtual const std::set<std::string>& GetTypes(void) const {
      return associatedTypes; 
    }
    virtual bool DidProcess(void) const {
      return didProcess;
    }
    // return true on success
    virtual bool HasEvent(void) const;
    virtual bool Init(RawEvent& event);
    virtual bool PreProcess(RawEvent &event);   
    virtual bool Process(RawEvent &event);   
    void EndProcess(void); // stop the process timer
    std::string GetName(void) const {
      return name;
    }
#ifdef useroot
    virtual bool AddBranch(TTree *tree);
    virtual void FillBranch(void);
#endif
};

#endif // __EVENTPROCESSOR_HPP_
