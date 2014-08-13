/** \file RootProcessor.cpp
 * \brief Implemenation of class to dump event info to a root tree
 * \author David Miller
 * \date Jan 2010
 */

#ifdef useroot

#include <algorithm>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "DetectorDriver.hpp"
#include "RootProcessor.hpp"

using std::cout;
using std::endl;

using namespace dammIds::vandle;

namespace dammIds {
	namespace vandle {
    }
}

/** Open a file for tree output */
//RootProcessor::RootProcessor(TFile* topFile_) : VandleProcessor(RP_OFFSET, RANGE, topFile_, "VandleROOT")
RootProcessor::RootProcessor(TFile* topFile_) : VandleProcessor(topFile_)
{
    std::cout << " RootProcessor: Initializing\n";
    name = "RootProcessor";
}

/** Add branches to the tree from the event processors in the driver */
bool RootProcessor::Init(RawEvent& rawev)
{
    DetectorDriver* driver = DetectorDriver::get();
    if(!topFile || !tree){
        cout << " RootProcessor: Failed to create ROOT objects\n";
        return false;
    }

    const vector<EventProcessor *>& drvProcess = driver->GetProcessors();

    for (vector<EventProcessor *>::const_iterator it = drvProcess.begin();
        it != drvProcess.end(); it++) {
        if ((*it)->AddBranch(tree)) {
        vecProcess.push_back(*it);
        set_union( (*it)->GetTypes().begin(), (*it)->GetTypes().end(),
            associatedTypes.begin(), associatedTypes.end(),
            inserter(associatedTypes, associatedTypes.begin()) );
        }	  
    } 
    return EventProcessor::Init(rawev);
}

/** Fill the tree for each event, saving to file occasionally */
bool RootProcessor::Process(RawEvent &event)
{
    if (!EventProcessor::Process(event))
	return false;

    //for (vector<EventProcessor *>::iterator it = vecProcess.begin();
	 //it != vecProcess.end(); it++) {
	//(*it)->FillBranch();	
    //}

    tree->Fill();
    if(tree->GetEntries() % 1000 == 0){ tree->AutoSave(); }
    if(tree->GetEntries() % 10000 == 0){ std::cout << " " << tree->GetEntries() << " entries\n"; }

    EndProcess();
    return true;
}

/** See if the detectors of interest have any events */
/*bool RootProcessor::HasEvent()
{
    return true;
}*/

#endif // useroot
