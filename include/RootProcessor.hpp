/** \file RootProcessor.hpp 
 * \brief Processor to dump data from events into a root tree
 *
 * This loops over other event processor to fill appropriate branches
 */

#ifndef useroot
#error USEROOT must be defined to use RootProcessor
#endif

#ifndef __ROOTPROCESSOR_HPP_
#define __ROOTPROCESSOR_HPP_

#include <vector>

#include "EventProcessor.hpp"

// forward declaration
class TFile;

using std::vector;

class RootProcessor : public VandleProcessor
{
 private:
    TFile *topFile; //< Pointer to master root file (DetectorDriver)
    TTree *tree; //< ROOT tree where event branches are filled

    /// All processors with AddBranch() information
    vector<EventProcessor *> vecProcess;
 public:
    RootProcessor(TFile*);
    virtual bool Init(RawEvent& rawev);
    virtual bool Process(RawEvent &event);
    //virtual bool HasEvent();
};

#endif // __ROOTPROCESSOR_HPP_
