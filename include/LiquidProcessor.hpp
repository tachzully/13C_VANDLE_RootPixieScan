/** \file LiquidProcessor.hpp
 *
 * Processor for liquid scintillator detectors
 */

#ifndef __LIQUIDPROCESSOR_HPP_
#define __LIQUIDPROCESSOR_HPP_

#include "EventProcessor.hpp"
#include "TimingInformation.hpp"
#include "Trace.hpp"

#include "TFile.h"
#include "TTree.h"

class LiquidProcessor : public EventProcessor
{
public:
    TFile *topFile; //< Pointer to master root file (DetectorDriver)
    TTree *tree; //< ROOT tree where event branches are filled

    LiquidProcessor(TFile*); // no virtual c'tors
    virtual bool PreProcess(RawEvent &event);
    virtual bool Process(RawEvent &event);
    virtual void DeclarePlots(void);
private:
    TimingInformation timeInfo;
};
#endif // __LIQUIDPROCSSEOR_HPP_
