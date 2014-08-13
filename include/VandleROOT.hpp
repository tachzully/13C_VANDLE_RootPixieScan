/** \file VandleROOT.hpp
 * \brief Root functionality for Vandle Processor
 */

#ifndef __VANDLEROOT_HPP_
#define __VANDLEROOT_HPP_

class TTree;

#include "VandleProcessor.hpp"

class VandleROOT : public VandleProcessor
{
 public: 
    VandleROOT(TFile*); //--- new constructor
    bool AddBranch(TTree *tree);
    bool Process(RawEvent &event);
    void FillBranch();
    void vmlMapCpy(const VMLMap vmlMap1);
    const VMLMap vmlMap2;
    TFile *topFile; //< Pointer to master root file (DetectorDriver)
    TTree *tree; //< ROOT tree where event branches are filled

 private:
    DataRoot vandle;

    bool isSmall;
    bool isBig;

    virtual void FillRoot(const VMLMap & endMap, const string &barType);
}; // class VandleROOT
#endif // __VANDLEROOT_HPP_
