/** \file VandleROOT.cpp
 *\brief Outputs VANDLE info to root - NOT UNDER DEVELOPMENT
 *
 *Writes certain information from the ScintProcessor into a ROOT tree.
 *This functionality is currently not supported, and may well give garbage.
 *
 *\author S. V. Paulauskas  
 *\date 09 May 2011
 */
#include <string>

#include "VandleROOT.hpp"
#include "DammPlotIds.hpp"
#include "VandleProcessor.hpp"

#include <TTree.h>

using std::string;

using namespace dammIds::vandle;

namespace dammIds {
	namespace vandle {
    }
}

//********** VandleROOT Constructor ********
//VandleROOT::VandleROOT(TFile *topFile_) : EventProcessor (VML_OFFSET, RANGE, topFile_, "VandleROOT")
VandleROOT::VandleROOT(TFile *topFile_) : VandleProcessor(topFile_)
{ 
    std::cout << " VandleROOT: Initializing\n";
    associatedTypes.insert("vandleSmall");
    associatedTypes.insert("vandleBig");
    associatedTypes.insert("tvandle");
}

//********** AddBranch **********
bool VandleROOT::AddBranch(TTree *tree)
{
    if (tree) {
        std::cout << " VandleROOT: Adding branches\n";
	string branchDef = "multiplicity/i:dummy/i:aveBaseline[10]/D:discrimination[10]/D:highResTime[10]/D:maxpos[10]/D:maxval[10]/D:phase[10]/D:stdDevBaseline[10]/D:tqdc[10]/D:location[10]/i";
	TBranch *vandleBranch  = tree->Branch("Vandle", &vandle, branchDef.c_str());
	/*TBranch *vandleBranch = tree->Branch("VandleSmallRight", &smallRight, branchDef.c_str()); 
	TBranch *vandleBranch1 = tree->Branch("VandleSmallLeft", &smallLeft, branchDef.c_str());
	TBranch *vandleBranch2  = tree->Branch("VandleBigRight", &bigRight, branchDef.c_str()); 
	TBranch *vandleBranch3 = tree->Branch("VandleBigLeft", &bigLeft, branchDef.c_str());*/
	
	return ( vandleBranch != NULL );
    } 
    return false;
}

bool VandleROOT::Process(RawEvent &event)
{
    if (!EventProcessor::Process(event))
	return false;

    FillBranch();	
    EndProcess();
    return true;
}

//********** vmlMapCpy **********
void VandleROOT::vmlMapCpy(const VMLMap vmlMap1)
{
    cout <<"VR1: "<<vmlMap1.size()<<endl;
//    vmlMap2 = vmlMap1;
}

//********** FillBranch **********
void VandleROOT::FillBranch()
{
    cout <<"VR2: "<<vmlMap.size()<<endl;
    FillRoot(vmlMap, "vml");
    /*FillRoot(bigMap, "big");
    FillRoot(smallMap, "small");
    
    if (!HasEvent())
	smallRight = smallLeft = bigRight = bigLeft = DataRoot();*/
}

//********** FillRoot **********
/* Perhaps the VMLMap does not possess many startBars from many different events. The event vectors
	might just record the pulses from one event. If so, those pulses are then put into a TDM and
	then correlated into a barEvent in the barMap. The effect would be to have a single barMap
	for one event, with multiple bars (paired pulses). If so, the original code would make much
	more sense, and the current code could easily manage it.
*/
void VandleROOT::FillRoot(const VMLMap &endMap, const string &barType)
{
    vandle = DataRoot(); //--- initiallizes data root variable arrays to 0, check to see if need new
    DataRoot *data; //--- creates DataRoot structure called data
    cout <<"VR2: "<<endMap.size()<<endl;
    for(VMLMap::const_iterator itTempA = endMap.begin();  //---loops over vml map, creating root structure
	itTempA != endMap.end(); itTempA++) {  

	data = & vandle; //--- sets data values to 0

	const vmlData &tempData = (*itTempA).second; //--- struct that absorbs values from current iteration

	/*data->tof[data->multiplicity] = tempData.tof;
	data->lqdc[data->multiplicity] = tempData.lqdc;
	data->rqdc[data->multiplicity] = tempData.rqdc;
	data->tsLow[data->multiplicity] = tempData.tsLow;
	data->tsHigh[data->multiplicity] = tempData.tsHigh;
	data->lMaxVal[data->multiplicity] = tempData.lMaxVal;
	data->rMaxVal[data->multiplicity] = tempData.rMaxVal;
	data->qdc[data->multiplicity] = tempData.qdc;
	data->energy[data->multiplicity] = tempData.energy;
	data->location[data->multiplicity] = (*itTempA).first;*/
	data->tof = tempData.tof;
	data->lqdc = tempData.lqdc;
	data->rqdc = tempData.rqdc;
	data->tsLow = tempData.tsLow;
	data->tsHigh = tempData.tsHigh;
	data->lMaxVal = tempData.lMaxVal;
	data->rMaxVal = tempData.rMaxVal;
	data->qdc = tempData.qdc;
	data->energy = tempData.energy;
	data->recoilEnergy = tempData.recoilEnergy;
    	data->recoilAngle = tempData.recoilAngle;
   	data->ejectAngle = tempData.ejectAngle;
    	data->exciteEnergy = tempData.exciteEnergy;

	data->location = (*itTempA).first;
	data->multiplicity++;
        cout <<"tempData.qdc = " << tempData.qdc << endl;
    }    
} 
