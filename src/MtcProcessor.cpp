/** \file MtcProcessor.cpp
 * \brief handling of mtc events
 *
 * derived from timeclass.cpp
 * doesn't handle old style NSCL correlations
 *
 * Start subtype corresponds to leading edge of tape move signal
 * Stop subtype corresponds to trailing edge of tape move signal
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "RawEvent.hpp"
#include "MtcProcessor.hpp"

using namespace std;
using namespace dammIds::mtc;

namespace dammIds {
    namespace mtc {
        const int D_TDIFF_BEAM_START   = 0;
        const int D_TDIFF_BEAM_STOP    = 1;
        const int D_TDIFF_MOVE_START   = 2;
        const int D_TDIFF_MOVE_STOP    = 3;
        const int D_MOVETIME           = 4;
        const int D_BEAMTIME           = 5;
        const int D_COUNTER            = 6;
        const int DD_TIME__DET_MTCEVENTS = 10;

        const int MOVE_START_BIN = 1;
        const int MOVE_STOP_BIN = 3;
        const int BEAM_START_BIN = 5;
        const int BEAM_STOP_BIN = 7;
    }
} // mtc namespace


MtcProcessor::MtcProcessor(void) : EventProcessor(OFFSET, RANGE), 
				   lastStartTime(NAN), lastStopTime(NAN)
{
    name = "mtc";

    associatedTypes.insert("timeclass"); // old detector type
    associatedTypes.insert("mtc");
}

void MtcProcessor::DeclarePlots(void)
{
    using namespace dammIds::mtc;
    
    const int counterBins = S3;
    const int timeBins = SA;

    DeclareHistogram1D(D_TDIFF_BEAM_START, timeBins, "Time diff btwn beam starts, 10 ms/bin");
    DeclareHistogram1D(D_TDIFF_BEAM_STOP, timeBins, "Time diff btwn beam stops, 10 ms/bin");
    DeclareHistogram1D(D_TDIFF_MOVE_START, timeBins, "Time diff btwn move starts, 10 ms/bin");
    DeclareHistogram1D(D_TDIFF_MOVE_STOP, timeBins, "Time diff btwn move stops, 10 ms/bin");
    DeclareHistogram1D(D_MOVETIME, timeBins, "Move time, 10 ms/bin");
    DeclareHistogram1D(D_BEAMTIME, timeBins, "Beam on time, 10 ms/bin");
    // Counter of events; see dammIds::mtc namespace for bin definition
    DeclareHistogram1D(D_COUNTER, counterBins, "MTC and beam counter");

    DeclareHistogram2D(DD_TIME__DET_MTCEVENTS, SF, S2, "MTC and beam events");
}

bool MtcProcessor::PreProcess(RawEvent &event)
{
    if (!EventProcessor::PreProcess(event))
        return false;

    const static DetectorSummary *mtcSummary = NULL;

    // plot with 10 ms bins
    const double mtcPlotResolution = 10e-3 / pixie::clockInSeconds;
    // for 2d plot of events 100ms / bin

    if (mtcSummary == NULL) {
        if ( sumMap.count("mtc") )
            mtcSummary = sumMap["mtc"];
        else if ( sumMap.count("timeclass") ) 
            mtcSummary = sumMap["timeclass"];
    }

    static const vector<ChanEvent*> &mtcEvents = mtcSummary->GetList();

    for (vector<ChanEvent*>::const_iterator it = mtcEvents.begin();
	 it != mtcEvents.end(); it++) {
        string subtype = (*it)->GetChanID().GetSubtype();
        double time   = (*it)->GetTime();	
        // Time of the first event
        static double t0 = time;
        string place = (*it)->GetChanID().GetPlaceName();

        const double eventsResolution = 100e-3 / pixie::clockInSeconds;
        const unsigned MTC_START = 0;
        const unsigned MTC_STOP = 1;
        const unsigned BEAM_START = 2;
        const unsigned BEAM_STOP = 3;
        double time_x = int((time - t0) / eventsResolution);

        if(place == "mtc_start_0") {

            double dt_start = time - 
                     TreeCorrelator::get()->place(place)->secondlast().time;
            TreeCorrelator::get()->place("TapeMove")->activate(time);
            TreeCorrelator::get()->place("Cycle")->deactivate(time);

            plot(D_TDIFF_MOVE_START, dt_start / mtcPlotResolution);
            plot(D_COUNTER, MOVE_START_BIN);
            plot(DD_TIME__DET_MTCEVENTS, time_x, MTC_START);

        } else if (place == "mtc_stop_0") {

            double dt_stop = time - 
                     TreeCorrelator::get()->place(place)->secondlast().time;
            double dt_move = time - 
                     TreeCorrelator::get()->place("mtc_start_0")->last().time;
            TreeCorrelator::get()->place("TapeMove")->deactivate(time);

            plot(D_TDIFF_MOVE_STOP, dt_stop / mtcPlotResolution);
            plot(D_MOVETIME, dt_move / mtcPlotResolution);
            plot(D_COUNTER, MOVE_STOP_BIN);
            plot(DD_TIME__DET_MTCEVENTS, time_x, MTC_STOP);

        } else if (place == "mtc_beam_start_0") {

            double dt_start = time -
                      TreeCorrelator::get()->place(place)->secondlast().time;
            TreeCorrelator::get()->place("Beam")->activate(time);
            TreeCorrelator::get()->place("Cycle")->activate(time);

            plot(D_TDIFF_BEAM_START, dt_start / mtcPlotResolution);
            plot(D_COUNTER, BEAM_START_BIN);
            plot(DD_TIME__DET_MTCEVENTS, time_x, BEAM_START);

        } else if (place == "mtc_beam_stop_0") {

            double dt_stop = time - 
                      TreeCorrelator::get()->place(place)->secondlast().time;
            double dt_beam = time - 
                      TreeCorrelator::get()->place("mtc_beam_start_0")->last().time;
            TreeCorrelator::get()->place("Beam")->deactivate(time);

            plot(D_TDIFF_BEAM_STOP, dt_stop / mtcPlotResolution);
            plot(D_BEAMTIME, dt_beam / mtcPlotResolution);
            plot(D_COUNTER, BEAM_STOP_BIN);
            plot(DD_TIME__DET_MTCEVENTS, time_x, BEAM_STOP);

        }
    }
    return true;
}

bool MtcProcessor::Process(RawEvent &event)
{
    if (!EventProcessor::Process(event))
        return false;
    using namespace dammIds::mtc;
    EndProcess(); // update processing time
    return true;
}
