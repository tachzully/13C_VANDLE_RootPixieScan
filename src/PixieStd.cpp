/** \file PixieStd.cpp
 *
 * \brief pixie_std provides the interface between the HRIBF scan 
 * and the C++ analysis
 *
 * This provides the interface between the HRIBF scan and the C++ analysis
 * and as such is not a class in its own right.  In this file the data
 * received from scan is first reassembled into a pixie16 spill and then
 * channel objects are made.
 *
 * The main program.  Buffers are passed to hissub_() and channel information
 * is extracted in ReadBuffData(). All channels that fired are stored as a
 * vector of pointers which is sorted based on time and then events are built
 * with each event being sent to the detector driver for processing.
 *
 * \author S. Liddick 
 * \date 20 July 2007 
 *
 * <strong> Modified : </strong> <br>
 * S. Liddick - 2-5-08 - Added in diagnostic spectra including:
 *   runtime, channel time difference in an event, time difference between
 *   events, length of event, and length of buffer <br>
 * 
 * S. Liddick - 5-14-08 - At SP's request, error message and termination occur if a
 *   module number is encountered in the data stream that is not included in
 *   the map.txt file <br>
 *
 * David Miller - 5-5-10 - Significant changes throughout for conciseness, 
 *   optimization, and better error checking of incoming buffers
 */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <cstring>
#include <ctime>

#include <unistd.h>
#include <sys/times.h>

#include "pixie16app_defs.h"

#include "DetectorDriver.hpp"
#include "DetectorLibrary.hpp"
#include "DetectorSummary.hpp"
#include "ChanEvent.hpp"
#include "RawEvent.hpp"
#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "Plots.hpp"
#include "PlotsRegister.hpp"
#include "TreeCorrelator.hpp"

using namespace std;
using namespace dammIds::raw;
using pixie::word_t;

/**
 * Contains event information, the information is filled in ScanList() and is 
 * referenced in DetectorDriver.cpp, particularly in ProcessEvent()
 */
RawEvent rawev;

enum HistoPoints {BUFFER_START, BUFFER_END, EVENT_START = 10, EVENT_CONTINUE};

// Function forward declarations
void ScanList(vector<ChanEvent*> &eventList, RawEvent& rawev);
void RemoveList(vector<ChanEvent*> &eventList);
void HistoStats(unsigned int, double, double, HistoPoints);

#ifdef newreadout
/**
 * \brief Extract channel information from the raw parameter array ibuf
 */
void hissub_sec(unsigned int *ibuf[],unsigned int *nhw);
bool MakeModuleData(const word_t *data, unsigned long nWords); 
#endif

int ReadBuffData(word_t *lbuf, unsigned long *BufLen,
		 vector<ChanEvent *> &eventList);
void Pixie16Error(int errornum);

const string scanMode = "scan";

/** \fn extern "C" void hissub_(unsigned short *ibuf[],unsigned short *nhw) 
 * \brief interface between scan and C++
 *
 * In a typical experiment, Pixie16 reads data from all modules when one module
 * has hit the maximum number of events which is programmed during experimental
 * setup.  This spill of data is then broken into smaller chunks for
 * transmitting across the network.  The hissub_ function takes the chunks
 * and reconstructs the spill.
 *
 * Summarizing the terminology:
 *  - Spill  - a readout of all Pixie16 modules
 *  - Buffer - the data from a specific Pixie16 module
 *  - Chunk  - the packet transferred from Pixie16 to the acquisition
 *    
 * The hissub_ function is passed a pointer to an array with data (ibuf) and
 * the number of half words (nhw) contained in it.  This function is used with
 * the new Pixie16 readout (which is the default).  If the old Pixie16 readout
 * is used, the code should be recompiled without the newreadout flag in which
 * case this particular function is not used.  
*/

#ifdef newreadout

// THIS SHOULD NOT BE SET LARGER THAN 1,000,000
//  this defines the maximum amount of data that will be received in a spill
const unsigned int TOTALREAD = 1000000;

#if defined(REVD) || defined(REVF)
const unsigned int maxWords = EXTERNAL_FIFO_LENGTH; //Revision D
#else
const unsigned int maxWords = IO_BUFFER_LENGTH; // Revision A
#endif

// Catch the exit call from scanor and clean up c++ objects CRT
extern "C" void cleanup_()
{
    std::cout << "\nCalling C++ destructors...\n";
    DetectorDriver* driver = DetectorDriver::get();
    delete driver; // Call c++ destructors
}

extern "C" void hissub_(unsigned short *sbuf[],unsigned short *nhw)
{
    const unsigned int maxChunks = 200;

    static word_t totData[TOTALREAD];
    // keep track of the number of bad spills
    static unsigned int spillInvalidCount = 0, spillValidCount = 0;
    static bool firstTime = true;
    // might take a few entries into this function to get all the buffers in a spill
    static unsigned int bufInSpill = 0;    
    static unsigned int dataWords = 0;
    
    /*Assign ibuf variable to local variable for use in function */
    word_t *buf=(word_t*)sbuf;
    
    /* Initialize variables */
    unsigned long totWords=0;
    word_t nWords=buf[0] / 4;
    word_t totBuf=buf[1];
    word_t bufNum=buf[2];
    static unsigned int lastBuf = U_DELIMITER;

    // Check to make sure the number of buffers is not excessively large 
    if (totBuf > maxChunks) {
        cout << "LARGE TOTNUM = " << bufNum << endl;
        return;
    }

    /* Find a starting point in a file immediately following the 5-word buffer
         which indicates the end of a spill
     */
    if(bufNum != 0 && firstTime) {
        do {
            if (buf[totWords] == U_DELIMITER) {
                cout << "  -1 DELIMITER, " 
                    << buf[totWords] << buf[totWords + 1] << endl;
                return;
            }
            nWords = buf[totWords] / 4;
            totBuf = buf[totWords+1];
            bufNum = buf[totWords+2];
            totWords += nWords+1;
            cout << "SKIP " << bufNum << " of " << totBuf << endl;
        } while(nWords != 5);
    }
    firstTime = false;
    
    do {
	do {
	    /*Determine the number of words, total number of buffers, and
	      current buffer number at this point in the chunk.  
	      Note: the total number of buffers is repeated for each 
	      buffer in the chunk */
	    if (buf[totWords] == U_DELIMITER) return;

	    nWords = buf[totWords] / 4;
 	    bufNum = buf[totWords+2]; 
	    // read total number of buffers later after we check if the last spill was good
	    if (lastBuf != U_DELIMITER && bufNum != lastBuf + 1) {
#ifdef VERBOSE
		cout << "Buffer skipped, Last: " << lastBuf << " of " << totBuf 
		     << " buffers read -- Now: " << bufNum << endl;
#endif
		// if we are only missing the vsn 9999 terminator, reconstruct it
		if (lastBuf + 2 == totBuf && bufInSpill == totBuf - 1) {
#ifdef VERBOSE
		    cout << "  Reconstructing final buffer " << lastBuf + 1 << "." << endl;
#endif		   
		    totData[dataWords++] = 2;
		    totData[dataWords++] = 9999;
		    
		    MakeModuleData(totData, dataWords);
		    spillValidCount++;
		    bufInSpill = 0; dataWords = 0; lastBuf = -1;
		} else if (bufNum == 0) {
#ifdef VERBOSE		    
		    cout << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
			 << "  INCOMPLETE BUFFER " << spillInvalidCount 
			 << "\n  " << spillValidCount << " valid spills so far."
			 << " Starting fresh spill." << endl;
#endif		   
		    spillInvalidCount++;
		    // throw away previous collected data and start fresh
		    bufInSpill = 0; dataWords = 0; lastBuf = -1;
		}
	    } // check that the chunks are in order
	    // update the total chunks only after the sanity checks above
	    totBuf = buf[totWords+1];
	    if (totBuf > maxChunks) {
#ifdef VERBOSE
		cout << "EEEEE LOST DATA: Total buffers = " << totBuf 
		     <<  ", word count = " << nWords << endl;
#endif
		return;
	    }
	    if (bufNum > totBuf - 1) {
#ifdef VERBOSE
		cout << "EEEEEEE LOST DATA: Buffer number " << bufNum
		     << " of total buffers " << totBuf << endl;
#endif
		return;
	    }
	    lastBuf = bufNum;

	    /* Increment the number of buffers in a spill*/
	    bufInSpill++;
	    if(nWords == 0) {
#ifdef VERBOSE
		cout << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE NWORDS 0" << endl;
#endif
		return;
	    }
	    
	    /* Extract this buffer information into the TotData array*/
	    memcpy(&totData[dataWords], &buf[totWords+3], (nWords - 3) * sizeof(int));
	    dataWords += nWords - 3;
	    
	    // Increment location in file 
	    // one extra word to pass over "-1" delimiter signalling end of buffer
	    totWords += nWords+1;
	    if (bufNum == totBuf - 1 && nWords != 5) {
		cout << "Strange final buffer " << bufNum << " of " << totBuf
		     << " with " << nWords << " words" << endl;
	    }
	    if (nWords == 5 && bufNum != totBuf - 1) {
#ifdef VERBOSE
		cout << "Five word buffer " << bufNum << " of " << totBuf
		     << " WORDS: " 
		     << hex << buf[3] << " " << buf[4] << dec << endl;
#endif		
	    }
	} while(nWords != 5 || bufNum != totBuf - 1);
	/* reached the end of a spill when nwords = 5 and last chunk in spill */

	/* make sure we retrieved all the chunks of the spill */
	if (bufInSpill != totBuf) {
#ifdef VERBOSE	  
	    cout << "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE  INCOMPLETE BUFFER "
		 << spillInvalidCount
		 << "\n I/B [  " << bufInSpill << " of " << totBuf << " : pos " << totWords 
		 << "    " << spillValidCount << " total spills"
		 << "\n| " << hex << buf[0] << " " << buf[1] << "  " 
		 << buf[2] << " " << buf[3]
		 << "\n| " << dec << buf[totWords] << " " << buf[totWords+1] << "  "
		 << buf[totWords+2] << " " << buf[totWords+3] << endl;
#endif
	    spillInvalidCount++; 
	} else {
	    spillValidCount++;
	    MakeModuleData(totData, dataWords);	    
	} // else the number of buffers is complete
	dataWords = 0; bufInSpill = 0; lastBuf = -1; // reset the number of buffers recorded
    } while (totWords < nhw[0] / 4);
}

/** \brief inserts a delimiter in between individual module data and at end of 
 * buffer. Data is then passed to hissub_sec() for processing.
 */
bool MakeModuleData(const word_t *data, unsigned long nWords)
{
    const unsigned int maxVsn = 14; // no more than 14 pixie modules per crate

    unsigned int inWords = 0, outWords = 0;
    
    static word_t modData[TOTALREAD];
    // create a constant pointer to this data block for passing to hissub_sec
    static word_t* dataPtr = modData; 

    do {
	word_t lenRec = data[inWords];	
        word_t vsn    = data[inWords+1];
	/* Check sanity of record length and vsn*/
	if(lenRec > maxWords || (vsn > maxVsn && vsn != 9999 && vsn != 1000)) { 
#ifdef VERBOSE
	    cout << "SANITY CHECK FAILED: lenRec = " << lenRec
		 << ", vsn = " << vsn << ", inWords = " << inWords
		 << " of " << nWords << ", outWords = " << outWords << endl;
#endif
	    // exit(EXIT_FAILURE);
	    return false;  
	}
	
	/*Extract the data from TotData and place into ModData*/	      
	memcpy(&modData[outWords], &data[inWords], lenRec * sizeof(word_t)); 
	inWords  += lenRec;
	outWords += lenRec;
	
	modData[outWords++]=U_DELIMITER;
    } while (inWords < nWords);
	    
    modData[outWords++]=U_DELIMITER;
    modData[outWords++]=U_DELIMITER;
	    
    if(nWords > TOTALREAD || inWords > TOTALREAD || outWords > TOTALREAD ) {
	cout << "Values of nn - " << nWords << " nk - "<< inWords  
	     << " mm - " << outWords << " TOTALREAD - " << TOTALREAD << endl;
	Pixie16Error(2); 
	return false;
    }

    //! shouldn't this be 4 * outWords
    unsigned int nhw = 8 * outWords; // calculate the number of half short ints

    hissub_sec(&dataPtr, &nhw);

    return true;
}
#endif


/**
 * If the new Pixie16 readout is used (default), this routine processes the
 * reconstructed buffer.  Specifically, it retrieves channel information
 * and places the channel information into a list of channels that triggered in
 * this spill.  The list of channels is sorted according to the event time
 * assigned to each channel by Pixie16 and the sorted list is passed to
 * ScanList() for raw event creation. 
 *
 * If the old pixie readout is used then this function is
 * redefined as hissub_.
 */
#ifdef newreadout
void hissub_sec(word_t *ibuf[],unsigned int *nhw)
#else
extern "C" void hissub_(unsigned short *ibuf[],unsigned short *nhw)
#endif
{
    static float hz = sysconf(_SC_CLK_TCK); // get the number of clock ticks per second
    static clock_t clockBegin; // initialization time
    static struct tms tmsBegin;

    vector<ChanEvent*> eventList; // vector to hold the events

    /* Pointer to singleton DetectorLibrary class */
    DetectorLibrary* modChan = DetectorLibrary::get();
    /* Pointer to singleton DetectorDriver class */
    DetectorDriver* driver = DetectorDriver::get();

    // local version of ibuf pointer
    word_t *lbuf;

    int retval = 0; // return value from various functions
    
    unsigned long bufLen;
    
    /*
      Various event counters
    */
    unsigned long numEvents = 0;
    static int counter = 0; // the number of times this function is called
    static int evCount;     // the number of times data is passed to ScanList
    static unsigned int lastVsn; // the last vsn read from the data
    time_t theTime = 0;

    /*
      Assign the local variable lbuf to the variable ibuf which is passed into
      the routine.  The difference between the new and old pixie16 readouts is
      the type of the variable and source of the variable ibuf.

      In the new readout ibuf is from a C++ function and is of type unsigned int*
      In the old readout ibuf is from a Fortran function and is of type
      unsigned short*

      This results in two different assignment statements depending on 
      the readout.
    */
#ifdef newreadout
    lbuf=(word_t *)ibuf[0];
#else
    lbuf=(word_t *)ibuf; //old readout
#endif

    /* Initialize the scan program before the first event */
    if (counter==0) {
        /* Retrieve the current time for use later to determine the total
	 * running time of the analysis.
	 */
        clockBegin = times(&tmsBegin);

	cout << "First buffer at " << clockBegin << " sys time" << endl;
        /* After completion the descriptions of all channels are in the modChan
	 * vector, the DetectorDriver and rawevent have been initialized with the
	 * detectors that will be used in this analysis.
	 */
    cout << "Using event width " << pixie::eventInSeconds * 1e6 << " us" << endl
         << "                  " << pixie::eventWidth
         << " in pixie16 clock tics." << endl;

	modChan->PrintUsedDetectors(rawev);
    if (verbose::MAP_INIT)
        modChan->PrintMap();

	driver->Init(rawev);
    
	/* Make a last check to see that everything is in order for the driver 
	 * before processing data
	 */
	if ( !driver->SanityCheck() ) {
	    cout << "Detector driver did not pass sanity check!" << endl;
	    exit(EXIT_FAILURE);
	}

        lastVsn=-1; // set last vsn to -1 so we expect vsn 0 first 	

	cout << "Init done at " << times(&tmsBegin) << " sys time." << endl;
    }
    counter++;
 
    unsigned int nWords=0;  // buffer counter, reset only for new buffer
 
    // true if the buffer being analyzed is split across a spill from pixie
    bool multSpill;

    do {
	word_t vsn = U_DELIMITER;
	bool fullSpill=false; //true if spill had all vsn's
        multSpill=false;  //assume all buffers are not split between spills    

        /* while the current location in the buffer has not gone beyond the end
         * of the buffer (ignoring the last three delimiters, continue reading
	 */
        while (nWords < (nhw[0]/2 - 6)) {
            /*
              Retrieve the record length and the vsn number
            */
            word_t lenRec = lbuf[nWords];
            vsn = lbuf[nWords+1];
            
            /* If the record length is -1 (after end of spill), increment the
	       location in the buffer by two and start over with the while loop
	     */
            if (lenRec == U_DELIMITER) {
                nWords += 2;  // increment two whole words and try again
                continue;                         
            }
	    // Buffer with vsn 1000 was inserted with the time for superheavy exp't
	    if (vsn == clockVsn) {
	      memcpy(&theTime, &lbuf[nWords+2], sizeof(time_t));
	      nWords += lenRec;
	    }
            /*
              If the record length is 6, this is an empty channel.
	      Skip this vsn and continue with the next
            */
	    //! Revision specific, so move to ReadBuffData
            if (lenRec==6) {
                nWords += lenRec+1; // one additional word for delimiter
                lastVsn=vsn;
                continue;
            }
            
            /* If both the current vsn inspected is within an acceptable 
	       range, begin reading the buffer.
            */
            if ( vsn < modChan->GetPhysicalModules()  ) {
	        if ( lastVsn != U_DELIMITER) {
		    // the modules should be read out cyclically
		    if ( ((lastVsn+1) % modChan->GetPhysicalModules() ) != vsn ) {
#ifdef VERBOSE
			cout << " MISSING BUFFER " << vsn << "/" 
                 << modChan->GetPhysicalModules()
			     << " -- lastVsn = " << lastVsn << "  " 
			     << ", length = " << lenRec << endl;
#endif
                        RemoveList(eventList);
                        fullSpill=true;
                    }
                }
                /* Read the buffer.  After read, the vector eventList will 
		   contain pointers to all channels that fired in this buffer
                */

                retval= ReadBuffData(&lbuf[nWords], &bufLen, eventList);

                
                /* If the return value is less than the error code, 
		   reading the buffer failed for some reason.  
		   Print error message and reset variables if necessary
                */
                if ( retval <= readbuff::ERROR ) {
		    cout << " READOUT PROBLEM " << retval 
			 << " in event " << counter << endl;
                    if ( retval == readbuff::ERROR ) {
			cout << "  Remove list " << lastVsn << " " << vsn << endl;
                        RemoveList(eventList); 	                        
                    }
                    return;
                } else if ( retval == 0 ) {
		    // empty buffers are regular in Rev. D data
		    // cout << " EMPTY BUFFER" << endl;
		  nWords += lenRec + 1;
		  lastVsn = vsn;
		  continue;
                  //  return;
                } else if ( retval > 0 ) {		
		  /* increment the total number of events observed */
		  numEvents += retval;
                }
                /* Update the variables that are keeping track of what has been
		   analyzed and increment the location in the current buffer
                */
                lastVsn = vsn;
                nWords += lenRec+1; // one extra word for delimiter
            } else {
		// bail out if we have lost our place,		
		//   (bad vsn) and process events     
		if (vsn != 9999 && vsn != clockVsn) {
#ifdef VERBOSE	    
		    cout << "UNEXPECTED VSN " << vsn << endl;
#endif
		}
		break;
	    }
        } // while still have words
	if (nWords > nhw[0] / 2 - 6) {
	    cout << "This actually happens!" << endl;	    
	}
        
        /* If the vsn is 9999 this is the end of a spill, signal this buffer
	   for processing and determine if the buffer is split between spills.
        */
        if ( vsn == 9999 || vsn == clockVsn ) {
            fullSpill = true;
            nWords += 3;//skip it
            if (lbuf[nWords+1] != U_DELIMITER) {
		cout << "this actually happens!" << endl;
                multSpill = true;
            }
            lastVsn=U_DELIMITER;
        }
        
        /* if there are events to process, continue */
        if( numEvents>0 ) {
	    if (fullSpill) { 	  // if full spill process events
		// sort the vector of pointers eventlist according to time
		double lastTimestamp = (*(eventList.rbegin()))->GetTime();

		sort(eventList.begin(),eventList.end(),CompareTime);
		driver->CorrelateClock(lastTimestamp, theTime);

		/* once the vector of pointers eventlist is sorted based on time,
		   begin the event processing in ScanList()
		*/
		ScanList(eventList, rawev);

		/* once the eventlist has been scanned, remove it from memory
		   and reset the number of events to zero and update the event
		   counter
		*/

		evCount++;		
		/*
		  every once in a while (when evcount is a multiple of 1000)
		  print the time elapsed doing the analysis
		*/
		if(evCount % 1000 == 0 || evCount == 1) {
		    tms tmsNow;
		    clock_t clockNow = times(&tmsNow);

		    if (theTime != 0) {
			cout << " data read up to poll status time " << ctime(&theTime);
		    }
		    cout << "   buffer = " << evCount << ", user time = " 
			 << (tmsNow.tms_utime - tmsBegin.tms_utime) / hz
			 << ", system time = " 
			 << (tmsNow.tms_stime - tmsBegin.tms_stime) / hz
			 << ", real time = "
			 << (clockNow - clockBegin) / hz 
			 << ", ts = " << lastTimestamp << endl;
		}		
		RemoveList(eventList);
		numEvents=0;
	    } // end fullSpill 
	    else {
		cout << "Spill split between buffers" << endl;
		return; //! this tosses out all events read into the vector so far
	    }	    
        }  // end numEvents > 0
        else if (retval != readbuff::STATS) {
	    cout << "bad buffer, numEvents = " << numEvents << endl;
            return;
        }
        
    } while (multSpill); // end while loop over multiple spills
    return;      
}


/** Remove events in list from memory when no longer needed */
void RemoveList(vector<ChanEvent*> &eventList)
{
    /*
      using the iterator and starting from the beginning and going to 
      the end of eventlist, delete the actual objects
    */
    for(vector<ChanEvent*>::iterator it = eventList.begin();
	it != eventList.end(); it++) {
        delete *it;
    }
    
    // once the actual objects are deleted, clear the vector eventList
    eventList.clear();   
}

/** \brief event by event analysis
 * 
 * ScanList() operates on the time sorted list of all channels that triggered in
 * a given spill.  Starting from the begining of the list and continuing to the
 * end, an individual channel event time is compared with the previous channel
 * event time to determine if they occur within a time period defined by the
 * diff_t variable (time is in units of 10 ns).  Depending on the answer,
 * different actions are performed:
 *   - yes - the two channels are grouped together as belonging to the same event
 *   and the current channel is added to the list of channels for the rawevent 
 *   - no - the previous rawevent is sent for processing and once finished, the
 *   rawevent is zeroed and the current channel placed inside it.
 */

void ScanList(vector<ChanEvent*> &eventList, RawEvent& rawev) 
{
    unsigned long chanTime, eventTime;

    DetectorLibrary* modChan = DetectorLibrary::get();
    DetectorDriver* driver = DetectorDriver::get();

    // local variable for the detectors used in a given event
    set<string> usedDetectors;
    
    vector<ChanEvent*>::iterator iEvent = eventList.begin();

    // local variables for the times of the current event, previous
    // event and time difference between the two
    double diffTime = 0;
    
    //set last_t to the time of the first event
    double lastTime = (*iEvent)->GetTime();
    double currTime = lastTime;
    unsigned int id = (*iEvent)->GetID();

    HistoStats(id, diffTime, lastTime, BUFFER_START);

    //loop over the list of channels that fired in this buffer
    for(; iEvent != eventList.end(); iEvent++) { 
        id = (*iEvent)->GetID();
        if (id == U_DELIMITER) {
            cout << "pattern 0 ignore" << endl;
            continue;
        }
        if ((*modChan)[id].GetType() == "ignore") {
            continue;
        }

        // this is a channel we're interested in
        chanTime  = (*iEvent)->GetTrigTime(); 
        eventTime = (*iEvent)->GetEventTimeLo();

        /* retrieve the current event time and determine the time difference 
        between the current and previous events. 
        */
        currTime = (*iEvent)->GetTime();
        diffTime = currTime - lastTime;

        /* if the time difference between the current and previous event is 
        larger than the event width, finalize the current event, otherwise
        treat this as part of the current event
        */
        if ( diffTime > pixie::eventWidth ) {
            if(rawev.Size() > 0) {
            /* detector driver accesses rawevent externally in order to
            have access to proper detector_summaries
            */
                driver->ProcessEvent(scanMode, rawev);
            }
    
            //after processing zero the rawevent variable
            rawev.Zero(usedDetectors);
            usedDetectors.clear();	    

            // Now clear all places in correlator (if resetable type)
            for (map<string, Place*>::iterator it =
                    TreeCorrelator::get()->places_.begin();
                it != TreeCorrelator::get()->places_.end(); ++it)
                if ((*it).second->resetable())
                    (*it).second->reset();

            HistoStats(id, diffTime, currTime, EVENT_START);
        } else HistoStats(id, diffTime, currTime, EVENT_CONTINUE);

        unsigned long dtimebin = 2000 + eventTime - chanTime;
        if (dtimebin < 0 || dtimebin > (unsigned)(SE)) {
            cout << "strange dtime for id " << id << ":" << dtimebin << endl;
        }
        driver->plot(D_TIME + id, dtimebin);

        usedDetectors.insert((*modChan)[id].GetType());
        rawev.AddChan(*iEvent);

        // update the time of the last event
        lastTime = currTime; 
    } //end loop over event list

    //process the last event in the buffer
    if (rawev.Size() > 0) {
        string mode;
        HistoStats(id, diffTime, currTime, BUFFER_END);

        driver->ProcessEvent(scanMode, rawev);
        rawev.Zero(usedDetectors);
    }
}

/**
 * At various points in the processing of data in ScanList(), HistoStats() is
 * called to increment some low level pixie16 informational and diagnostic
 * spectra.  The list of spectra filled includes runtime in second and
 * milliseconds, the deadtime, time between events, and time width of an event.
 */
void HistoStats(unsigned int id, double diff, double clock, HistoPoints event)
{
    static const int specNoBins = SE;

    static double start, stop;
    static int count;
    static double firstTime = 0.;
    static double bufStart;

    double runTimeSecs   = (clock - firstTime) * pixie::clockInSeconds;
    int    rowNumSecs    = int(runTimeSecs / specNoBins);
    double remainNumSecs = runTimeSecs - rowNumSecs * specNoBins;

    double runTimeMsecs   = runTimeSecs * 1000;
    int    rowNumMsecs    = int(runTimeMsecs / specNoBins);
    double remainNumMsecs = runTimeMsecs - rowNumMsecs * specNoBins;

    static double bufEnd = 0, bufLength = 0;
    // static double deadTime = 0 // not used
    DetectorDriver* driver = DetectorDriver::get();

    if (firstTime > clock) {
        cout << "Backwards clock jump detected: prior start " << firstTime
            << ", now " << clock << endl;
        // detect a backwards clock jump which occurs when some of the
        //   last buffers of a previous run sneak into the beginning of the 
        //   next run, elapsed time of last buffers is usually small but 
        //   just in case make some room for it
        double elapsed = stop - firstTime;
        // make an artificial 10 second gap by 
        //   resetting the first time accordingly
        firstTime = clock - 10 / pixie::clockInSeconds - elapsed;
        cout << elapsed*pixie::clockInSeconds << " prior seconds elapsed "
            << ", resetting first time to " << firstTime << endl;	
    }

    switch (event) {
        case BUFFER_START:
            bufStart = clock;
            if(firstTime == 0.) {
                firstTime = clock;
            } else if (bufLength != 0.){
                //plot time between buffers as a function of time - dead time spectrum	    
                // deadTime += (clock - bufEnd)*pixie::clockInSeconds;
                // plot(DD_DEAD_TIME_CUMUL,remainNumSecs,rownum,int(deadTime/runTimeSecs));	    	    
                driver->plot(dammIds::raw::DD_BUFFER_START_TIME, remainNumSecs,rowNumSecs, (clock-bufEnd)/bufLength*1000.);	    
            }
            break;
        case BUFFER_END:
            driver->plot(D_BUFFER_END_TIME, (stop - bufStart) * pixie::clockInSeconds * 1000);
            bufEnd = clock;
            bufLength = clock - bufStart;
        case EVENT_START:
            driver->plot(D_EVENT_LENGTH, stop - start);
            driver->plot(D_EVENT_GAP, diff);
            driver->plot(D_EVENT_MULTIPLICITY, count);
            
            start = stop = clock; // reset the counters      
            count = 1;
            break;
        case EVENT_CONTINUE:
            count++;
            if(diff > 0.) {
                driver->plot(D_SUBEVENT_GAP, diff + 100);
            }
            stop = clock;
            break;
        default:
            cout << "Unexpected type " << event << " given to HistoStats" << endl;
    }

    //fill these spectra on all events, id plots and runtime.
    // Exclude event type 0/1 since it will also appear as an
    // event type 11
    if ( event != BUFFER_START && event != BUFFER_END ){      
        driver->plot(DD_RUNTIME_SEC, remainNumSecs, rowNumSecs);
        driver->plot(DD_RUNTIME_MSEC, remainNumMsecs, rowNumMsecs);
        //fill scalar spectrum (per second) 
        driver->plot(D_HIT_SPECTRUM, id);
        driver->plot(D_SCALAR + id, runTimeSecs);
    }
}


/** \brief pixie16 scan error handling.
 *
 * Print out an error message and terminate program depending on value
 * of errorNum. 
 */
void Pixie16Error(int errorNum)
{
  //print error message depending on errornum
  switch (errorNum) {
      case 1:
	  cout << endl;
	  cout << " **************  SCAN ERROR  ****************** " << endl;
	  cout << "There appears to be more modules in the data " << endl;
	  cout << "stream than are present in the map.txt file. " << endl;
	  cout << "Please verify that the map.txt file is correct " << endl;
	  cout << "This is a fatal error, program terminating" << endl;
	  exit(EXIT_FAILURE);
      case 2:
	  cout << endl;
	  cout << "***************  SCAN ERROR  ******************* "<<endl;
	  cout << "One of the variables named nn, nk, or mm" << endl;
	  cout << "have exceeded the value of TOTALREAD. The value of" << endl;
	  cout << "TOTALREAD MUST NEVER exceed 1000000 for correct " << endl;
	  cout << "opertation of code between 32-bit and 64-bit architecture " << endl;
	  cout << "Either these variables have not been zeroed correctly or" << endl;
	  cout << "the poll program controlling pixie16 is trying to send too " << endl;
	  cout << "much data at once" << endl;
	  cout << "This is a fatal error, program terminating " << endl;
	  exit(EXIT_FAILURE);
  }
}
