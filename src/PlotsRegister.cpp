/** \file PlotsRegister.cpp
 * \brief Declaration of singleton register
 */

#include "PlotsRegister.hpp"
#include <iostream>
#include <cstdlib>

using namespace std;

PlotsRegister* PlotsRegister::instance = NULL;

/** Instance is created upon first call */
PlotsRegister* PlotsRegister::get() {
    if (!instance) {
        instance = new PlotsRegister();
    }
    return instance;
}

bool PlotsRegister::CheckRange (int min, int max) const
{
    bool exists = false;
    
    unsigned sz = reg.size();
    for (unsigned id = 0; id < sz; ++id) {
        if ( (min < reg[id].first && max < reg[id].first) ||
             (min > reg[id].second && max > reg[id].second) )
            continue;
        else {
            exists = true;
            break;
        }
    }
    return exists;
}

bool PlotsRegister::Add (int offset, int range) //--- offset and range found in DammPlotIds
{
    // Special case: empty Plots list
    if (offset == 0 && range == 0)
        return true;
    
    int min = offset;
    int max = offset + range - 1;
    
    if (max < min) {
        cerr << "  PlotsRegister: Attempt to register incorrect histogram ids range: " << min << " to" << max << endl;
        exit(1);
    }
    
    if (min < 1 || max > 7999) {
        cerr << "  PlotsRegister: Attempt to register histogram ids: " << min << " to " << max << endl;
        cerr << "  PlotsRegister: Valid range is 1 to 7999" << endl;
        exit(1);
    }
    
    if (CheckRange(min, max)) {
        cerr << "  PlotsRegister: Attempt to register histogram ids: " << min << " to " << max << endl;
        cerr << "  PlotsRegister: This range is already registered." << endl;
        exit(EXIT_FAILURE);
    }
    
    reg.push_back( std::pair<int, int>(min, max) );

    //cout << "  PlotsRegister: Histogram ids " << min << " to " << max << " registered." << endl;
    return true;        
}
