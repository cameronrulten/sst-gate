#ifndef SST_GATE_GLOBALS_H
#define SST_GATE_GLOBALS_H

#include "AOpticsManager.h"

//define useful units: basically the system works in cm and the other units convert to centimeters
static const double cm = AOpticsManager::cm(); //centimeters
static const double mm = AOpticsManager::mm(); //millimeters
static const double um = AOpticsManager::um(); //micrometers
static const double nm = AOpticsManager::nm(); //nanometers
static const double  m = AOpticsManager::m();  //meters
static const double  km = AOpticsManager::km();//kilometers

struct clonable {
    virtual ~clonable() {}
    virtual clonable* clone() const = 0;
};

#endif
