#ifndef UEGO_H
#define UEGO_H

////////////////////////////////////////////////////////////
// $Id: uego.h,v 2.6 1998/04/01 21:05:18 jelasity Exp $
// uego.h
// collects different headers for simplification of usage;
// should be included from all modules
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 02 20 ini.h -> uegoini.h
//		98 03 10 nreal.h nbin.h
////////////////////////////////////////////////////////////
#include "functions/Extras.h"
#define UEGO_VERSION Extras::getName().c_str() //__DATE__ " " __TIME__ //"2.6"

// --- for sending errors and warnings, must be defined in a user interface
void message(const char *, short);

// --- defines for message
#define MSG_INFORMATION 2
#define MSG_ERROR 1
#define MSG_NOTHING 0

#include <cstdint>
// --- define to keep source clean
#define INI (Master::ini())

#include <math.h>
#include <stdio.h>

class SearchSpElement;
class SpeciesList;
class Master;

#include "master.h"
#include "model/Molecule.h"
#include "uegoini.h"
#include "uegorand.h"

#include "nmol.h"
#include "searchsp.h"
#include "speclist.h"
#include <pthread.h>

#endif
