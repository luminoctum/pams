#ifndef INCLUDE
#define INCLUDE
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ads/timer.h>
#include "patch_variable.hh"
#include "time_stepping.hh"

/* define domain decomposition */
#define     NTILES_IN_X     2
#define     NTILES_IN_Y     2
/* define interpolation order */
#define     O_INTERP        2
/* define constants */
#define     PI              3.1415926

/* typedef of major variables */
#include "module_shallow_water.hh"
typedef PatchVariable<O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> Patch;
typedef std::vector<Patch> StateVector;
typedef ShallowWater<O_INTERP, NTILES_IN_X, NTILES_IN_Y> ForwardModel;



#endif
