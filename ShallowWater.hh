#ifndef SHALLOWWATER
#define SHALLOWWATER
#include <vector>
#include "Dynamics.hh"

class ShallowWater{
protected:
    ArrayXXf uflux;
    ArrayXXf vflux;
    ArrayXXf ght;
    std::vector<StagGridx> uflux_d;
    std::vector<StagGridy> vflux_d;
    std::vector<Grid> ght_d;
}
#endif 
