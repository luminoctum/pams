#ifndef DYNAMICS
#define DYNAMICS
#include "Grid.hh"

class Dynamics{
protected:
    FiniteDifference<1> diff;
    FiniteInterpolation<2> half;
public:
    inline void advection(
            StagGridx &uwind, 
            StagGridy &vwind, 
            Grid &result){
        result.main_t() += (
            diff.x(uwind.main() * result.mainx(), 1) 
            + diff.y(vwind.main() * result.mainy(), 1)
            );
    }
    inline void self_advection(
            StagGridx &uwind, 
            StagGridy &vwind){
        uwind.main_t() += (
                diff.x(uwind.mainx() * uwind.mainx(), 1) 
                + diff.y(uwind.mainy() * vwind.mainx(), 1)
                );
        vwind.main_t() += (
                diff.y(vwind.mainy() * vwind.mainy(), 1) 
                + diff.x(uwind.mainy() * vwind.mainx(), 1)
                );
    }
    /*
    void coriolis(
            StagGridx &uwind,
            StagGridy &vwind){
        uwind.mainx() += f * vwind.mainx();
        vwind.mainx() += f * uwind.mainx();
    }
    */
};

#endif
