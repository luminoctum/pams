#ifndef DYNAMICS
#define DYNAMICS
#include "Grid.hh"

class Dynamics{
protected:
    FiniteDifference<1> diff;
    FiniteInterpolation<2> half;
public:
    inline void advection(
            float dtdx,
            StagGridx &uwind, 
            StagGridy &vwind, 
            Grid &result){
        result.main_t() += dtdx * (
            diff.x(uwind.main() * result.mainx()) + diff.y(vwind.main() * result.mainy())
            );
    }
    inline void self_advection(
            float dtdx,
            StagGridx &uwind, 
            StagGridy &vwind){
        uwind.main_t() += dtdx * (
                diff.x(uwind.mainx() * uwind.mainx()) 
                + diff.y(uwind.mainy() * vwind.mainx())
                );
        vwind.main_t() += dtdx * (
                diff.y(vwind.mainy() * vwind.mainy()) 
                + diff.x(uwind.mainy() * vwind.mainx())
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
