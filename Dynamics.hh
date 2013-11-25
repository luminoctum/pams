#ifndef DYNAMICS
#define DYNAMICS
#include "Grid.hh"

template<int order = 2>
class Dynamics{
typedef StagGridx<order/2> HalfGridx;
typedef StagGridx<order/2> HalfGridy;
typedef Grid<order/2> IntGrid;
protected:
    FiniteDifference<1> diff;
public:
    inline void advection(
            StagGridx<1> &uwind, 
            StagGridy<1> &vwind, 
            Grid<1> &result){
        result.main_t() += (
            diff.x(uwind.main() * result.mainx(), 1) 
            + diff.y(vwind.main() * result.mainy(), 1)
            );
    }
    inline void self_advection(
            StagGridx<1> &uwind, 
            StagGridy<1> &vwind){
        uwind.main_t() += (
                diff.x(uwind.mainx() * uwind.mainx(), 1) 
                + diff.y(uwind.mainy() * vwind.mainx(), 1)
                );
        vwind.main_t() += (
                diff.y(vwind.mainy() * vwind.mainy(), 1) 
                + diff.x(uwind.mainy() * vwind.mainx(), 1)
                );
    }
    inline void coriolis(
            float f,
            StagGridx<1> &uwind,
            StagGridy<1> &vwind){
        uwind.main_t() += f * vwind.mainq();
        vwind.main_t() += - f * uwind.mainq();
    }
    inline void gradient(
            Grid<1> &phi,
            StagGridx<1> &vx,
            StagGridy<1> &vy){
        vx.main_t() += diff(phi.mainx(), 1);
        vy.main_t() += diff(phi.mainy(), 1);
    }
};

#endif
