#ifndef DYNAMICS
#define DYNAMICS
#include "Grid.hh"

template<int order>
class Dynamics{
typedef StagGridx<order/2> HalfGridx;
typedef StagGridy<order/2> HalfGridy;
typedef Grid<order/2> MainGrid;
protected:
    FiniteDifference<1> diff;
    float dx, dy, f;
public:
    Dynamics(){
        dx  =   1.;
        dy  =   1.;
        f   =   2.;
    }
    inline void advection(
            const HalfGridx &uwind, 
            const HalfGridy &vwind, 
            MainGrid &result) const{
        result.main_t() += (
            diff.x(uwind.main() * result.mainx(), dx) 
            + diff.y(vwind.main() * result.mainy(), dy)
            );
    }
    inline void self_advection(
            HalfGridx &uwind, 
            HalfGridy &vwind) const{
        uwind.main_t() += (
                diff.x(uwind.mainx() * uwind.mainx(), dx) 
                + diff.y(uwind.mainy() * vwind.mainx(), dy)
                );
        vwind.main_t() += (
                diff.y(vwind.mainy() * vwind.mainy(), dy) 
                + diff.x(uwind.mainy() * vwind.mainx(), dx)
                );
    }
    inline void coriolis(
            HalfGridx &uwind,
            HalfGridy &vwind) const{
        uwind.main_t() += f * vwind.mainq();
        vwind.main_t() += - f * uwind.mainq();
    }
    inline void gradient(
            const MainGrid &phi,
            HalfGridx &vx,
            HalfGridy &vy) const{
        vx.main_t() += diff.x(phi.extendx(), dx);
        vy.main_t() += diff.y(phi.extendy(), dy);
    }
};

#endif
