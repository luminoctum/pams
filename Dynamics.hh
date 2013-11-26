#ifndef DYNAMICS
#define DYNAMICS
#include "GridVariable.hh"

template<int order>
class Dynamics{
typedef GridVariable<order/2> Grid;
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
            const Grid &uwind, 
            const Grid &vwind, 
            Grid &result) const{
        result.main_t() += (
            diff.x(uwind.main() * result.mainx(), dx) 
            + diff.y(vwind.main() * result.mainy(), dy)
            );
    }
    inline void self_advection(
            Grid &uwind, 
            Grid &vwind) const{
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
            Grid &uwind,
            Grid &vwind) const{
        uwind.main_t() += f * vwind.v_to_u();
        vwind.main_t() += - f * uwind.u_to_v();
    }
    inline void gradient(
            const Grid &phi,
            Grid &vx,
            Grid &vy) const{
        vx.main_t() += diff.x(phi.extendx(), dx);
        vy.main_t() += diff.y(phi.extendy(), dy);
    }
};

#endif
