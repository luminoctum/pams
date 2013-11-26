#ifndef DYNAMICS
#define DYNAMICS
#include "TileVariable.hh"

template<int order>
class Dynamics{
typedef TileVariable<order/2> Tile;
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
            const Tile &uwind, 
            const Tile &vwind, 
            Tile &result) const{
        result.main_t() += (
            diff.x(uwind.main() * result.mainx(), dx) 
            + diff.y(vwind.main() * result.mainy(), dy)
            );
    }
    inline void self_advection(
            Tile &uwind, 
            Tile &vwind) const{
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
            Tile &uwind,
            Tile &vwind) const{
        uwind.main_t() += f * vwind.v_to_u();
        vwind.main_t() += - f * uwind.u_to_v();
    }
    inline void gradient(
            const Tile &phi,
            Tile &vx,
            Tile &vy) const{
        vx.main_t() += diff.x(phi.extendx(), dx);
        vy.main_t() += diff.y(phi.extendy(), dy);
    }
};

#endif
