#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "GridVariable.hh"
#include <vector>

template <int order, int ntile_x, int ntile_y>
class ShallowWater{
typedef GridVariable<order/2> Grid;
typedef PatchVariable<GridVariable, order/2, ntile_x, ntile_y> Patch;
typedef std::vector<Patch> StateVector;
protected:
    FiniteDifference<1> diff;
    float dx, dy, f;
public:
    ShallowWater(){
        dx = 1.0;
        dy = 1.0;
        f = 1.0;
    }
    inline void advection(
            const Grid &uflux, 
            const Grid &vflux, 
            Grid &result) const{
        result.main_t() += diff.x(uflux.main(), dx) + diff.y(vflux.main(), dy);
    }
    inline void self_advection(
            const Grid &phi,
            Grid &uflux, 
            Grid &vflux) const{
        uflux.main_t() += (
                diff.x(uflux.mainx() * uflux.mainx() / phi.extendx(), dx) 
                + diff.y(uflux.mainy() * vflux.mainx() / phi.mainq(), dy)
                );
        vflux.main_t() += (
                diff.y(vflux.mainy() * vflux.mainy() / phi.extendy(), dy) 
                + diff.x(uflux.mainy() * vflux.mainx() / phi.mainq(), dx)
                );
    }
    inline void coriolis(
            float f,
            Grid &uflux,
            Grid &vflux) const{
        uflux.main_t() += f * vflux.v_to_u();
        vflux.main_t() += - f * uflux.u_to_v();
    }
    inline void gradient(
            const Grid &phi,
            Grid &vx,
            Grid &vy) const{
        vx.main_t() += diff.x(0.5 * phi.extendx() * phi.extendx(), dx);
        vy.main_t() += diff.y(0.5 * phi.extendy() * phi.extendy(), dy);
    }
    void operator()(StateVector &state){
        // 0 : phi
        // 1 : uwind
        // 2 : vwind
        #pragma omp parallel for
        for (int i = 0; i < NTILES_IN_X * NTILES_IN_Y; i++){
            advection(state[1](i), state[2](i), state[0](i));
            self_advection(state[0](i), state[1](i), state[2](i));
            gradient(state[0](i), state[1](i), state[2](i));
            coriolis(f, state[1](i), state[2](i));
        }
    }
};

#endif
