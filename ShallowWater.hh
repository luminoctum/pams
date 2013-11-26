#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "Grid.hh"
#include <vector>

template <int order, class patch>
class ShallowWater{
typedef StagGridx<order/2> HalfGridx;
typedef StagGridy<order/2> HalfGridy;
typedef Grid<order/2> MainGrid;
protected:
    FiniteDifference<1> diff;
    float dx, dy, f;
public:
    inline void advection(
            const HalfGridx &uflux, 
            const HalfGridy &vflux, 
            MainGrid &result) const{
        result.main_t() += diff.x(uflux.main(), dx) + diff.y(vflux.main(), dy);
    }
    inline void self_advection(
            const MainGrid &phi,
            HalfGridx &uflux, 
            HalfGridy &vflux) const{
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
            HalfGridx &uflux,
            HalfGridy &vflux) const{
        uflux.main_t() += f * vflux.mainq();
        vflux.main_t() += - f * uflux.mainq();
    }
    inline void gradient(
            const MainGrid &phi,
            HalfGridx &vx,
            HalfGridy &vy) const{
        vx.main_t() += diff.x(0.5 * phi.extendx() * phi.extendx(), dx);
        vy.main_t() += diff.y(0.5 * phi.extendy() * phi.extendy(), dy);
    }
    void forward(std::vector<patch> &state){
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
