#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "module_base.hh"
#include "tile_variable.hh"

template <int order, int ntile_x, int ntile_y>
class ShallowWater : public ModuleBase{
typedef TileVariable<order/2> Tile;
protected:
    FiniteDifference<1> diff;
    HyperDifference<4> dissip;
    float f;
public:
    ShallowWater() : ModuleBase(){
        f = 1E-4;
    }
    inline void advection(
            const Tile &uflux, 
            const Tile &vflux, 
            Tile &result) const{
        result.main_t() -= diff.x(uflux.main(), dx) + diff.y(vflux.main(), dy);
    }
    inline void self_advection(
            const Tile &phi,
            Tile &uflux, 
            Tile &vflux) const{
        uflux.main_t() -= (
                diff.x(uflux.mainx() * uflux.mainx() / phi.extendx(), dx) 
                + diff.y(uflux.mainy() * vflux.mainx() / phi.mainq(), dy)
                );
        vflux.main_t() -= (
                diff.y(vflux.mainy() * vflux.mainy() / phi.extendy(), dy) 
                + diff.x(uflux.mainy() * vflux.mainx() / phi.mainq(), dx)
                );
    }
    inline void coriolis(
            float f,
            Tile &uflux,
            Tile &vflux) const{
        uflux.main_t() += f * vflux.v_to_u();
        vflux.main_t() += - f * uflux.u_to_v();
    }
    inline void gradient(
            const Tile &phi,
            Tile &vx,
            Tile &vy) const{
        vx.main_t() -= diff.x(0.5 * phi.extendx() * phi.extendx(), dx);
        vy.main_t() -= diff.y(0.5 * phi.extendy() * phi.extendy(), dy);
    }
    inline void diffusion(
            Tile &phi,
            Tile &uflux,
            Tile &vflux) const{
        phi.main_t() += 1E-4 * (dissip.x(phi.main()) + dissip.y(phi.main()));
        uflux.main_t() += 1E-4 * (dissip.x(uflux.main()) + dissip.y(uflux.main()));
        vflux.main_t() += 1E-4 * (dissip.x(vflux.main()) + dissip.y(vflux.main()));
    }
    inline void tracer(
            const Tile &phi,
            const Tile &uflux,
            const Tile &vflux,
            Tile &tracer) const {
        tracer.main_t() -= (
                diff.x(uflux.main() * tracer.mainx() / phi.mainx(), dx)
                + diff.y(vflux.main() * tracer.mainy() / phi.mainy(), dy)
                );
    }

    inline void total_energy(
            const Tile &phi,
            const Tile &uwind,
            const Tile &vwind,
            Tile &energy){
        energy.main() = 0.5 * (
            phi.main() * phi.main()
            + uwind.mainx(-1) * uwind.mainx(-1) / phi.main()
            + vwind.mainy(-1) * vwind.mainy(-1) / phi.main()
            );
    }
    template<typename StateType>
    void operator()(StateType &state){
        // 0 : phi
        // 1 : uwind
        // 2 : vwind
        // 3 : tracer
        #pragma omp parallel for
        for (int i = 0; i < ntile_x * ntile_y; i++){
            advection(state[1](i), state[2](i), state[0](i));
            self_advection(state[0](i), state[1](i), state[2](i));
            gradient(state[0](i), state[1](i), state[2](i));
            coriolis(f, state[1](i), state[2](i));
            tracer(state[0](i), state[1](i), state[2](i), state[3](i));
            //diffusion(state[0](i), state[1](i), state[2](i));
        }
    }
    template<typename StateType>
    void check_energy(StateType &state){
        #pragma omp parallel for
        for (int i = 0; i < ntile_x * ntile_y; i++){
            total_energy(state[0](i), state[1](i), state[2](i), state[3](i));
        }
    }
};

#endif
