#ifndef TIMESTEPPING
#define TIMESTEPPING
#include "ShallowWater.hh"
#include <array>

template<int order, class Forward, class StateVector>
class Runge_Kutta{
protected:
    std::array<StateVector, order> q;
    float wgt1[order], wgt2[order];
public:
    Runge_Kutta(){}
    void operator()(StateVector &state, float start, float end, float dt){
        for (float time = start; time < end; time += dt){
            q[0] = state;
            ALARM("aa");
            ALARM(q[0][0].value);
            for (int i = 0; i < order - 1; i++){
                //Forward(q[i]);
                q[i + 1] = q[i];
                for (int j = 0; j < state.size(); j++){
                    ALARM(q[i][j].value);
                    //q[i + 1][j].update(wgt1[i] * dt);
                    //state[j].value_t += wgt2[j] * q[i][j].value_t;
                }
            }
            for (int j = 0; j < state.size(); j++){
                state[j].update(dt);
            }
        }
    }
};

#endif
