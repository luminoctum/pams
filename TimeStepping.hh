#ifndef TIMESTEPPING
#define TIMESTEPPING
#include <array>

template<int order, class ForwardModel, class StateVector>
class Runge_Kutta{
protected:
    std::array<StateVector, order> q;
    float wgt1[order], wgt2[order];
    ForwardModel forward;
public:
    Runge_Kutta(){}
    void operator()(StateVector &state, float start, float end, float dt){
        for (float time = start; time < end; time += dt){
            q[0] = state;
            for (int i = 0; i < order - 1; i++){
                forward(q[i]);
                q[i + 1] = q[i];
                for (int j = 0; j < state.size(); j++){
                    q[i + 1][j].update(wgt1[i] * dt);
                    state[j].value_t += wgt2[j] * q[i][j].value_t;
                }
            }
            for (int j = 0; j < state.size(); j++) state[j].update(dt);
        }
    }
};

#endif
