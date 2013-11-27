#ifndef TIMESTEPPING
#define TIMESTEPPING
#include <array>

template<int order, class StateType>
class Runge_Kutta{
protected:
    std::array<StateType, order> q;
    float wgt1[order], wgt2[order];
public:
    Runge_Kutta(){}
    template<class ForwardModel>
    void do_step(ForwardModel forward, StateType &state, float dt){
        q[0] = state;
        for (int i = 0; i < order - 1; i++){
            forward(q[i]);
            q[i + 1] = q[i];
            for (size_t j = 0; j < state.size(); j++){
                q[i + 1][j].update(wgt1[i] * dt);
                state[j].value_t += wgt2[j] * q[i][j].value_t;
            }
        }
        for (size_t j = 0; j < state.size(); j++) state[j].update(dt);
    }
};

#endif
