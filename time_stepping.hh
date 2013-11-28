#ifndef TIMESTEPPING
#define TIMESTEPPING

template<int order, class StateType>
class Runge_Kutta{
protected:
    StateType buffer;
    float wgt1[order - 1], wgt2[order];
public:
    Runge_Kutta() : 
        wgt1({0.5, 0.5, 1.}),
        wgt2({1./6., 1./3., 1./3., 1./6.}){}
    template<class ForwardModel>
    void do_step(ForwardModel forward, StateType &state, float dt){
        buffer = state;
        for (size_t i = 0; i < state.size(); i++) buffer[i].reset_ptr();
        for (int i = 0; i < order - 1; i++){
            forward(buffer);
            for (size_t j = 0; j < state.size(); j++){
                state[j].value_t += wgt2[i] * buffer[j].value_t;
                buffer[j].value = state[j].value;
                buffer[j].update(wgt1[i] * dt);
            }
        }
        // final step
        forward(buffer);
        for (size_t j = 0; j < state.size(); j++){
            state[j].value_t += wgt2[order - 1] * buffer[j].value_t;
            state[j].update(dt);
        }
    }
};

#endif
