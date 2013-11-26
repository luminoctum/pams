#ifndef Forward
#define Forward
#include "ShallowWater.hh"

typedef PatchVariable<Grid, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> MainPatch;
typedef PatchVariable<StagGridx, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> StagPatchx;
typedef PatchVariable<StagGridy, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> StagPatchy;

template<int order, class patch>
class Runge_Kutta : public {
protected:
    std::vector<patch> q[4];
public:
    void operator()(std::vector<patch> &state, float start, float end, float dt){
        for (float time = start; time < end; time += dt){
            q[0] = state;
            for (int i = 0; i < order; i++){
                Forward(state);
                q[i + 1] = q[i];
                for (int j = 0; j < state.size(); j++){
                    q[i + 1][j].update(op[i] * dt);
                }
            }
            for (int j = 0; j < state.size(); j++){
                for (int i = 0; i < order; i++) 
                    state[j].value_t += op2[i] * q[i][j].value_t;
                state[j].update(dt);
            }
        }
    }
}

#endif
