#ifndef ODEINT
#define ODEINT

template<int order, typename StateVector>
class RungeKutta{
protected:
    StateVector q[order];
public:
    RungeKutta(){
        float op1 = {0.5, 0.5, 1.}
        float op2 = {1./6., 1./3., 1./3., 1./3.}
    }
    integrate(MainState m, StagState ufloat start, float end, float dt){
        for (float time = start; time < end; time += dt){
            q[0] = var;
            for (int i = 0; i < order; i++){
                forward(q[i]);
                q[i + 1] = q[i];
                q[i + 1].update(op[i] * dt);
            }
            for (int i = 0; i < order; i++){
                var.value_t += op2[i] * q[i].value_t;
            }
            var.update(dt);
        }
    }
};

#endif
