#ifndef ODEINT
#define ODEINT

template<int order>
class RungeKutta{
    // var is patch
    float op1 = {0.5, 0.5, 1.}
    float op2 = {1./6., 1./3., 1./3., 1./3.}
    for (int i = 0; i < order; i++){
        forward(q[i]);
        q[i + 1] = q[i];
        q[i + 1].update(op[i] * dt);
    }
    for (int i = 0; i < order; i++){
        var.value_t += op2[i] * q[i].value_t;
    }
    var.update(dt);
};

#endif
