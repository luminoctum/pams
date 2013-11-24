#ifndef ODEINT
#define ODEINT

template<int order>
class RungeKutta{
    // var is patch
    q1 = var;
    q2 = var;
    q3 = var;
    q4 = var;
    int op1 = {0.5, 0.5, 1.}
    int op2 = {1., 2., 2., 1.}
    for (int i = 0; i < order; i++){
        forward(q[i], q[i + 1]);
        q[i + 1].update(op[i] * dt);
    }
    var.tendency = 1./6.*(q1.tendency + 2. * q2.tendency + 2. * q3.tendency + q4.tendency);
    var.update(dt);
};

#endif
