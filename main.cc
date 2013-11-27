#include <iostream>
#include <fstream>
#include <vector>
#include "patch_variable.hh"
#include "time_stepping.hh"
#include "module_shallow_water.hh"
#include "include.hh"
using namespace std;

int main(){
    Patch a(6, 6, 'i');
    Patch u(7, 6, 'x');
    Patch v(6, 7, 'y');
    ForwardModel forward;
    StateVector state;
    Runge_Kutta<4, StateVector> stepper;

    u.setLeftRightZero();
    v.setBottomTopZero();
    state.push_back(a);
    state.push_back(u);
    state.push_back(v);
    for (size_t i = 0; i < state.size(); i++){
        state[i].reset_ptr();
        cout << state[i].value << endl << endl;
    }
    cout << state[0].value.sum() << endl << endl;

    stepper.do_step(forward, state, 0.5);

    for (size_t i = 0; i < state.size(); i++){
        state[i].update(1.0);
        cout << state[i].value << endl << endl;
    }
    cout << state[0].value.sum() << endl << endl;


}
