#define EIGEN_RUNTIME_NO_MALLOC
#include <iostream>
#include <fstream>
#include <vector>
#include "PatchVariable.hh"
#include "Dynamics.hh"
#include "ShallowWater.hh"
#include "TimeStepping.hh"
//#include "Advection.hh"

using namespace std;

void test_finite_difference(){
    FiniteDifference<1> diff;
    FiniteInterpolation<4> half;
    HyperDifference<4> hypr;
    ArrayXXf a(5, 5);
    a << 1,2,3,4,1,
        4,5,6,4,2,
        7,8,9,4,3,
        7,2,1,4,2,
        2,5,1,4,6;
    cout << a << endl;
    cout << diff.x(a, 2) << endl;
    cout << diff.y(a, 2) << endl;
    cout << half.x(a) << endl;
    cout << half.y(a) << endl;
    cout << hypr.x(a, 2) << endl;
    cout << hypr.y(a, 2) << endl;
}

void test_nonstaggered_grid(){
    PatchVariable<GridVariable,1,2,2> b(6, 6, 'i');
    cout << b.value << endl << endl;
    cout << b.tile[0][0].main() << endl << endl;
    cout << b.tile[0][1].main() << endl << endl;
    cout << b.tile[1][0].main() << endl << endl;
    cout << b.tile[1][1].main() << endl << endl;
    cout << b.tile[0][0].mainx() << endl << endl;
    cout << b.tile[0][0].mainy() << endl << endl;
}

void test_staggered_grid(){
    PatchVariable<GridVariable,1,2,2> a(7, 6, 'x');
    cout << a.value << endl << endl;
    cout << a.value_t << endl << endl;
    cout << a.tile[0][0].main() << endl << endl;
    cout << a.tile[0][1].main() << endl << endl;
    cout << a.tile[1][0].main() << endl << endl;
    cout << a.tile[1][1].main() << endl << endl;
    cout << a.tile[0][0].mainx() << endl << endl;
    cout << a.tile[0][0].mainy() << endl << endl;
    cout << a.tile[0][0].u_to_v() << endl << endl;
    a.tile[0][0].main() += 2;
    cout << a.value << endl << endl;
}

void test_dynamics(){
    PatchVariable<GridVariable, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> a(6, 6, 'i');
    PatchVariable<GridVariable, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> u(7, 6, 'x');
    PatchVariable<GridVariable, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> v(6, 7, 'y');
    Dynamics<O_INTERP> dyn;
    u.setLeftRightZero();
    v.setBottomTopZero();

    cout << "==================== initial value ==================== " <<  endl;
    cout << a.value << endl << endl;
    cout << u.value << endl << endl;
    cout << v.value << endl << endl;
    cout << "Total Sum " << a.value.sum() << endl;


    cout << "==================== domain decomposition ==================== " <<  endl;
    
    //a.clean_t(); u.clean_t(); v.clean_t();
    for (int i = 0; i < NTILES_IN_X * NTILES_IN_Y; i++){
        dyn.advection(u(i), v(i), a(i));
        dyn.self_advection(u(i), v(i));
        dyn.coriolis(u(i), v(i));
        dyn.gradient(a(i), u(i), v(i));
    }
    a.update(1.); u.update(1.); v.update(1.);
    cout << a.value << endl << endl;
    cout << u.value << endl << endl;
    cout << v.value << endl << endl;
    cout << "Total Sum " << a.value.sum() << endl;
}

int main(){
    //internal::set_is_malloc_allowed(false);
    //test_finite_difference();
    //test_nonstaggered_grid();
    //test_staggered_grid();
    //test_dynamics();
    //
    typedef PatchVariable<GridVariable, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> Patch;
    typedef std::vector<Patch> StateVector;
    Patch a(6, 6, 'i');
    Patch u(7, 6, 'x');
    Patch v(6, 7, 'y');
    u.setLeftRightZero();
    v.setBottomTopZero();
    StateVector state;
    state.push_back(a);
    state.push_back(u);
    state.push_back(v);
    for (int i = 0; i < state.size(); i++){
        state[i].reset_ptr();
        cout << state[i].value << endl << endl;
    }
    cout << state[0].value.sum() << endl << endl;
    ShallowWater<O_INTERP, NTILES_IN_X, NTILES_IN_Y> forward;
    forward(state);
    for (int i = 0; i < state.size(); i++){
        state[i].update(1.0);
        cout << state[i].value << endl << endl;
    }
    cout << state[0].value.sum() << endl << endl;
    Runge_Kutta<4, 
        ShallowWater<O_INTERP, NTILES_IN_X, NTILES_IN_Y>,
        StateVector> stepper;
    stepper(state, 0, 1, 0.5);

    /*
    TimeStepper<4, StateVector, ShallowWater> stepper;
    for (int i = 0; i < 5; i++){
        cout << u.value << endl << endl;
        cout << v.value << endl << endl;
        cout << a.value << endl << endl;
        cout << "Total Sum: " << a.value.sum() << endl;
        forward(u, v, a);
        cout << "Forward step #" << i << " :"<< endl;
        cout << u.value << endl << endl;
        cout << v.value << endl << endl;
        cout << a.value << endl << endl;
        cout << "Total Sum: " << a.value.sum() << endl;
    }
    */
}
