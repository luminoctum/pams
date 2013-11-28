#include "include.hh"
using namespace std;

int main(){
    float energy;
    Patch a(6, 6, 'i', PERIODIC);
    Patch u(7, 6, 'x', PERIODIC);
    Patch v(6, 7, 'y', PERIODIC);
    Patch e(6, 6, 'i', PERIODIC);
    /*
    Patch a(6, 6, 'i', DIRICHLET);
    Patch u(7, 6, 'x', DIRICHLET);
    Patch v(6, 7, 'y', DIRICHLET);
    Patch e(6, 6, 'i', DIRICHLET);
    a.value = 10 + a.value;
    u.setLeftRightZero();
    v.setBottomTopZero();
    */

    a.value = 10 + a.value;
    ForwardModel forward;
    StateVector state{a,u,v,e};
    Runge_Kutta<4, StateVector> stepper;
    for (size_t i = 0; i < state.size(); i++){
        state[i].reset_ptr();
        state[i].update(0.);
        cout << state[i].main() << endl << endl;
    }
    cout << "Mass = " << state[0].main().sum() << endl;
    energy = 0.5 * (
            state[1](0).mainx(-1) * state[1](0).mainx(-1) / state[0](0).main() 
            + state[2](0).mainy(-1) * state[2](0).mainy(-1) / state[0](0).main()
            + state[0](0).main() * state[0](0).main()
            ).sum();
    cout << "Energy = " << energy << endl << endl;
    //cout << "Energy = " << state[3].main().sum() << endl << endl;
    cout << "aa" << endl;
    cout << forward.dx << endl;
    cout << forward.dy << endl;
    cout << forward.nrows << endl;
    cout << forward.ncols << endl;
    cout << forward.start << endl;
    cout << forward.end << endl;
    cout << forward.step << endl;
    cout << forward.ncVar["phi"]<< endl;

    /*
    for (int i = 0; i < 10; i++){
        stepper.do_step(forward, state, 0.01);
        for (size_t i = 0; i < state.size(); i++){
            state[i].update(1.);
            cout << state[i].main() << endl << endl;
        }
        cout << "Mass = " << state[0].main().sum() << endl;
        cout << "Energy = " << state[3].main().sum() << endl << endl;
    }*/
}
