#include "include.hh"
using namespace std;

int main(){
    float time;
    ForwardModel model;
    Runge_Kutta<4, StateVector> stepper;
    Patch a("phi", model.nrows, model.ncols, 'i', PERIODIC);
    Patch u("uwind", model.nrows + 1, model.ncols, 'x', PERIODIC);
    Patch v("vwind", model.nrows, model.ncols + 1, 'y', PERIODIC);
    /*
    Patch a("phi", model.nrows, model.ncols, 'i', DIRICHLET);
    Patch u("uwind", model.nrows + 1, model.ncols, 'x', DIRICHLET);
    Patch v("vwind", model.nrows, model.ncols + 1, 'y', DIRICHLET);
    */
    StateVector state{a,u,v};
    for (size_t i = 0; i < state.size(); i++){
        state[i].reset_ptr();
        state[i].value = model.ncVar[state[i].name];
        state[i].update(0.);
    }
    for (time = model.start; time < model.end; time += model.step){
        cout << time << endl;
        model.ncwrite(state, time);
        /*
        model(state);
        for (size_t i = 0; i < state.size(); i++){
            state[i].update(model.step);
        }*/
        stepper.do_step(model, state, model.step);
    }
    model.ncwrite(state, time);
}
