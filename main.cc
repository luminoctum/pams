#include "include.hh"
using namespace std;

int main(){
    float time;
    ads::Timer timer;
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
        state[i].value = model.ncVar[state[i].name];
        state[i].update(0.);
    }
    timer.tic();
    for (time = model.start; time < model.end; time += model.step){
        if (model.current % model.frame == 0){
            cout << setw(8) << left << model.current
                << setw(15) << left << setprecision(5) << time
                << setw(15) << left << setprecision(5) << timer.toc() 
                << endl;
            model.ncwrite(state, time);
        }
        stepper.do_step(model, state, model.step);
    }
    model.ncwrite(state, time);
}
