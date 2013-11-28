#include "include.hh"
using namespace std;

int main(){
    float time;
    ads::Timer timer;
    ForwardModel model;
    Runge_Kutta<4, StateVector> stepper;
    int nrows = model.nrows + O_INTERP;
    int ncols = model.ncols + O_INTERP;
    Patch a("phi", nrows, ncols, 'i', false, PERIODIC);
    Patch u("uwind", nrows + 1, ncols, 'x', true, PERIODIC);
    Patch v("vwind", nrows, ncols + 1, 'y', true, PERIODIC);
    Patch c("tracer", nrows, ncols, 'i', true, PERIODIC);
    /*
    Patch a("phi", model.nrows, model.ncols, 'i', DIRICHLET);
    Patch u("uwind", model.nrows + 1, model.ncols, 'x', DIRICHLET);
    Patch v("vwind", model.nrows, model.ncols + 1, 'y', DIRICHLET);
    */
    StateVector state{a,u,v,c}; // Zeroth variable has to be mass variable
    for (size_t i = 0; i < state.size(); i++){
        if (state[i].scale_by_mass){
            state[i].main() = state[0].main_stag(state[i].stag) * model.ncVar[state[i].name];
        } else {
            state[i].main() = model.ncVar[state[i].name];
        }
        state[i].update(0.);
    }
    cout << model << endl;
    printf("%-8s%-16s%-16s%-16s\n", "Steps", "Total Mass", "Model Time (s)", "Elapsed Time (s)");
    timer.tic();
    for (time = model.start; time < model.end; time += model.step){
        stepper.do_step(model, state, model.step);
        if (model.current % model.frame == 0) { 
            model.ncwrite(state, time + model.step);
        }
        if (model.current % (10 * model.frame) == 0){
            printf("%-8d%-16.3E%-16.3E%-16.3f\n", 
                    model.current, state[0].main().sum(), time, timer.toc()
            );
        }
    }
}
