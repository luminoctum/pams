#include "include.hh"
#include "getopt.h"
//#include "unistd.h"
using namespace std;

int main(int argc, char *argv[]){
    int op = 0;
    long restart = 0;
    while ((op = getopt(argc, argv, "c:")) != -1)
        switch (op){
            case 'c':
                restart = atol(optarg);
                break;
            default:
                abort();
        }
    
    float time;
    ads::Timer timer;
    ForwardModel model;
    Runge_Kutta<4, StateVector> stepper;
    int nrows = model.nrows + O_INTERP;
    int ncols = model.ncols + O_INTERP;
    long steps = 0;
    Patch a("phi", nrows, ncols, 'i', false, PERIODIC, DIRICHLET);
    Patch u("uwind", nrows + 1, ncols, 'x', true, PERIODIC, DIRICHLET);
    Patch v("vwind", nrows, ncols + 1, 'y', true, PERIODIC, DIRICHLET);
    Patch c("tracer", nrows, ncols, 'i', true, PERIODIC);
    Patch e("energy", nrows, ncols, 'i', false, PERIODIC);
    StateVector state{a,u,v,c,e}; // Zeroth variable has to be mass variable
    for (size_t i = 0; i < state.size(); i++){
        if (state[i].scale_by_mass){
            state[i].main() = state[0].main_stag(state[i].stag) * model.ncVar[state[i].name];
        } else {
            state[i].main() = model.ncVar[state[i].name];
        }
        state[i].set_halo();
    }
    //cout << state[2].value << endl;
    cout << model << endl;
    printf("%-8s%-16s%-16s%-16s%-16s\n", 
            "Steps", "Total Mass", "Total Energy", "Model Time (s)", "Elapsed Time (s)"
    );
    timer.tic();
    // a bug in the restart run, model end is not updated
    for (time = model.start; time < model.end + restart; time += model.step){
        stepper.do_step(model, state, model.step);
        steps++;
        if (steps % model.frame == 0)
            model.ncwrite(state, time + model.step);
        if (steps % (10 * model.frame) == 0){
            model.check_energy(state);
            printf("%-8d%-16.3E%-16.3E%-16.3E%-16.3f\n", 
                    model.current, 
                    state[0].main().sum(), state[4].main().sum(),
                    time, timer.toc()
            );
        }
    }
}
