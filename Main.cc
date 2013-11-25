#define EIGEN_RUNTIME_NO_MALLOC
#include <iostream>
#include <fstream>
#include <vector>
#include "PatchVariable.hh"
#include "Dynamics.hh"
#include "ForwardStateVector.hh"
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
    PatchVariable<Grid,1,2,2> b(6, 6);
    cout << b.value << endl << endl;
    cout << b.tile[0][0].main() << endl << endl;
    cout << b.tile[0][1].main() << endl << endl;
    cout << b.tile[1][0].main() << endl << endl;
    cout << b.tile[1][1].main() << endl << endl;
    cout << b.tile[0][0].mainx() << endl << endl;
    cout << b.tile[0][0].mainy() << endl << endl;
}

void test_staggered_grid(){
    PatchVariable<StagGridx,1,2,2> a(7, 6);
    cout << a.value << endl << endl;
    cout << a.value_t << endl << endl;
    cout << a.tile[0][0].main() << endl << endl;
    cout << a.tile[0][1].main() << endl << endl;
    cout << a.tile[1][0].main() << endl << endl;
    cout << a.tile[1][1].main() << endl << endl;
    cout << a.tile[0][0].mainx() << endl << endl;
    cout << a.tile[0][0].mainy() << endl << endl;
    cout << a.tile[0][0].mainq() << endl << endl;
    a.tile[0][0].main() += 2;
    cout << a.value << endl << endl;
}

void test_dynamics(){
    PatchVariable<Grid, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> a(6, 6);
    PatchVariable<StagGridx, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> u(7, 6);
    PatchVariable<StagGridy, O_INTERP/2, NTILES_IN_X, NTILES_IN_Y> v(6, 7);
    Dynamics<O_INTERP> dyn;
    u.setLeftRightZero();
    v.setBottomTopZero();
    /*
    u.value.row(0).setZero(); 
    u.value.row(u.nrows - 1).setZero();
    u.value.row(1).setZero(); 
    u.value.row(u.nrows - 2).setZero();
    v.value.col(0).setZero(); 
    v.value.col(v.ncols - 1).setZero();
    v.value.col(1).setZero(); 
    v.value.col(v.ncols - 2).setZero();
    */

    cout << "==================== initial value ==================== " <<  endl;
    cout << a.value << endl << endl;
    cout << u.value << endl << endl;
    cout << v.value << endl << endl;
    cout << "Total Sum " << a.value.sum() << endl;
    //for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++)
    //    cout << u.tile[i][j].main() << endl << endl;

    /*
    cout << "==================== original domain ==================== " <<  endl;
    dyn.advection(pu, pv, pa);
    dyn.self_advection(pu, pv);
    cout << "sum = " << pa.main_t().sum() << endl << endl;
    cout << pa.main_t()<< endl << endl;
    cout << pu.main_t()<< endl << endl;
    cout << pv.main_t()<< endl << endl;
    */

    cout << "==================== domain decomposition ==================== " <<  endl;
    
    //a.clean_t(); u.clean_t(); v.clean_t();
    /*
    for (int i = 0; i < NTILES_IN_X; i++) for (int j = 0; j < NTILES_IN_Y; j++){
        dyn.advection(u.tile[i][j], v.tile[i][j], a.tile[i][j]);
        dyn.self_advection(u.tile[i][j], v.tile[i][j]);
        dyn.coriolis(u.tile[i][j], v.tile[i][j]);
        dyn.gradient(a.tile[i][j], u.tile[i][j], v.tile[i][j]);
    }*/
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
    MainPatch a(6, 6);
    StagPatchx u(7, 6);
    StagPatchy v(6, 7);
    u.setLeftRightZero();
    v.setBottomTopZero();
    ForwardStateVector<ShallowWater> forward;
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
}
