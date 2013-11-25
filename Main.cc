#define EIGEN_RUNTIME_NO_MALLOC
#include <iostream>
#include <fstream>
#include <vector>
#include "Grid.hh"
#include "Dynamics.hh"
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
}

void test_staggered_grid(){
    PatchVariable<StagGridy,1,2,2> a(6, 7);
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
    PatchVariable<Grid,1,3,3> a(11, 11);
    PatchVariable<StagGridx,1,3,3> u(12, 11);
    PatchVariable<StagGridy,1,3,3> v(11, 12);
    Dynamics<2> dyn;
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
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++){
        dyn.advection(u.tile[i][j], v.tile[i][j], a.tile[i][j]);
        dyn.self_advection(u.tile[i][j], v.tile[i][j]);
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
    test_dynamics();
}
