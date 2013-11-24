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
    cout << diff.x(a) << endl;
    cout << diff.y(a) << endl;
    cout << half.x(a) << endl;
    cout << half.y(a) << endl;
    cout << hypr.x(a) << endl;
    cout << hypr.y(a) << endl;
}

void test_nonstaggered_grid(){
    ProgVariable b(10, 10);
    Grid p1(b, 1, 4, 1, 4), 
         p2(b, 4, 7, 1, 4), 
         p3(b, 1, 7, 1, 4);
    cout << p1.main() << endl << endl;
    cout << p2.main() << endl << endl;
    cout << p3.main() << endl << endl;
    cout << p1.mainx() << endl << endl;
    cout << p2.mainx() << endl << endl;
    cout << p3.mainx() << endl << endl;
}

void test_staggered_grid(){
    ProgVariable a(5, 6);
    a.value << 0,2,3,4,1,0,
        0,5,6,4,2,0,
        0,8,9,4,3,0,
        0,2,1,4,2,0,
        0,5,1,4,6,0;
    StagGridy x1(a, 1, 4, 1, 5);
    cout << a.value << endl << endl;
    cout << a.tendency << endl << endl;
    //Grid p1;
    //p1.domain = &a;
    //p1.istart = 1; p1.iend = 3;
    //p1.jstart = 1; p1.jend = 4;
    //cout << p1.main() << endl << endl;
    //cout << p1.mainq() << endl << endl;
    //x1.main() += 2;
    //cout << x1.mainy() << endl << endl;
    cout << x1.main() << endl << endl;
    cout << x1.mainy() << endl << endl;
    cout << x1.mainq() << endl << endl;
}

void test_dynamics(){
    ProgVariable a(10, 10);
    ProgVariable u(11, 10);
    ProgVariable v(10, 11);
    Dynamics dyn;
    u.value.row(0).setZero(); 
    u.value.row(u.nrows - 1).setZero();
    u.value.row(1).setZero(); 
    u.value.row(u.nrows - 2).setZero();
    v.value.col(0).setZero(); 
    v.value.col(v.ncols - 1).setZero();
    v.value.col(1).setZero(); 
    v.value.col(v.ncols - 2).setZero();

    cout << "==================== initial value ==================== " <<  endl;
    Grid pa(a, 1, 9, 1, 9);
    StagGridx pu(u, 1, 10, 1, 9);
    StagGridy pv(v, 1, 9, 1, 10);
    cout << pa.main() << endl << endl;
    cout << pu.main() << endl << endl;
    cout << pv.main() << endl << endl;

    cout << "==================== original domain ==================== " <<  endl;
    dyn.advection(1, pu, pv, pa);
    dyn.self_advection(1, pu, pv);
    cout << "sum = " << pa.main_t().sum() << endl << endl;
    cout << pa.main() + pa.main_t()<< endl << endl;
    cout << pu.main() + pu.main_t()<< endl << endl;
    cout << pv.main() + pv.main_t()<< endl << endl;

    cout << "==================== domain decomposition ==================== " <<  endl;
    
    Grid p1a(a, 1, 5, 1, 9);
    Grid p2a(a, 5, 9, 1, 9);
    StagGridx p1u(u, 1, 6, 1, 9);
    StagGridy p1v(v, 1, 5, 1, 10);
    StagGridx p2u(u, 5, 10, 1, 9);
    StagGridy p2v(v, 5, 9, 1, 10);
    dyn.advection(1, p1u, p1v, p1a);
    dyn.advection(1, p2u, p2v, p2a);
    dyn.self_advection(1, p1u, p1v);
    dyn.self_advection(1, p2u, p2v);
    cout << p1a.main() + p1a.main_t() << endl;
    cout << p2a.main() + p2a.main_t() << endl << endl;
    cout << p1u.main() + p1u.main_t() << endl;
    cout << p2u.main() + p2u.main_t() << endl << endl;
    cout << p1v.main() + p1v.main_t()<< endl;
    cout << p2v.main() + p2v.main_t()<< endl << endl;

    cout << a.tendency << endl;
    cout << u.tendency << endl;
}

int main(){
    //internal::set_is_malloc_allowed(false);
    //test_finite_difference();
    //test_nonstaggered_grid();
    //test_staggered_grid();
    test_dynamics();
}
