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
    ArrayXXf b(10, 10);
    Grid p1(&b, 1, 4, 1, 4), 
         p2(&b, 4, 7, 1, 4), 
         p3(&b, 1, 7, 1, 4);
    b.setRandom();
    cout << b << endl << endl;
    cout << p1.main() << endl << endl;
    cout << p2.main() << endl << endl;
    cout << p3.main() << endl << endl;
    cout << p1.mainx() << endl << endl;
    cout << p2.mainx() << endl << endl;
    cout << p3.mainx() << endl << endl;
}

void test_staggered_grid(){
    Grid p1;
    ArrayXXf a(5, 6);
    a << 0,2,3,4,1,0,
        0,5,6,4,2,0,
        0,8,9,4,3,0,
        0,2,1,4,2,0,
        0,5,1,4,6,0;
    StagGridy x1(&a, 1, 4, 1, 5);
    //p1.domain = &a;
    //p1.istart = 1; p1.iend = 3;
    //p1.jstart = 1; p1.jend = 4;
    cout << a << endl << endl;
    //cout << p1.main() << endl << endl;
    //cout << p1.mainq() << endl << endl;
    //x1.main() += 2;
    //cout << x1.mainy() << endl << endl;
    cout << x1.main() << endl << endl;
    cout << x1.mainy() << endl << endl;
    cout << x1.mainq() << endl << endl;
}

void test_dynamics(){
    ArrayXXf a(10, 10);
    ArrayXXf u(11, 10);
    ArrayXXf v(10, 11);
    Dynamics dyn;
    a.setRandom();
    u.setRandom(); 
    u.row(0).setZero(); u.row(u.rows() - 1).setZero();
    u.row(1).setZero(); u.row(u.rows() - 2).setZero();
    v.setRandom(); 
    v.col(0).setZero(); v.col(v.cols() - 1).setZero();
    v.col(1).setZero(); v.col(v.cols() - 2).setZero();

    Grid pa(&a, 1, 9, 1, 9);
    StagGridx pu(&u, 1, 10, 1, 9);
    StagGridy pv(&v, 1, 9, 1, 10);
    //cout << v << endl << endl;
    //cout << pa.main() << endl << endl;
    //cout << "sum = " << pa.main().sum() << endl << endl;
    cout << pu.main() << endl << endl;
    cout << pv.main() << endl << endl;
    //dyn.advection(1, pu, pv, pa);
    //pa.update();
    //cout << a << endl << endl;
    //cout << "new sum = " << pa.main().sum() << endl << endl;
    /*
    dyn.self_advection(1, pu, pv);
    pu.update(); pv.update();
    cout << pu.main() << endl << endl;
    cout << pv.main() << endl << endl;
    */
    Grid p1a(&a, 1, 5, 1, 9);
    StagGridx p1u(&u, 1, 6, 1, 9);
    StagGridy p1v(&v, 1, 5, 1, 10);
    Grid p2a(&a, 5, 9, 1, 9);
    StagGridx p2u(&u, 5, 10, 1, 9);
    StagGridy p2v(&v, 5, 9, 1, 10);
    //dyn.advection(1, p1u, p1v, p1a);
    //dyn.advection(1, p2u, p2v, p2a);
    //p1a.update(); p2a.update();
    //cout << "new sum = " << pa.main().sum() << endl << endl;
    //cout << a << endl << endl;
    dyn.self_advection(1, p1u, p1v);
    dyn.self_advection(1, p2u, p2v);
    p1u.update(); p1v.update();
    p2u.update(); p2v.update();
    cout << pu.main() << endl << endl;
    cout << pv.main() << endl << endl;
}

int main(){
    //internal::set_is_malloc_allowed(false);
    //test_finite_difference();
    //test_nonstaggered_grid();
    //test_staggered_grid();
    //test_dynamics();
}
