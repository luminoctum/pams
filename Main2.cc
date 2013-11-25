#define EIGEN_RUNTIME_NO_MALLOC
#include <iostream>
#include <fstream>
#include "Grid.hh"
//#include "Dynamics.hh"
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
    PatchVariable b;
    b.value.setRandom(10, 10);
    b.tendency.setZero(10, 10);
    cout << b.value << endl << endl;
    cout << b.tendency << endl << endl;
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

/*
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
*/

int main(){
    //internal::set_is_malloc_allowed(false);
    //test_finite_difference();
    test_nonstaggered_grid();
    //test_staggered_grid();
    //test_dynamics();
}
