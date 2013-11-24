#ifndef GRID
#define GRID
#include "NumericalMethod.hh"

class Grid{
protected:
    FiniteInterpolation<OINTERP> half;
    ArrayXXf **domain, buffer;
    std::string name;
    int istart, jstart, iend, jend; // start and end index, end is not included
    int is, ie, js, je;
    int halo; // halo width, at least 1

public:
    Grid(){};
    Grid(ArrayXXf **_domain, int _istart, int _iend, int _jstart, int _jend){
        halo        = OINTERP / 2;
        domain[0]   = _domain[0];
        domain[1]   = _domain[1];
        istart      = _istart;
        iend        = _iend;
        jstart      = _jstart;
        jend        = _jend;
        is          = MAX2(istart - halo, 0);
        ie          = MIN2(iend + halo, domain[0]->rows());
        js          = MAX2(jstart - halo, 0);
        je          = MIN2(jend + halo, domain[0]->cols());
    }
    inline Block<ArrayXXf> main(){
        return domain[0]->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline Block<ArrayXXf> mainx(){
        buffer = half.x(domain[0]->block(
                is, jstart, 
                ie - is, jend - jstart
                ));
        return buffer.block(istart - is - 1, 0, iend - istart + 1, jend - jstart);
    }
    inline Block<ArrayXXf> mainy(){
        buffer = half.y(domain[0]->block(
                istart, js, 
                iend - istart, je - js
                ));
        return buffer.block(0, jstart - js - 1, iend - istart, jend - jstart + 1);
    }
    inline Block<ArrayXXf> mainq(){
        buffer = half.x(half.y(domain[0]->block(
                is, js, 
                ie - is, je - js
                )));
        return buffer.block(istart - is - 1, jstart - js - 1, iend - istart + 1, jend - jstart + 1);
    }
};

class StagGridx : public Grid{
public:
    StagGridx() : Grid(){}
    StagGridx(ArrayXXf **_domain, int _istart, int _iend, int _jstart, int _jend):
        Grid(_domain, _istart, _iend, _jstart, _jend){}
    inline Block<ArrayXXf> mainq(){
        buffer = half.x(half.y(domain[0]->block(
                is, js, 
                ie - is, je - js
                )));
        return buffer.block(istart - is, jstart - js - 1, iend - istart - 1, jend - jstart + 1);
    }
};

class StagGridy : public Grid{
public:
    StagGridy() : Grid(){}
    StagGridy(ArrayXXf **_domain, int _istart, int _iend, int _jstart, int _jend):
        Grid(_domain, _istart, _iend, _jstart, _jend){}
    inline Block<ArrayXXf> mainq(){
        buffer = half.x(half.y(domain[0]->block(
                is, js, 
                ie - is, je - js
                )));
        return buffer.block(istart - is - 1, jstart - js, iend - istart + 1, jend - jstart - 1);
    }
};

#endif
