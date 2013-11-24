#ifndef GRID
#define GRID
#include "NumericalMethod.hh"
#include "ProgVariable.hh"

class Grid{
protected:
    FiniteInterpolation<OINTERP> half;
    ArrayXXf *patch, *patch_t;
    ArrayXXf buffer;
    std::string name;
    int istart, jstart, iend, jend; // start and end index, end is not included
    int is, ie, js, je;
    int halo; // halo width, at least 1

public:
    Grid(){};
    Grid(ProgVariable &var, int _istart, int _iend, int _jstart, int _jend){
        halo        = OINTERP / 2;
        patch       = &var.value;
        patch_t     = &var.tendency;
        istart      = _istart;
        iend        = _iend;
        jstart      = _jstart;
        jend        = _jend;
        is          = MAX2(istart - halo, 0);
        ie          = MIN2(iend + halo, patch->rows());
        js          = MAX2(jstart - halo, 0);
        je          = MIN2(jend + halo, patch->cols());
    }
    inline Block<ArrayXXf> main(){
        return patch->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline Block<ArrayXXf> main_t(){
        return patch_t->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline Block<ArrayXXf> mainx(){
        buffer = half.x(patch->block(
                is, jstart, 
                ie - is, jend - jstart
                ));
        return buffer.block(istart - is - 1, 0, iend - istart + 1, jend - jstart);
    }
    inline Block<ArrayXXf> mainy(){
        buffer = half.y(patch->block(
                istart, js, 
                iend - istart, je - js
                ));
        return buffer.block(0, jstart - js - 1, iend - istart, jend - jstart + 1);
    }
    inline Block<ArrayXXf> mainq(){
        buffer = half.x(half.y(patch->block(
                is, js, 
                ie - is, je - js
                )));
        return buffer.block(istart - is - 1, jstart - js - 1, iend - istart + 1, jend - jstart + 1);
    }
};

class StagGridx : public Grid{
public:
    StagGridx() : Grid(){}
    StagGridx(ProgVariable &var, int _istart, int _iend, int _jstart, int _jend):
        Grid(var, _istart, _iend, _jstart, _jend){}
    inline Block<ArrayXXf> mainq(){
        buffer = half.x(half.y(patch->block(
                is, js, 
                ie - is, je - js
                )));
        return buffer.block(istart - is, jstart - js - 1, iend - istart - 1, jend - jstart + 1);
    }
};

class StagGridy : public Grid{
public:
    StagGridy() : Grid(){}
    StagGridy(ProgVariable &var, int _istart, int _iend, int _jstart, int _jend):
        Grid(var, _istart, _iend, _jstart, _jend){}
    inline Block<ArrayXXf> mainq(){
        buffer = half.x(half.y(patch->block(
                is, js, 
                ie - is, je - js
                )));
        return buffer.block(istart - is - 1, jstart - js, iend - istart + 1, jend - jstart - 1);
    }
};

#endif
