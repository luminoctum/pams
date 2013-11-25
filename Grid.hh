#ifndef GRID
#define GRID
#include "NumericalMethod.hh"

template <int halo = 1>
class Grid{
protected:
    FiniteInterpolation<2 * halo> half;
    ArrayXXf *patch, *patch_t;
    std::string name;
    int istart, jstart, iend, jend; // start and end index, end is not included
    int is, ie, js, je;
    int offset_x, offset_y;

public:
    Grid(){
        offset_x    = 0;
        offset_y    = 0;
    };
    Grid(ArrayXXf *_patch, ArrayXXf *_patch_t, 
            int _istart, int _jstart, int _ilen, int _jlen){
        set(_patch, _patch_t, _istart, _jstart, _ilen, _jlen);
    }
    void set(ArrayXXf *_patch, ArrayXXf *_patch_t, 
            int _istart, int _jstart, int _ilen, int _jlen){
        patch       = _patch;
        patch_t     = _patch_t;
        istart      = _istart;
        jstart      = _jstart;
        iend        = _istart + _ilen;
        jend        = _jstart + _jlen;
        is          = MAX2(istart - halo, 0);
        ie          = MIN2(iend + halo, patch->rows());
        js          = MAX2(jstart - halo, 0);
        je          = MIN2(jend + halo, patch->cols());
        offset_x    = 0;
        offset_y    = 0;
    }
    inline Block<ArrayXXf> main(){
        return patch->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline Block<ArrayXXf> main() const{
        return patch->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline Block<ArrayXXf> extendx() const{
        return patch->block(
                istart - 1, jstart,
                iend - istart + 2, jend - jstart
                );
    }
    inline Block<ArrayXXf> extendy() const{
        return patch->block(
                istart, jstart - 1,
                iend - istart, jend - jstart + 2
                );
    }
    inline Block<ArrayXXf> main_t(){
        return patch_t->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline ArrayXXf mainx() const{
        ArrayXXf buffer;
        buffer = half.x(patch->block(
                is, jstart, 
                ie - is, jend - jstart
                )).block(istart - is - 1, 0, iend - istart + 1, jend - jstart);
        return buffer;
    }
    inline ArrayXXf mainy() const{
        ArrayXXf buffer;
        buffer = half.y(patch->block(
                istart, js, 
                iend - istart, je - js
                )).block(0, jstart - js - 1, iend - istart, jend - jstart + 1);
        return buffer;
    }
    inline ArrayXXf mainq() const{
        ArrayXXf buffer;
        buffer = half.x(half.y(patch->block(
                is, js, 
                ie - is, je - js
                ))).block(istart - is - 1, jstart - js - 1, iend - istart + 1, jend - jstart + 1);
        return buffer;
    }
    inline int get_offset_x() const{
        return offset_x;
    };
    inline int get_offset_y() const{
        return offset_y;
    }
};

template <int halo = 1>
class StagGridx : public Grid<halo>{
public:
    StagGridx() : Grid<halo>(){
        this->offset_x    = 1;
        this->offset_y    = 0;
    }
    StagGridx(ArrayXXf *_patch, ArrayXXf *_patch_t, int _istart, int _iend, int _jstart, int _jend):
        Grid<halo>(_patch, _patch_t, _istart, _iend, _jstart, _jend){
            this->offset_x    = 1;
            this->offset_y    = 0;
    }
    inline ArrayXXf mainq() const{
        ArrayXXf buffer;
        buffer = this->half.x(this->half.y(this->patch->block(
                this->is, this->js, 
                this->ie - this->is, this->je - this->js
                ))).block(
                this->istart - this->is, this->jstart - this->js - 1, 
                this->iend - this->istart - 1, this->jend - this->jstart + 1
                );
        return buffer;
    }
};

template <int halo = 1>
class StagGridy : public Grid<halo>{
public:
    StagGridy() : Grid<halo>(){
        this->offset_x    = 0;
        this->offset_y    = 1;
    }
    StagGridy(ArrayXXf *_patch, ArrayXXf *_patch_t, int _istart, int _iend, int _jstart, int _jend):
        Grid<halo>(_patch, _patch_t, _istart, _iend, _jstart, _jend){
        this->offset_x    = 0;
        this->offset_y    = 1;
    }
    inline ArrayXXf mainq() const{
        ArrayXXf buffer;
        buffer = this->half.x(this->half.y(this->patch->block(
                this->is, this->js, 
                this->ie - this->is, this->je - this->js
                ))).block(
                this->istart - this->is - 1, this->jstart - this->js, 
                this->iend - this->istart + 1, this->jend - this->jstart - 1
                );
        return buffer;
    }
};

#endif
