#ifndef TILEVARIABLE
#define TILEVARIABLE
#include "finite_method.hh"

template <int halo>
class TileVariable{
protected:
    FiniteInterpolation<2 * halo> half;
    ArrayXXf *patch, *patch_t;
    int istart, jstart, iend, jend; // start and end index, end is not included
    int is, ie, js, je;
public:
    TileVariable(){}
    TileVariable(ArrayXXf *_patch, ArrayXXf *_patch_t, 
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
    inline Block<ArrayXXf> main_t() const{
        return patch_t->block(
                istart, jstart,
                iend - istart, jend - jstart
                );
    }
    inline ArrayXXf mainx(int ext = 1) const{
        ArrayXXf buffer;
        buffer = half.x(patch->block(
                is, jstart, 
                ie - is, jend - jstart
                )).block(istart - is - ext, 0, iend - istart + ext, jend - jstart);
        return buffer;
    }
    inline ArrayXXf mainy(int ext = 1) const{
        ArrayXXf buffer;
        buffer = half.y(patch->block(
                istart, js, 
                iend - istart, je - js
                )).block(0, jstart - js - ext, iend - istart, jend - jstart + ext);
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
    inline ArrayXXf u_to_v() const{
        ArrayXXf buffer;
        buffer = half.x(half.y(patch->block(
                is, js, 
                ie - is, je - js
                ))).block(istart - is, jstart - js - 1, iend - istart - 1, jend - jstart + 1);
        return buffer;
    }
    inline ArrayXXf v_to_u() const{
        ArrayXXf buffer;
        buffer = half.x(half.y(patch->block(
                is, js, 
                ie - is, je - js
                ))).block(istart - is - 1, jstart - js, iend - istart + 1, jend - jstart - 1);
        return buffer;
    }
};

#endif
