#ifndef PATCH_VARIABLE
#define PATCH_VARIABLE
#include "tile_variable.hh"

enum Boundary{DIRICHLET, NEUMANN, PERIODIC};

template <template <int>class TileType, int halo, int ntile_x, int ntile_y>
struct PatchVariable{
    std::string name;
    ArrayXXf value;
    ArrayXXf value_t;
    TileType<halo> tile[ntile_x][ntile_y];
    int patch_rows, patch_cols;
    int tile_rows, tile_cols;
    int offset_x, offset_y;
    char stag;
    Boundary bnd[4];

    PatchVariable(){};
    PatchVariable(std::string _name, int _nrows, int _ncols, char _stag,
            Boundary bnd1 = DIRICHLET, Boundary bnd2 = DIRICHLET,
            Boundary bnd3 = DIRICHLET, Boundary bnd4 = DIRICHLET){
        name = _name;
        value.setRandom(_nrows, _ncols);
        value_t.setZero(_nrows, _ncols);
        switch (_stag){
            case 'i':
                offset_x = 0;
                offset_y = 0;
                break;
            case 'x':
                offset_x = 1;
                offset_y = 0;
                break;
            case 'y':
                offset_x = 0;
                offset_y = 1;
                break;
            default:
                ASSERT_NOT_SUPPORTED;
        }
        int tx, ty;
        tx = _nrows - 2 * halo + offset_x * (ntile_x - 1);
        ty = _ncols - 2 * halo + offset_y * (ntile_y - 1);
        if ((tx % ntile_x) || (ty % ntile_y)){ASSERT_NOT_SUPPORTED;}
        tile_rows = tx / ntile_x;
        tile_cols = ty / ntile_y;
        for (int i = 0; i < ntile_x; i++)
            for (int j = 0; j < ntile_y; j++)
                tile[i][j].set(&value, &value_t, 
                        i * tile_rows + halo - i * offset_x, 
                        j * tile_cols + halo - j * offset_y,
                        tile_rows, tile_cols, _stag
                        );
        patch_rows = _nrows;
        patch_cols = _ncols;
        stag       = _stag;
        if (bnd2 == 0 && bnd3 == 0 && bnd4 == 0){
            for (int i = 0; i < 4; i++) bnd[i] = bnd1;
        } else if (bnd3 == 0 && bnd4 == 0){
            for (int i = 0; i < 2; i++) bnd[i] = bnd1;
            for (int i = 2; i < 4; i++) bnd[i] = bnd2;
        } else{
            bnd[0] = bnd1; bnd[1] = bnd2; bnd[2] = bnd3; bnd[3] = bnd4;
        }
        if (bnd[0] == PERIODIC){if (stag == 'x') { setLeftRightSame(); }}
        if (bnd[2] == PERIODIC){if (stag == 'y') { setBottomTopSame(); }}
    }
    inline TileType<halo>& operator()(int index){
        return tile[index / ntile_y][index % ntile_y];
    }
    inline Block<ArrayXXf> main(){
        return value.block(halo, halo, 
                patch_rows - 2 * halo, patch_cols - 2 * halo);
    }
    inline void reset_ptr(){
        for (int i = 0; i < ntile_x; i++)
            for (int j = 0; j < ntile_y; j++)
                tile[i][j].set(&value, &value_t, 
                        i * tile_rows + halo - i * offset_x, 
                        j * tile_cols + halo - j * offset_y,
                        tile_rows, tile_cols, stag
                        );
    }
    inline void setLeftRightZero(){
        for (int i = 0; i <= halo; i++){
            value.row(i).setZero();
            value.row(patch_rows - 1 - i).setZero();
        }
    }
    inline void setLeftRightSame(){
        ArrayXXf buffer = 
            0.5 * (value.row(halo) + value.row(patch_rows -  2 * halo));
        value.row(halo) = buffer;
        value.row(patch_rows - 2 * halo) = buffer;
    }
    inline void setBottomTopZero(){
        for (int i = 0; i <= halo; i++){
            value.col(i).setZero();
            value.col(patch_cols - 1 - i).setZero();
        }
    }
    inline void setBottomTopSame(){
        ArrayXXf buffer = 
            0.5 * (value.col(halo) + value.col(patch_cols -  2 * halo));
        value.col(halo) = buffer;
        value.col(patch_cols - 2 * halo) = buffer;
    }
    void update(float dt){
        // for now, the patch is the whole domain, halo update is merged into update
        if (stag == 'x'){
            for (int i = 1; i < ntile_x; i++)
                value_t.row(i * tile_rows + halo - i * offset_x) /= 2;
            if (bnd[0] == DIRICHLET){ value_t.row(halo).setZero(); }
            if (bnd[1] == DIRICHLET){ value_t.row(patch_rows - halo - 1).setZero(); }
        } else if (stag == 'y'){
            for (int i = 1; i < ntile_y; i++)
                value_t.col(i * tile_cols + halo - i * offset_y) /= 2;
            if (bnd[2] == DIRICHLET){ value_t.col(halo).setZero(); }
            if (bnd[3] == DIRICHLET){ value_t.col(patch_cols - halo - 1).setZero(); }
        }
        value = value + value_t * dt;
        value_t.setZero();
        if (bnd[0] == NEUMANN){
            for (int i = 0; i < halo; i++) value.row(i) = value.row(halo);
        } else if (bnd[0] == PERIODIC){
            for (int i = 0; i < halo; i++) 
                value.row(i) = value.row(patch_rows + i - 2 * halo - offset_x);
        } 
        if (bnd[1] == NEUMANN){
            for (int i = 0; i < halo; i++) 
                value.row(patch_rows - 1 - i) = value.row(patch_rows - 1 - halo);
        } else if (bnd[1] == PERIODIC){
            for (int i = 0; i < halo; i++) 
                value.row(patch_rows - 1 - i) = value.row(2 * halo - 1 - i + offset_x);
        }
        if (bnd[2] == NEUMANN){
            for (int i = 0; i < halo; i++) value.col(i) = value.col(halo);
        } else if (bnd[2] == PERIODIC){
            for (int i = 0; i < halo; i++) 
                value.col(i) = value.col(patch_cols + i - 2 * halo - offset_y);
        } 
        if (bnd[3] == NEUMANN){
            for (int i = 0; i < halo; i++) 
                value.col(patch_cols - 1 - i) = value.col(patch_cols - 1 - halo);
        } else if (bnd[3] == PERIODIC){
            for (int i = 0; i < halo; i++) 
                value.col(patch_cols - 1 - i) = value.col(2 * halo - 1 - i + offset_y);
        } 
    }
};

#endif
