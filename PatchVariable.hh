#ifndef PROGVARIABLE
#define PROGVARIABLE

template <template <int>class gridtype, int halo = 1, int ntile_x = 1, int ntile_y = 1>
struct PatchVariable{
    ArrayXXf value;
    ArrayXXf value_t;
    gridtype<halo> tile[ntile_x][ntile_y];
    int patch_rows, patch_cols;
    int tile_rows, tile_cols;
    int offset_x, offset_y;

    PatchVariable(){};
    PatchVariable(int _nrows, int _ncols){
        int tx, ty;
        value.setRandom(_nrows, _ncols);
        value_t.setZero(_nrows, _ncols);
        offset_x = tile[0][0].offset_x;
        offset_y = tile[0][0].offset_y;
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
                        tile_rows, tile_cols
                        );
        patch_rows = _nrows;
        patch_cols = _ncols;
    }
    void update(float dt){
        if (offset_x == 1){
            for (int i = 1; i < ntile_x; i++)
                value_t.row(i * tile_rows + halo - i * offset_x) /= 2;
        }
        if (offset_y == 1){
            for (int i = 1; i < ntile_y; i++)
                value_t.col(i * tile_cols + halo - i * offset_y) /= 2;
        }
        value = value + value_t * dt;
        value_t.setZero();
    }
    inline void clean_t(){
        value_t.setZero();
    }
    inline void setLeftRightZero(){
        value.row(halo).setZero();
        value.row(patch_rows -1 - halo).setZero();
    }
    inline void setBottomTopZero(){
        value.col(halo).setZero();
        value.col(patch_cols -1 - halo).setZero();
    }
};

#endif
