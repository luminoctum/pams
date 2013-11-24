#ifndef PROGVARIABLE
#define PROGVARIABLE

struct ProgVariable{
    ArrayXXf value;
    ArrayXXf tendency;
    int nrows, ncols;

    ProgVariable(){};
    ProgVariable(int _nrows, int _ncols){
        nrows = _nrows;
        ncols = _ncols;
        value.setRandom(nrows, ncols);
        tendency.setZero(nrows, ncols);
    }
    inline void resize(int nrows, int ncols){
        value.resize(nrows, ncols);
        tendency.resize(nrows, ncols);
        value.setRandom();
        tendency.setZero();
    };
    inline void update(){
        value = value + tendency;
        tendency.setZero();
    }
};

#endif
