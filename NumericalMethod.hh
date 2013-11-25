#ifndef NUMERICAL_METHOD
#define NUMERICAL_METHOD
#include <Eigen/Dense>
#include <cmath>
#include "Include.hh"

using namespace Eigen;

template <int order>
class FiniteDifference{
    public:
        template<typename Derived>
        ArrayXXf x(const ArrayBase<Derived> &var, float dx) const{
            ArrayXXf result;
            int nrows = var.rows(), ncols = var.cols();
            if (order == 2){
                result.resize(nrows, ncols);
                result.row(0) = (var.row(1) - var.row(0)) / dx;
                for (size_t i = 1; i < nrows - 1; i++)
                    result.row(i) = (var.row(i + 1) - var.row(i - 1)) / (2. * dx);
                result.row(nrows - 1) = (var.row(nrows - 1) - var.row(nrows - 2)) / dx;
            } else if (order == 1){
                result.resize(nrows - 1, ncols);
                for (size_t i = 0; i < nrows - 1; i++)
                    result.row(i) = (var.row(i + 1) - var.row(i)) / dx;
            }
            return result;
        }
        template<typename Derived>
        ArrayXXf y(const ArrayBase<Derived> &var, float dy) const{
            ArrayXXf result;
            int nrows = var.rows(), ncols = var.cols();
            if (order == 2){
                result.resize(nrows, ncols);
                result.col(0) = (var.col(1) - var.col(0)) / dy;
                for (size_t i = 1; i < ncols - 1; i++)
                    result.col(i) = (var.col(i + 1) - var.col(i - 1)) / (2. * dy);
                result.col(ncols - 1) = (var.col(ncols - 1) - var.col(ncols - 2)) / dy;
            } else if (order == 1){
                result.resize(nrows, ncols - 1);
                for (size_t i = 0; i < ncols - 1; i++)
                    result.col(i) = (var.col(i + 1) - var.col(i)) / dy;
            }
            return result;
        }
};

template <int order>
class HyperDifference{
    // calculate higher order difference
protected:
    static const float coeff[4][7];
public:
    template<typename Derived>
    ArrayXXf x(const ArrayBase<Derived> &var, float dx) const{
        int r, nrows = var.rows(), ncols = var.cols();
        ArrayXXf result = ArrayXXf::Zero(nrows, ncols);
        for (size_t i = 0; i < nrows; i++){
            r = MIN3(i, nrows - i - 1, order / 2);
            for (int j = -r; j < r + 1; j++)
                result.row(i) += coeff[r][j + 3] * var.row(i + j) * std::pow(dx, -order);
        }
        return result;
    }
    template<typename Derived>
    ArrayXXf y(const ArrayBase<Derived> &var, float dy) const{
        int r, nrows = var.rows(), ncols = var.cols();
        ArrayXXf result = ArrayXXf::Zero(nrows, ncols);
        for (size_t i = 0; i < ncols; i++){
            r = MIN3(i, ncols - i - 1, order / 2);
            for (int j = -r; j < r + 1; j++)
                result.col(i) += coeff[r][j + 3] * var.col(i + j) * std::pow(dy, -order);
        }
        return result;
    }
};

template <int  order>
const float HyperDifference<order>::coeff[4][7] = {
    {0,  0,  0,   0,  0,  0, 0},
    {0,  0,  1,  -2,  1,  0, 0},
    {0, -1,  4,  -6,  4, -1, 0},
    {1, -6, 15, -20, 15, -6, 1}
};


template <int order>
class FiniteInterpolation{
protected:
    ArrayXXf wsignx, wsigny;
    static const float coeff[6][6];
public:
    template<typename Derived>
    ArrayXXf x(const ArrayBase<Derived> &var) const{
        int r, nrows = var.rows(), ncols = var.cols();
        ArrayXXf result = ArrayXXf::Zero(nrows - 1, ncols);
        for (size_t i = 0; i < nrows - 1; i++){
            r = MIN3(i + 1, nrows - i - 1, order / 2);
            for (int j = -r + 1; j < r + 1; j++)
                result.row(i) += coeff[2 * r - 1][j + 2] * var.row(i + j);
        }
        return result;
    }
    template<typename Derived>
    ArrayXXf y(const ArrayBase<Derived> &var) const{
        int r, nrows = var.rows(), ncols = var.cols();
        ArrayXXf result = ArrayXXf::Zero(nrows, ncols - 1);
        for (size_t i = 0; i < ncols - 1; i++){
            r = MIN3(i + 1, ncols - i - 1, order / 2);
            for (int j = -r + 1; j < r + 1; j++)
                result.col(i) += coeff[2 * r - 1][j + 2] * var.col(i + j);
        }
        return result;
    }
    template<typename Deriveda, typename Derivedb>
    ArrayXXf x1(const ArrayBase<Deriveda> &wind, const ArrayBase<Derivedb> &var) const{
        // upwind transport
        int r, nrows = nrows, ncols = ncols;
        ArrayXXf result = ArrayXXf::Zero(nrows - 1, ncols);
        wsignx = (wind.block(1,0,nrows - 1,ncols) > 0).select(
                ArrayXXf::Zero(nrows - 1, ncols) + 1., 
                ArrayXXf::Zero(nrows - 1, ncols) - 1.);
        for (size_t i = 0; i < nrows - 1; i++){
            r = MIN3(i + 1, nrows - i - 1, order / 2);
            for (int j = -r + 1; j < r + 1; j++)
                result.row(i) += coeff[2 * r - 1][j + 2] * var.row(i + j) +
                    wsignx.row(i) * coeff[2 * r - 2][j + 2] * var.row(i + j);
        }
        return result;
    }
    template<typename Deriveda, typename Derivedb>
    ArrayXXf y1(const ArrayBase<Deriveda> &wind, const ArrayBase<Derivedb> &var) const{
        int r, nrows = nrows, ncols = ncols;
        ArrayXXf result = ArrayXXf::Zero(nrows, ncols - 1);
        wsigny = (wind.block(0,1,nrows,ncols - 1) > 0).select(
                ArrayXXf::Zero(nrows, ncols - 1) + 1., 
                ArrayXXf::Zero(nrows, ncols - 1) - 1.);
        for (size_t i = 0; i < ncols - 1; i++){
            r = MIN3(i + 1, ncols - i - 1, order / 2);
            for (int j = -r + 1; j < r + 1; j++)
                result.col(i) += coeff[2 * r - 1][j + 2] * var.col(i + j) +
                    wsigny.col(i) * coeff[2 * r - 2][j + 2] * var.col(i + j);
        }
        return result;
    }
};

template <int order>
const float FiniteInterpolation<order>::coeff[6][6] = {
        { 0   ,  0     , 0.5    , -0.5    ,  0      ,  0     },
        { 0   ,  0     , 0.5    ,  0.5    ,  0      ,  0     },
        { 0   , -1./12., 1./4.  , -1./4.  ,  1./12. ,  0     },
        { 0   , -1./12., 7./12. ,  7./12  , -1./12. ,  0     },
        {1./60, -1./12., 1./6.  , -1./6.  ,  1./12. , -1./60.},
        {1./60, -2./15., 37./60.,  37./60., -2./15. ,  1./60.}
};

#endif

