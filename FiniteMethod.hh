#define FINITEMETHOD
#define FINITEMETHOD

#include "Include.hh"
#include "AuxGrid.hh"
#define MIDDLE(order, var, i) ( \
        order == 1 ? - 0.5 * (var(i + 1) - var(i)) : ( \
        order == 2 ? \
            static_cast<MainGrid>(0.5 * (var(i + 1) + var(i))) : ( \
        order == 3 ? \
            static_cast<MainGrid>(1./12. * (var(i + 2) - var(i - 1)) \
                - 1./4. * (var(i + 1) - var(i))) : ( \
        order == 4 ? \
            static_cast<MainGrid>(7./12. * (var(i + 1) + var(i)) \
                - 1./12. * (var(i + 2) + var(i - 1))) : ( \
        order == 5 ? \
            static_cast<MainGrid>(-1./60. * (var(i + 3) - var(i - 2)) \
                + 1./12. * (var(i + 2) - var(i - 1)) \
                - 1./6. * (var(i + 1) - var(i))) : \
        /* order == 6 */ \
            static_cast<MainGrid>(37./60. * (var(i + 1) + var(i)) \
                - 2./15. * (var(i + 2) + var(i - 1)) \
                + 1./60. * (var(i + 3) + var(i - 2))))))) \
        )
#define DIFFN(order, var, i) ( \
        order == 2 ? var(i-1)-2*var(i)+var(i+1) : ( \
        order == 4 ? \
            static_cast<MainGrid>(-var(i-2)+4*var(i-1)-6*var(i)+4*var(i+1)-var(i+2)) : \
        /* order == 6 */ \
            static_cast<MainGrid>(var(i-3)-6*var(i-2)+15*var(i-1)-20*var(i)+15*var(i+1)-6*var(i+2)+var(i+3))) \
        )

template <int order>
class Difference{
    /* Calculate finite difference */
public:
    MainGrid x(const MainGrid &var){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        if (order == 2){
            result.resize(nrows, ncols);
            result.row(0) = var.row(1) - var.row(0);
            for (size_t i = 1; i < nrows - 1; i++)
                result.row(i) = 0.5 * (var.row(i + 1) - var.row(i - 1));
            result.row(nrows - 1) = var.row(nrows - 1) - var.row(nrows - 2);
        } else if (order == 1){
            result.resize(nrows - 1, ncols);
            for (size_t i = 0; i < nrows - 1; i++)
                result.row(i) = var.row(i + 1) - var.row(i);
        } 
        return result;
    }
    MainGrid y(const MainGrid &var){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        if (order == 2){
            result.resize(nrows, ncols);
            result.col(0) = var.col(1) - var.col(0);
            for (size_t i = 1; i < var.cols() - 1; i++)
                result.col(i) = 0.5 * (var.col(i + 1) - var.col(i - 1));
            result.col(ncols - 1) = var.col(ncols - 1) - var.col(ncols - 2);
        } else if (order == 1){
            result.resize(nrows, ncols - 1);
            for (size_t i = 0; i < ncols - 1; i++)
                result.col(i) = var.col(i + 1) - var.col(i);
        } 
        return result;
    }
    MainGrid x(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        if (order == 2){
            result.resize(nrows, ncols);
            result = this->x(var);
            result.row(0) = (var.row(1) - aux.left.row(0)) / 2.;
            result.row(nrows - 1) = (aux.right.row(0) - var.row(nrows - 2)) / 2.;
        } else if (order == 1){
            result.resize(nrows + 1, ncols);
            result.block(1, 0, nrows - 1, ncols) = this->x(var);
            result.row(0) = var.row(0) - aux.left.row(0);
            result.row(nrows) = aux.right.row(0) - var.row(nrows - 1);
        } 
        return result;
    }
    MainGrid y(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        if (order == 2){
            result.resize(nrows, ncols);
            result = this->y(var);
            result.col(0) = (var.col(1) - aux.bottom.col(0)) / 2.;
            result.col(ncols - 1) = (aux.top.col(0) - var.col(ncols - 2)) / 2.;
        } else if (order == 1){
            result.resize(nrows, ncols + 1);
            result.block(0, 1, nrows, ncols - 1) = this->y(var);
            result.col(0) = var.col(0) - aux.bottom.col(0);
            result.col(ncols) = aux.top.col(0) - var.col(ncols - 1);
        } 
        return result;
    }
};

template <int order>
class DifferenceN{
    /* calculate higher order difference */
public:
    MainGrid x(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows, ncols);
        result.row(0) = var.row(1) - 2 * var.row(0) + aux.left.row(0);
        for (size_t i = 1; i < nrows - 1; i++)
            result.row(i) = DIFFN(2*MIN(i, nrows-1-i, order/2), var.row, i);
        result.row(nrows - 1) = var.row(nrows - 2) - 2 * var.row(nrows - 1) + aux.right.row(0);
        return result;
    }
    MainGrid y(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows, ncols);
        result.col(0) = var.col(1) - 2 * var.col(0) + aux.bottom.col(0);
        for (size_t i = 1; i < ncols - 1; i++)
            result.col(i) = DIFFN(2*MIN(i, ncols-1-i, order/2), var.col, i);
        result.col(ncols - 1) = var.col(ncols - 2) - 2 * var.col(ncols - 1) + aux.top.col(0);
        return result;
    }
};

template <int order = 0>
class Integral{
    /* integrate over an axis, problem exist when you differential it on one axis
     * and integrate it over another axis */
public:
    MainGrid xup(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows - 1, ncols);
        result.row(0) = aux.left.row(0) + var.row(0);
        for (size_t i = 1; i < nrows - 1; i++)
            result.row(i) = result.row(i - 1) + var.row(i);
        return result;
    }
    MainGrid xdown(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows - 1, ncols);
        result.row(nrows - 2) = aux.right.row(0) - var.row(nrows - 1);
        for (int i = nrows - 2; i > 0; i--)
            result.row(i - 1) = result.row(i) - var.row(i);
        return result;
    }
    MainGrid yup(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        /*
        result.resize(nrows, ncols - 1);
        result.col(0) = aux.bottom.col(0) + var.col(0);
        for (size_t i = 1; i < ncols - 1; i++)
            result.col(i) = result.col(i - 1) + var.col(i);
        */
        result.resize(nrows, ncols + 1);
        result.col(0) = aux.bottom.col(0);
        for (size_t i = 0; i < ncols; i++)
            result.col(i + 1) = result.col(i) + var.col(i);
        return result;
    }
    MainGrid ydown(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        /*
        result.resize(nrows, ncols - 1);
        result.col(ncols - 2) = aux.top.col(0) - var.col(ncols - 1);
        for (int i = ncols - 2; i > 0; i--)
            result.col(i - 1) = result.col(i) - var.col(i);
        */
        result.resize(nrows, ncols + 1);
        result.col(ncols) = aux.top.col(0);
        for (int i = ncols; i > 0; i--)
            result.col(i - 1) = result.col(i) - var.col(i - 1);
        return result;
    }
};

template <int order>
class Interpolate{
    /* Make an interpolation to auxf grid */
protected:
    MainGrid wsignx, wsigny, buffer;
public:
    MainGrid x(const MainGrid &var){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows - 1, ncols);
        for (size_t i = 0; i < nrows - 1; i++)
            result.row(i) = MIDDLE(2*MIN(i+1,nrows-1-i,order/2), var.row, i);
        return result;
    }
    MainGrid y(const MainGrid &var){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows, ncols - 1);
        for (size_t i = 0; i < ncols - 1; i++)
            result.col(i) = MIDDLE(2*MIN(i+1,ncols-1-i,order/2), var.col, i);
        return result;
    }
    MainGrid x(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows + 1, ncols);
        result.block(1, 0, nrows - 1, ncols) = this->x(var);
        result.row(0) = (var.row(0) + aux.left.row(0)) / 2.;
        result.row(nrows) = (var.row(nrows - 1) + aux.right.row(0)) / 2.;
        return result;
    }
    MainGrid y(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows, ncols + 1);
        result.block(0, 1, nrows, ncols - 1) = this->y(var);
        result.col(0) = (var.col(0) + aux.bottom) / 2.;
        result.col(ncols) = (var.col(ncols - 1) + aux.top) / 2.;
        return result;
    }
    MainGrid quad(const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        buffer.resize(nrows + 2, ncols + 2);
        buffer.block(1, 1, nrows, ncols) = var;
        buffer.block(0, 1, 1, ncols) = aux.left;
        buffer.block(nrows + 1, 1, 1, ncols) = aux.right;
        buffer.block(1, 0, nrows, 1) = aux.bottom;
        buffer.block(1, ncols + 1, nrows, 1) = aux.top;
        buffer(0, 0) = (buffer(0, 1) + buffer(1, 0)) / 2.;
        buffer(0, ncols + 1) = (buffer(0, ncols) + buffer(1, ncols + 1)) / 2.;
        buffer(nrows + 1, 0) = (buffer(nrows, 0) + buffer(nrows + 1, 1)) / 2.;
        buffer(nrows + 1, ncols + 1) = (buffer(nrows + 1, ncols) + buffer(nrows, ncols + 1)) / 2.;
        result = this->y(this->x(buffer));
        return result;
    }
    MainGrid x1(const MainGrid &wind, const MainGrid &var){
        // upwind transport
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows - 1, ncols);
        wsignx = (wind.block(1,0,nrows - 1,ncols) > 0).select(
                ZERO2(nrows - 1, ncols) + 1., 
                ZERO2(nrows - 1, ncols) - 1.);
        for (size_t i = 0; i < nrows - 1; i++)
            result.row(i) = MIDDLE(2*MIN(i+1, nrows-1-i, order/2), var.row, i) +
                wsignx.row(i) * MIDDLE(2*MIN(i+1, nrows-1-i, order/2) - 1, var.row, i);
        return result;
    }
    MainGrid y1(const MainGrid &wind, const MainGrid &var){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows, ncols - 1);
        wsigny = (wind.block(0,1,nrows,ncols - 1) > 0).select(
                ZERO2(nrows, ncols - 1) + 1., 
                ZERO2(nrows, ncols - 1) - 1.);
        for (size_t i = 0; i < ncols - 1; i++)
            result.col(i) = MIDDLE(2*MIN(i+1, ncols-1-i, order/2), var.col, i) +
                wsigny.col(i) * MIDDLE(2*MIN(i+1, ncols-1-i, order/2) - 1, var.col, i);
        return result;
    }
    MainGrid x1(const MainGrid &wind, const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows + 1, ncols);
        result.block(1, 0, nrows - 1, ncols) = this->x1(wind, var);
        result.row(0) = (var.row(0) + aux.left) / 2. +
            (wind.row(0) > 0).select(ZERO2(1, ncols) + 1., ZERO2(1, ncols) - 1.)
            * 0.5 * (aux.left - var.row(0));
        result.row(nrows) = (var.row(nrows - 1) + aux.right) / 2. +
            (wind.row(nrows) > 0).select(ZERO2(1, ncols) + 1., ZERO2(1, ncols) - 1.)
            * 0.5 * (var.row(nrows - 1) - aux.right);
        return result;
    }
    MainGrid y1(const MainGrid &wind, const MainGrid &var, const AuxGrid &aux){
        int nrows = var.rows(), ncols = var.cols();
        MainGrid result;
        result.resize(nrows, ncols + 1);
        result.block(0, 1, nrows, ncols - 1) = this->y1(wind, var);
        result.col(0) = (var.col(0) + aux.bottom) / 2. + 
            (wind.col(0) > 0).select(ZERO2(nrows, 1) + 1., ZERO2(nrows, 1) - 1.) 
            * 0.5 * (aux.bottom - var.col(0));
        result.col(ncols) = (var.col(ncols - 1) + aux.top) / 2. + 
            (wind.col(ncols) > 0).select(ZERO2(nrows, 1) + 1., ZERO2(nrows, 1) - 1.) 
            * 0.5 * (var.col(ncols - 1) - aux.top);
        return result;
    }
};

#endif
