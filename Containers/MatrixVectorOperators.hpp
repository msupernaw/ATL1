/* 
 * File:   MatrixVectorOperators.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 3:32 PM
 */

#ifndef MATRIXVECTOROPERATORS_HPP
#define	MATRIXVECTOROPERATORS_HPP
#include "MatrixExpressionBase.hpp"

namespace atl {

    template< class T, class A>
    struct MatrixVectorExpression {
        typedef T RET_TYPE;

        const A & Cast() const {
            return static_cast<const A&> (*this);
        }

        inline const size_t Size(const int32_t & dimension) const {
            return Cast().Size(dimension);
        }

        inline const size_t Dimensions() const {
            return Cast(). Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return Cast().operator()(i, j);
        }

    };

    template< class LHS, class RHS>
    struct MatrixVectorAdd : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixVectorAdd<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixVectorAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(1) == rhs_m.Size(0));

#endif

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() > rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) + rhs_m(i);
        }




    };

    template< class LHS, class RHS>
    struct VectorMatrixAdd : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixAdd<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(0) == rhs_m.Size(1));

#endif

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() > rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i) + rhs_m(i, j);
        }

    };

    template< class LHS, class RHS>
    struct MatrixVectorSubtract : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixVectorSubtract<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixVectorSubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(1) == rhs_m.Size(0));

#endif

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() > rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) - rhs_m(i);
        }




    };

    template< class LHS, class RHS>
    struct VectorMatrixSubtract : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixSubtract<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixSubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(0) == rhs_m.Size(1));

#endif

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() > rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i) - rhs_m(i, j);
        }

    };

    template< class LHS, class RHS>
    struct MatrixVectorMultiply : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixVectorMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m; //matrix
        const RHS& rhs_m; //vector

        inline explicit MatrixVectorMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(1) == rhs_m.Size(0));

#endif

        }

        inline const size_t Dimensions() const {
            return 1;
        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return lhs_m.Size(0) < rhs_m.Size(0) ? lhs_m.Size(0) : rhs_m.Size(0);
                case 1:
                    return 1;
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret;
            for (int k = 0; k < rhs_m.Size(0); k++) {
                ret += rhs_m(k) * lhs_m(i, k);
            }
            return ret;
        }




    };

    template< class LHS, class RHS>
    struct VectorMatrixMultiply : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(0) == rhs_m.Size(1));

#endif

        }

        inline const size_t Dimensions() const {
            return 1;
        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return rhs_m.Size(0) < lhs_m.Size(0) ? rhs_m.Size(0) : lhs_m.Size(0);
                case 1:
                    return 1;
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret;
            for (int k = 0; k < rhs_m.Size(1); k++) {
                ret += lhs_m(k) * rhs_m(i, k);
            }
            return ret;
        }

    };

    template< class LHS, class RHS>
    struct MatrixVectorDivide : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixVectorDivide<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixVectorDivide(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(1) == rhs_m.Size(0));

#endif

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() > rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) / rhs_m(i);
        }




    };

    template< class LHS, class RHS>
    struct VectorMatrixDivide : MatrixVectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixDivide<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixDivide(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Size(0) == rhs_m.Size(1));

#endif

        }
        
        
        
        

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() > rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i) / rhs_m(i, j);
        }

    };

    template <class LHS, class RHS>
    inline const MatrixVectorAdd< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixVectorAdd< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorMatrixAdd< LHS, RHS> operator+(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorMatrixAdd< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixVectorSubtract< LHS, RHS> operator-(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixVectorSubtract< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorMatrixSubtract< LHS, RHS> operator-(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorMatrixSubtract< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixVectorMultiply< LHS, RHS> operator*(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixVectorMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorMatrixMultiply< LHS, RHS> operator*(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorMatrixMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixVectorDivide< LHS, RHS> operator/(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixVectorDivide< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorMatrixDivide< LHS, RHS> operator/(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorMatrixDivide< LHS, RHS > (a.Cast(), b.Cast());
    }

    template<class T2, class A>
    std::ostream& operator<<(std::ostream& out, const atl::MatrixVectorExpression< T2, A> &expr) {

        int dims = expr.Dimensions();
        std::stringstream ss;

        out << ss.str();
        for (int i = 0; i < expr.Size(0); i++) {
            for (int j = 0; j < expr.Size(1); j++) {

                out << expr(i, j) << " ";
            }
            out << "\n";
        }


        return out;
    }



}

#endif	/* MATRIXVECTOROPERATORS_HPP */

