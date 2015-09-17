/* 
 * File:   MatrixVectorOperators.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 3:32 PM
 */

#ifndef MATRIXVECTOROPERATORS_HPP
#define	MATRIXVECTOROPERATORS_HPP
#include "MatrixExpressionBase.hpp"
#include "../Traits/Primitive.hpp"
#include "Vector.hpp"
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
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef BASE_TYPE RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixVectorAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixVectorSubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixSubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m; //matrix
        const RHS& rhs_m; //vector

        inline explicit MatrixVectorMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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
            RET_TYPE ret = static_cast<RET_TYPE> (0.0);
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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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
            RET_TYPE ret = static_cast<RET_TYPE> (0.0);
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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixVectorDivide(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixDivide(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

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

    template< class LHS, class RHS>
    struct VectorMatrixVectorAdd : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorAdd<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixVectorAdd(const LHS& lhs, const RHS &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i) + rhs_m(i, j);
            return ret;
        }

    };

    template<class LHS, class RHS>
    const VectorMatrixVectorAdd<LHS, RHS> operator+(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
        return VectorMatrixVectorAdd<LHS, RHS>(lhs.Cast(), rhs.Cast());
    }

    template< class LHS, class RHS>
    struct VectorMatrixVectorSubtract : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorSubtract<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixVectorSubtract(const LHS& lhs, const RHS &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i) - rhs_m(i, j);
            return ret;
        }

    };

    template<class LHS, class RHS>
    const VectorMatrixVectorSubtract<LHS, RHS> operator-(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
        return VectorMatrixVectorSubtract<LHS, RHS>(lhs.Cast(), rhs.Cast());
    }
    
    template< class LHS, class RHS>
    struct VectorMatrixVectorMultiply : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixVectorMultiply(const LHS& lhs, const RHS &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i) * rhs_m(i, j);
            return ret;
        }

    };

    template<class LHS, class RHS>
    const VectorMatrixVectorMultiply<LHS, RHS> operator*(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
        return VectorMatrixVectorMultiply<LHS, RHS>(lhs.Cast(), rhs.Cast());
    }
    
    template< class LHS, class RHS>
    struct VectorMatrixVectorDivide : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorDivide<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixVectorDivide(const LHS& lhs, const RHS &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i) / rhs_m(i, j);
            return ret;
        }

    };

    template<class LHS, class RHS>
    const VectorMatrixVectorDivide<LHS, RHS> operator/(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
        return VectorMatrixVectorDivide<LHS, RHS>(lhs.Cast(), rhs.Cast());
    }
    

   
   template< class LHS, class RHS>
    struct VectorMatrixAddVector : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixAddVector<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixAddVector(const LHS& lhs, const RHS &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i,j) + rhs_m(i);
            return ret;
        }

    };

    template<class LHS, class RHS>
    const VectorMatrixAddVector<LHS, RHS> operator+(const atl::MatrixVectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::VectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
        return VectorMatrixAddVector<LHS, RHS>(lhs.Cast(), rhs.Cast());
    }
//
//    template< class LHS, class RHS>
//    struct VectorMatrixVectorSubtract : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorSubtract<LHS, RHS> > {
//        typedef typename LHS::RET_TYPE RET_TYPEL;
//        typedef typename RHS::RET_TYPE RET_TYPER;
//
//        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
//        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
//
//        const LHS& lhs_m;
//        const RHS& rhs_m;
//
//        inline explicit VectorMatrixVectorSubtract(const LHS& lhs, const RHS &rhs)
//        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
//
//        }
//
//        inline const size_t Dimensions() const {
//            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
//        }
//
//        inline const size_t Size(const int32_t & dimension) const {
//            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
//        }
//
//        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
//            RET_TYPE ret = lhs_m(i) - rhs_m(i, j);
//            return ret;
//        }
//
//    };
//
//    template<class LHS, class RHS>
//    const VectorMatrixVectorSubtract<LHS, RHS> operator-(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
//        return VectorMatrixVectorSubtract<LHS, RHS>(lhs.Cast(), rhs.Cast());
//    }
//    
//    template< class LHS, class RHS>
//    struct VectorMatrixVectorMultiply : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorMultiply<LHS, RHS> > {
//        typedef typename LHS::RET_TYPE RET_TYPEL;
//        typedef typename RHS::RET_TYPE RET_TYPER;
//
//        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
//        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
//
//        const LHS& lhs_m;
//        const RHS& rhs_m;
//
//        inline explicit VectorMatrixVectorMultiply(const LHS& lhs, const RHS &rhs)
//        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
//
//        }
//
//        inline const size_t Dimensions() const {
//            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
//        }
//
//        inline const size_t Size(const int32_t & dimension) const {
//            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
//        }
//
//        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
//            RET_TYPE ret = lhs_m(i) * rhs_m(i, j);
//            return ret;
//        }
//
//    };
//
//    template<class LHS, class RHS>
//    const VectorMatrixVectorMultiply<LHS, RHS> operator*(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
//        return VectorMatrixVectorMultiply<LHS, RHS>(lhs.Cast(), rhs.Cast());
//    }
//    
//    template< class LHS, class RHS>
//    struct VectorMatrixVectorDivide : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorMatrixVectorDivide<LHS, RHS> > {
//        typedef typename LHS::RET_TYPE RET_TYPEL;
//        typedef typename RHS::RET_TYPE RET_TYPER;
//
//        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
//        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
//
//        const LHS& lhs_m;
//        const RHS& rhs_m;
//
//        inline explicit VectorMatrixVectorDivide(const LHS& lhs, const RHS &rhs)
//        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
//
//        }
//
//        inline const size_t Dimensions() const {
//            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
//        }
//
//        inline const size_t Size(const int32_t & dimension) const {
//            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
//        }
//
//        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
//            RET_TYPE ret = lhs_m(i) / rhs_m(i, j);
//            return ret;
//        }
//
//    };
//
//    template<class LHS, class RHS>
//    const VectorMatrixVectorDivide<LHS, RHS> operator/(const atl::VectorExpression< typename LHS::RET_TYPE, LHS>& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) {
//        return VectorMatrixVectorDivide<LHS, RHS>(lhs.Cast(), rhs.Cast());
//    }
//     
   
    
    template< class T, class RHS>
    struct VectorMatrixScalarAdd : MatrixVectorExpression<typename atl::PromoteType< T, typename RHS::RET_TYPE >::return_type, VectorMatrixScalarAdd<T, RHS> > {
        typedef T RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        T lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixScalarAdd(T lhs, const RHS &rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m + rhs_m(i, j);
            return ret;
        }

    };

#define MAKE_SCALAR_VM_PLUS(TYPE) \
    template< class RHS> \
    const VectorMatrixScalarAdd<TYPE, RHS> operator+(const TYPE& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) { \
        return VectorMatrixScalarAdd<TYPE, RHS>(lhs, rhs.Cast()); \
    } \
    

    MAKE_SCALAR_VM_PLUS(atl::Variable<long double>)
    MAKE_SCALAR_VM_PLUS(atl::Variable<double>)
    MAKE_SCALAR_VM_PLUS(atl::Variable<float>)
    MAKE_SCALAR_VM_PLUS(long double)
    MAKE_SCALAR_VM_PLUS(double)
    MAKE_SCALAR_VM_PLUS(float)
    MAKE_SCALAR_VM_PLUS(int)
    MAKE_SCALAR_VM_PLUS(uint32_t)
    MAKE_SCALAR_VM_PLUS(uint64_t)
    MAKE_SCALAR_VM_PLUS(long)
    MAKE_SCALAR_VM_PLUS(short)

    template< class T, class RHS>
    struct VectorMatrixScalarSubtract : MatrixVectorExpression<typename atl::PromoteType< T, typename RHS::RET_TYPE >::return_type, VectorMatrixScalarSubtract<T, RHS> > {
        typedef T RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        T lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixScalarSubtract(T lhs, const RHS &rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m - rhs_m(i, j);
            return ret;
        }

    };

#define MAKE_VM_SCALAR_MINUS(TYPE) \
    template<class RHS> \
    const VectorMatrixScalarSubtract<TYPE, RHS> operator-(const TYPE& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) { \
        return VectorMatrixScalarSubtract<TYPE, RHS>(lhs, rhs.Cast()); \
    }\
    

    MAKE_VM_SCALAR_MINUS(atl::Variable<long double>)
    MAKE_VM_SCALAR_MINUS(atl::Variable<double>)
    MAKE_VM_SCALAR_MINUS(atl::Variable<float>)
    MAKE_VM_SCALAR_MINUS(long double)
    MAKE_VM_SCALAR_MINUS(double)
    MAKE_VM_SCALAR_MINUS(float)
    MAKE_VM_SCALAR_MINUS(uint32_t)
    MAKE_VM_SCALAR_MINUS(int)
    MAKE_VM_SCALAR_MINUS(uint64_t)
    MAKE_VM_SCALAR_MINUS(long)
    MAKE_VM_SCALAR_MINUS(short)

    template< class T, class RHS>
    struct VectorMatrixScalarMultiply : MatrixVectorExpression<typename atl::PromoteType< T, typename RHS::RET_TYPE >::return_type, VectorMatrixScalarMultiply<T, RHS> > {
        typedef T RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        T lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixScalarMultiply(T lhs, const RHS &rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m * rhs_m(i, j);
            return ret;
        }

    };

#define MAKE_VM_SCALAR_MULTIPLY(TYPE) \
    template< class RHS> \
    const VectorMatrixScalarMultiply<TYPE, RHS> operator*(const TYPE& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) { \
        return VectorMatrixScalarMultiply<TYPE, RHS>(lhs, rhs.Cast()); \
    }\
    

    MAKE_VM_SCALAR_MULTIPLY(atl::Variable<long double>)
    MAKE_VM_SCALAR_MULTIPLY(atl::Variable<double>)
    MAKE_VM_SCALAR_MULTIPLY(atl::Variable<float>)
    MAKE_VM_SCALAR_MULTIPLY(long double)
    MAKE_VM_SCALAR_MULTIPLY(double)
    MAKE_VM_SCALAR_MULTIPLY(float)
    MAKE_VM_SCALAR_MULTIPLY(uint32_t)
    MAKE_VM_SCALAR_MULTIPLY(int)
    MAKE_VM_SCALAR_MULTIPLY(uint64_t)
    MAKE_VM_SCALAR_MULTIPLY(long)
    MAKE_VM_SCALAR_MULTIPLY(short)

    template< class T, class RHS>
    struct VectorMatrixScalarDivide : MatrixVectorExpression<typename atl::PromoteType< T, typename RHS::RET_TYPE >::return_type, VectorMatrixScalarDivide<T, RHS> > {
        typedef T RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        T lhs_m;
        const RHS& rhs_m;

        inline explicit VectorMatrixScalarDivide(T lhs, const RHS &rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()) {

        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m / rhs_m(i, j);
            return ret;
        }

    };

#define MAKE_VM_SCALAR_DIVIDE(TYPE) \
    template< class RHS>\
    const VectorMatrixScalarDivide<TYPE, RHS> operator/(const TYPE& lhs, const atl::MatrixVectorExpression< typename RHS::RET_TYPE, RHS> &rhs) { \
        return VectorMatrixScalarDivide<TYPE, RHS>(lhs, rhs.Cast()); \
    } \


    MAKE_VM_SCALAR_DIVIDE(atl::Variable<long double>)
    MAKE_VM_SCALAR_DIVIDE(atl::Variable<double>)
    MAKE_VM_SCALAR_DIVIDE(atl::Variable<float>)
    MAKE_VM_SCALAR_DIVIDE(long double)
    MAKE_VM_SCALAR_DIVIDE(double)
    MAKE_VM_SCALAR_DIVIDE(float)
    MAKE_VM_SCALAR_DIVIDE(uint32_t)
    MAKE_VM_SCALAR_DIVIDE(int)
    MAKE_VM_SCALAR_DIVIDE(uint64_t)
    MAKE_VM_SCALAR_DIVIDE(long)
    MAKE_VM_SCALAR_DIVIDE(short)

    template< class LHS, class T>
    struct VectorMatrixAddScalar : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, T >::return_type, VectorMatrixAddScalar<LHS, T> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef T RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T >::return_type RET_TYPE;

        const LHS& lhs_m;
        T rhs_m;

        inline explicit VectorMatrixAddScalar(const LHS& lhs, const T &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i, j) + rhs_m;
            return ret;
        }

    };

#define MAKE_VM_ADD_SCALAR(TYPE) \
    template<class LHS> \
    const VectorMatrixAddScalar<LHS, TYPE> operator+(const atl::MatrixVectorExpression< typename LHS::RET_TYPE, LHS> &lhs, const TYPE& rhs) { \
        return VectorMatrixAddScalar<LHS, TYPE>(lhs.Cast(), rhs);\
    }\
    

    MAKE_VM_ADD_SCALAR(atl::Variable<long double>)
    MAKE_VM_ADD_SCALAR(atl::Variable<double>)
    MAKE_VM_ADD_SCALAR(atl::Variable<float>)
    MAKE_VM_ADD_SCALAR(long double)
    MAKE_VM_ADD_SCALAR(double)
    MAKE_VM_ADD_SCALAR(float)
    MAKE_VM_ADD_SCALAR(uint32_t)
    MAKE_VM_ADD_SCALAR(int)
    MAKE_VM_ADD_SCALAR(uint64_t)
    MAKE_VM_ADD_SCALAR(long)
    MAKE_VM_ADD_SCALAR(short)

    template< class LHS, class T>
    struct VectorMatrixSubtractScalar : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, T >::return_type, VectorMatrixSubtractScalar<LHS, T> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef T RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T >::return_type RET_TYPE;

        const LHS& lhs_m;
        T rhs_m;

        inline explicit VectorMatrixSubtractScalar(const LHS& lhs, const T &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i, j) - rhs_m;
            return ret;
        }

    };

#define MAKE_VM_SUBTRACT_SCALAR(TYPE)\
    template<class LHS>\
    const VectorMatrixSubtractScalar<LHS, TYPE> operator-(const atl::MatrixVectorExpression< typename LHS::RET_TYPE, LHS> &lhs, const TYPE& rhs) {\
        return VectorMatrixSubtractScalar<LHS, TYPE>(lhs.Cast(), rhs);\
    }\


    MAKE_VM_SUBTRACT_SCALAR(atl::Variable<long double>)
    MAKE_VM_SUBTRACT_SCALAR(atl::Variable<double>)
    MAKE_VM_SUBTRACT_SCALAR(atl::Variable<float>)
    MAKE_VM_SUBTRACT_SCALAR(long double)
    MAKE_VM_SUBTRACT_SCALAR(double)
    MAKE_VM_SUBTRACT_SCALAR(float)
    MAKE_VM_SUBTRACT_SCALAR(uint32_t)
    MAKE_VM_SUBTRACT_SCALAR(int)
    MAKE_VM_SUBTRACT_SCALAR(uint64_t)
    MAKE_VM_SUBTRACT_SCALAR(long)
    MAKE_VM_SUBTRACT_SCALAR(short)

    template< class LHS, class T>
    struct VectorMatrixMultiplyScalar : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, T >::return_type, VectorMatrixMultiplyScalar<LHS, T> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef T RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T >::return_type RET_TYPE;

        const LHS& lhs_m;
        T rhs_m;

        inline explicit VectorMatrixMultiplyScalar(const LHS& lhs, const T &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i, j) * rhs_m;
            return ret;
        }

    };

#define MAKE_VM_MULTIPLY_SCALAR(TYPE) \
    template<class LHS> \
    const VectorMatrixMultiplyScalar<LHS, TYPE> operator*(const atl::MatrixVectorExpression< typename LHS::RET_TYPE, LHS> &lhs, const TYPE& rhs) { \
        return VectorMatrixMultiplyScalar<LHS, TYPE>(lhs.Cast(), rhs); \
    } \


    MAKE_VM_MULTIPLY_SCALAR(atl::Variable<long double>)
    MAKE_VM_MULTIPLY_SCALAR(atl::Variable<double>)
    MAKE_VM_MULTIPLY_SCALAR(atl::Variable<float>)
    MAKE_VM_MULTIPLY_SCALAR(long double)
    MAKE_VM_MULTIPLY_SCALAR(double)
    MAKE_VM_MULTIPLY_SCALAR(float)
    MAKE_VM_MULTIPLY_SCALAR(uint32_t)
    MAKE_VM_MULTIPLY_SCALAR(int)
    MAKE_VM_MULTIPLY_SCALAR(uint64_t)
    MAKE_VM_MULTIPLY_SCALAR(long)
    MAKE_VM_MULTIPLY_SCALAR(short)

    template< class LHS, class T>
    struct VectorMatrixDivideScalar : MatrixVectorExpression<typename atl::PromoteType< typename LHS::RET_TYPE, T >::return_type, VectorMatrixDivideScalar<LHS, T> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef T RET_TYPER;

        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T >::return_type BASE_TYPE;
        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T >::return_type RET_TYPE;

        const LHS& lhs_m;
        T rhs_m;

        inline explicit VectorMatrixDivideScalar(const LHS& lhs, const T &rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {

        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i, j) / rhs_m;
            return ret;
        }

    };

#define MAKE_VM_DIVIDE_SCALAR(TYPE) \
    template<class LHS> \
    const VectorMatrixDivideScalar<LHS, TYPE> operator/(const atl::MatrixVectorExpression< typename LHS::RET_TYPE, LHS> &lhs, const TYPE& rhs) { \
        return VectorMatrixDivideScalar<LHS, TYPE>(lhs.Cast(), rhs); \
    } \

    MAKE_VM_DIVIDE_SCALAR(atl::Variable<long double>)
    MAKE_VM_DIVIDE_SCALAR(atl::Variable<double>)
    MAKE_VM_DIVIDE_SCALAR(atl::Variable<float>)
    MAKE_VM_DIVIDE_SCALAR(long double)
    MAKE_VM_DIVIDE_SCALAR(double)
    MAKE_VM_DIVIDE_SCALAR(float)
    MAKE_VM_DIVIDE_SCALAR(uint32_t)
    MAKE_VM_DIVIDE_SCALAR(int)
    MAKE_VM_DIVIDE_SCALAR(uint64_t)
    MAKE_VM_DIVIDE_SCALAR(long)
    MAKE_VM_DIVIDE_SCALAR(short)    

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

