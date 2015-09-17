/* 
 * File:   VectorArrayOperators.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 3:32 PM
 */

#ifndef VECTORARRAYOPERATORS_HPP
#define	VECTORARRAYOPERATORS_HPP

#include "VectorExpressionBase.hpp"
#include "ArrayExpressionBase.hpp"

namespace atl {
#ifdef ATL_CONTAINER_EXPERIMENTAL
    template< class T, class A>
    struct VectorArrayExpression {
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

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return Cast().operator()(i);
        }

    };

    template< class LHS, class RHS>
    struct VectorArrayAdd : VectorArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorArrayAdd<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorArrayAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            //            for (int i = 0; i < 7; i++) {
            assert(lhs_m.Dimensions() == rhs_m.Dimensions());
            assert(lhs_m.Size(0) == rhs_m.Size(0));
            //            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return lhs_m(i) + rhs_m(i);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }


    };

    template< class LHS, class RHS>
    struct VectorArraySubtract : VectorArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorArraySubtract<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorArraySubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING

            assert(lhs_m.Dimensions() == rhs_m.Dimensions());
            assert(lhs_m.Size(0) == rhs_m.Size(0));
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return lhs_m(i) - rhs_m(i);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }



    };

    template< class LHS, class RHS>
    struct VectorArrayMultiply : VectorArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorArrayMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorArrayMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            //            for (int i = 0; i < 7; i++) {
            //                assert(lhs_m.Size(i) == rhs_m.Size(i));
            //            }
            assert(lhs_m.Dimensions() == rhs_m.Dimensions());
            assert(lhs_m.Size(0) == rhs_m.Size(0));
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return lhs_m(i) * rhs_m(i);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }
    };

    template< class LHS, class RHS>
    struct VectorArrayDivide : VectorArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorArrayDivide<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorArrayDivide(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            //            for (int i = 0; i < 7; i++) {
            //                assert(lhs_m.Size(i) == rhs_m.Size(i));
            //            }
            assert(lhs_m.Dimensions() == rhs_m.Dimensions());
            assert(lhs_m.Size(0) == rhs_m.Size(0));
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return lhs_m(i) / rhs_m(i);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }
    };

    template <class LHS, class RHS>
    inline const VectorArrayAdd< LHS, RHS> operator+(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArrayAdd< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArrayAdd< LHS, RHS> operator+(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArrayAdd< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArraySubtract< LHS, RHS> operator-(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArraySubtract< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArraySubtract< LHS, RHS> operator-(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArraySubtract< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArrayMultiply< LHS, RHS> operator*(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArrayMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArrayMultiply< LHS, RHS> operator*(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArrayMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArrayDivide< LHS, RHS> operator/(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArrayDivide< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const VectorArrayDivide< LHS, RHS> operator/(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorArrayDivide< LHS, RHS > (a.Cast(), b.Cast());
    }
#endif
}


#endif	/* VECTORARRAYOPERATORS_HPP */

