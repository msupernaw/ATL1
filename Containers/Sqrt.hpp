/* 
 * File:   Sqrt.hpp
 * Author: matthewsupernaw
 *
 * Created on May 28, 2015, 7:12 AM
 */

#ifndef SQRT_HPP
#define	SQRT_HPP

#include "VectorExpressionBase.hpp"


namespace atl {

    template< class EXP>
    struct VectorSqrt : VectorExpression<typename EXP::RET_TYPE, VectorSqrt<EXP> > {
        const EXP& expr_m;
        typedef typename EXP::RET_TYPE RET_TYPE;
        typedef typename EXP::BASE_TYPE BASE_TYPE;

        VectorSqrt(const EXP& exp) : expr_m(exp.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            return expr_m.Size(dimension);
        }

        //        inline const size_t Dimensions() const {
        //            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        //        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return std::sqrt(expr_m(i));
        }

        /*
         *
         *Returns the first valid index for this vector. 
         */
        inline const size_t IndexMin() const {
            return expr_m.IndexMin();
        }

        /*
         *
         *Returns the last valid index for this vector. 
         */
        inline const size_t IndexMax() const {
            return expr_m.IndexMax();
        }

        /**
         * Get a value based on the raw index for the underlying 
         * data. valid index is 0 - (length -1).
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return std::sqrt(expr_m(i));
        }

        inline void IsAliased(bool& aliased, void* ptr) const {
            expr_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            expr_m.ExpressionLength(length);
        }


    };

    template<typename EXP>
    const VectorSqrt< EXP> sqrt(const VectorExpression<typename EXP::RET_TYPE, EXP>& exp) {
        return VectorSqrt<EXP> (exp.Cast());
    }

}

#endif	/* SQRT_HPP */

