/* 
 * File:   MatrixExpressionBase.hpp
 * Author: matthewsupernaw
 *
 * Created on October 21, 2014, 10:01 AM
 */

#ifndef MATRIXEXPRESSIONBASE_HPP
#define	MATRIXEXPRESSIONBASE_HPP
#include "../Traits/Promote.hpp"
#include "../AutoDiff/AutoDiff.hpp"

namespace atl {

    template< class T, class A>
    struct MatrixExpression {
        typedef T RET_TYPE;

        const A & Cast() const {
            return static_cast<const A&> (*this);
        }

        inline const size_t Size(const int32_t & dimension) const {
            return Cast().Size(dimension);
        }

        inline const size_t Dimensions() const {
            return 2;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return Cast().operator()(i, j);
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return Cast().IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return Cast().IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return Cast().AtRaw(i, j);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            Cast().IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            Cast().ExpressionLength(length);
        }
    };



}

#endif	/* MATRIXEXPRESSIONBASE_HPP */

