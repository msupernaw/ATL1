/* 
 * File:   ArrayExpressionBase.hpp
 * Author: matthewsupernaw
 *
 * Created on October 21, 2014, 9:57 AM
 */

#ifndef ARRAYEXPRESSIONBASE_HPP
#define	ARRAYEXPRESSIONBASE_HPP

#include "../Traits/Promote.hpp"
#include "../AutoDiff/AutoDiff.hpp"
namespace atl {

    template< class T, class A>
    struct ArrayExpression {
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
        //

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return Cast().operator()(i);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return Cast().operator()(i, j);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return Cast().operator()(i, j, k);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return Cast().operator()(i, j, k, l);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return Cast().operator()(i, j, k, l, m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return Cast().operator()(i, j, k, l, m, n);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return Cast().operator()(i, j, k, l, m, n, o);
        }

        inline void ExpressionLength(uint32_t& length)const {
            Cast().ExpressionLength(length);
        }

    };


}


#endif	/* ARRAYEXPRESSIONBASE_HPP */

