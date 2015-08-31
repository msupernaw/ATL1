/* 
 * File:   VectorDotProduct.hpp
 * Author: matthewsupernaw
 *
 * Created on May 20, 2015, 4:25 PM
 */

#ifndef VECTORDOTPRODUCT_HPP
#define	VECTORDOTPRODUCT_HPP

#include "VectorExpressionBase.hpp"
#include "ArrayTraits.hpp"

namespace atl {

    template <class LHS, class RHS>
    inline const typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type Dot(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type R_TYPE;

        R_TYPE sum;

        int lmin = a.IndexMin();
        int lmax = a.IndexMax();
        int rmin = b.IndexMin();
        int rmax = b.IndexMax();

        int min = std::max(lmin, rmin);
        int max = std::min(lmax, rmax);

        for (int i = min; i <= max; i++) {
            sum += a(i) * b(i);
        }
        return sum;
    }



}

#endif	/* VECTORDOTPRODUCT_HPP */

