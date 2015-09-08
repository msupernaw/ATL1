/* 
 * File:   Expression.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:16 AM
 */

#ifndef EXPRESSION_HPP
#define	EXPRESSION_HPP

//#define ET4AD_USE_SSE


#include "GradientStructure.hpp"
#include <string>
#include <sstream>
#include <memory>
#include <complex>

#ifdef ET4AD_USE_SSE

#include <x86intrin.h>
#include "Traits.hpp"

#endif

#define OVERLOAD_STD

namespace atl {

    /**
     * Base class for expression types.
     */
    template<class REAL_T, class A>
    struct ExpressionBase {
        typedef REAL_T BASE_TYPE;
        typedef A DIFF_EXPRESSION;

        ExpressionBase() {

        }

        /**
         * Cast this expression template to it's child
         * representation.
         * 
         * @return 
         */
        inline const A & Cast() const {
            return static_cast<const A&> (*this);
        }

        inline void VariableCount(uint32_t& count) const {
            Cast().VariableCount(count);
        }

        inline const REAL_T GetValue() const {
            return Cast().GetValue();
        }
        //
//
//        operator REAL_T() {
//            return Cast().GetValue();
//        }
//
//        operator REAL_T()const {
//            return Cast().GetValue();
//        }

      

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            Cast().PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t a) const {
            return Cast().EvaluateDerivative(a);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return Cast().EvaluateDerivative(a, b);
        }

        const ExpressionBase& operator=(const ExpressionBase & exp) const {
            return *this;
        }
    };

    template <class REAL_T, class T>
    std::ostream & operator<<(std::ostream &vout, const atl::ExpressionBase<REAL_T, T> &e) {
        vout << e.GetValue();
        return vout;
    }

    template<class REAL_T, class T, class TT>
    inline const int operator==(const atl::ExpressionBase<REAL_T, T>& lhs, const atl::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() == rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator!=(const atl::ExpressionBase<REAL_T, T>& lhs, const atl::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() != rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator<(const atl::ExpressionBase<REAL_T, T>& lhs, const atl::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() < rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator>(const atl::ExpressionBase<REAL_T, T>& lhs, const atl::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() > rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator<=(const atl::ExpressionBase<REAL_T, T>& lhs, const atl::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() <= rhs.GetValue();
    }

    template<class REAL_T, class T, class TT>
    inline const int operator>=(const atl::ExpressionBase<REAL_T, T>& lhs, const atl::ExpressionBase<REAL_T, TT>& rhs) {
        return lhs.GetValue() >= rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator==(const REAL_T &lhs, const atl::ExpressionBase<REAL_T, T>& rhs) {
        return lhs == rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator!=(const REAL_T &lhs, const atl::ExpressionBase<REAL_T, T>& rhs) {
        return lhs != rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator<(const REAL_T &lhs, const atl::ExpressionBase<REAL_T, T>& rhs) {
        return lhs < rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator>(const REAL_T &lhs, const atl::ExpressionBase<REAL_T, T>& rhs) {
        return lhs > rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator<=(const REAL_T &lhs, const atl::ExpressionBase<REAL_T, T>& rhs) {
        return lhs <= rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator>=(const REAL_T &lhs, const atl::ExpressionBase<REAL_T, T>& rhs) {
        return lhs >= rhs.GetValue();
    }

    template<class REAL_T, class T>
    inline const int operator==(const atl::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() == rhs;
    }

    template<class REAL_T, class T>
    inline const int operator!=(const atl::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() != rhs;
    }

    template<class REAL_T, class T>
    inline const int operator<(const atl::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() <= rhs;
    }

    template<class REAL_T, class T>
    inline const int operator>(const atl::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() > rhs;
    }

    template<class REAL_T, class T>
    inline const int operator<=(const atl::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() <= rhs;
    }

    template<class REAL_T, class T>
    inline const int operator>=(const atl::ExpressionBase<REAL_T, T>& lhs, const REAL_T &rhs) {
        return lhs.GetValue() >= rhs;
    }
}
namespace std {

    template<class T, class A>
    bool isfinite(atl::ExpressionBase<T, A>& v) {
        return false;
    }

    template<class T, class A>
    bool isinf(atl::ExpressionBase<T, A>& v) {
        return false;
    }
}

#endif	/* EXPRESSION_HPP */

