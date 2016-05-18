/* 
 * File:   Expression.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:16 AM
 */


/**
 *
 * @author  Matthew R. Supernaw
 *
 * Public Domain Notice
 * National Oceanic And Atmospheric Administration
 *
 * This software is a "United States Government Work" under the terms of the
 * United States Copyright Act.  It was written as part of the author's official
 * duties as a United States Government employee and thus cannot be copyrighted.
 * This software is freely available to the public for use. The National Oceanic
 * And Atmospheric Administration and the U.S. Government have not placed any
 * restriction on its use or reproduction.  Although all reasonable efforts have
 * been taken to ensure the accuracy and reliability of the software and data,
 * the National Oceanic And Atmospheric Administration and the U.S. Government
 * do not and cannot warrant the performance warrant the performance or results
 * that may be obtained by using this  software or data. The National Oceanic
 * And Atmospheric Administration and the U.S. Government disclaim all
 * warranties, express or implied, including warranties of performance,
 * merchantability or fitness for any particular purpose.
 *
 * Please cite the author(s) in any work or product based on this material.
 *
 */


#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP

//#define ET4AD_USE_SSE



#include "GradientStructure.hpp"
#include "DynamicExpression.hpp"
#include "../Utilities/Platform.hpp"
#include "../Utilities/BigFloat.hpp"
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

        //        operator REAL_T() {
        //            return Cast().GetValue();
        //        }
        //
        //        operator REAL_T()const {
        //            return Cast().GetValue();
        //        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            Cast().PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            Cast().PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            Cast().PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            Cast().PushAdjoints(adjoints, coefficient);
        }

        bool IsNonlinear()const {
            return Cast().IsNonlinear();
        }

        bool IsNonFunction()const {
            return Cast().IsNonFunction();
        }

        inline void MakeNLInteractions(bool b = false)const {
            Cast().MakeNLInteractions(b);
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            Cast().PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t a) const {
            return Cast().EvaluateDerivative(a);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return Cast().EvaluateDerivative(a, b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return Cast().EvaluateDerivative(x, y, z);
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return Cast().GetDynamicExpession();
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

#endif /* EXPRESSION_HPP */

