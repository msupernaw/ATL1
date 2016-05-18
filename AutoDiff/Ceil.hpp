/* 
 * File:   Ceil.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:56 PM
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


#ifndef ET4AD_CEIL_HPP
#define ET4AD_CEIL_HPP

#include <cmath>
#include "Expression.hpp"


namespace atl {

    /**
     * Expression template for handling the ceiling of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Ceil;
}

namespace std {

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Ceil<REAL_T, EXPR> ceil(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template for handling the ceiling of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Ceil : public ExpressionBase<REAL_T, Ceil<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Ceil(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::ceil(expr_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            expr_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr_m.PushIds(ids);
        }

         bool IsNonFunction()const {
             return expr_m.IsNonFunction();
        }
        
        bool IsNonlinear()const {
            return expr_m.IsNonLinear();
        }

        inline void MakeNLInteractions(bool b = false)const {
            expr_m.MakeNLInteractions(b);
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return 0.0; //expr_m.EvaluateDerivative(id) * GetValue();
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return 0.0; //expr_m.EvaluateDerivative(a, b) * GetValue();

        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return 0.0;
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicCeil<REAL_T>(expr_m.GetDynamicExpession());
        }


    private:
        const EXPR& expr_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Ceil<REAL_T, EXPR> ceil(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Ceil<REAL_T, EXPR > (expr.Cast());
    }


}

namespace std {

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Ceil<REAL_T, EXPR> ceil(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Ceil<REAL_T, EXPR > (expr.Cast());
    }
}


#endif /* CEIL_HPP */

