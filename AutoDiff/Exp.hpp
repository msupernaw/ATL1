/* 
 * File:   Exp.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:37 PM
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


#ifndef ET4AD_EXP_HPP
#define ET4AD_EXP_HPP

#include <cmath>
#include "Expression.hpp"
#define EXP_OF_B REAL_T(114200738981568423454048256.0)

namespace atl {
    /**
     * Expression template to compute e raised to a expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Exp;
}
namespace std {

    /**
     * Override for the exp function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Exp<REAL_T, EXPR> exp(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}
namespace atl {

    /**
     * Expression template to compute e raised to a expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Exp : public ExpressionBase<REAL_T, Exp<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;
        const EXPR& expr_m;
        const REAL_T value_m;
        const REAL_T evalue_m;

        Exp(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(Compute(expr.GetValue())),evalue_m(expr.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline const REAL_T Compute(const REAL_T &x) {
            return std::exp(x);
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

        bool IsNonFunction() const {
            return true;
        }

        bool IsNonlinear()const {
            return true;
        }

        inline void MakeNLInteractions(bool b = false)const {
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            //                        expr_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * std::log(M_E) * this->GetValue();
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            //            REAL_T fx = value_m;
            //            return value_m*(expr_m.EvaluateDerivative(a))*(expr_m.EvaluateDerivative(b)) + value_m*(expr_m.EvaluateDerivative(a, b));

            return value_m * std::pow(std::log(M_E), 2.0)*(expr_m.EvaluateDerivative(a))*(expr_m.EvaluateDerivative(b)) + value_m * std::log(M_E)*(expr_m.EvaluateDerivative(a, b));

        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            //            return myexp(evalue_m)*(expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))*(expr_m.EvaluateDerivative(z)) + myexp(evalue_m)*(expr_m.EvaluateDerivative(x, y))*(expr_m.EvaluateDerivative(z)) +
            //                    myexp(evalue_m)*(expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y, z)) + myexp(evalue_m)*(expr_m.EvaluateDerivative(x, z))*(expr_m.EvaluateDerivative(y)) + myexp(evalue_m)*
            //                    (expr_m.EvaluateDerivative(x, y, z));

            return value_m * std::pow(std::log(M_E), 3.0) * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))*(expr_m.EvaluateDerivative(z)) + value_m * std::pow(std::log(M_E), 2.0) * (expr_m.EvaluateDerivative(x, y))*
                    (expr_m.EvaluateDerivative(z)) + value_m * std::pow(std::log(M_E), 2.0) * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y, z)) + value_m * std::pow(std::log(M_E), 2.0)* (expr_m.EvaluateDerivative(x, z))*(expr_m.EvaluateDerivative(y))
                    + value_m * std::log(M_E)*(expr_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicExp<REAL_T>(expr_m.GetDynamicExpession());
        }

        std::string ToString() const {
            std::stringstream ss;
            ss << "Exp(" << expr_m.ToString() << ")";
            return ss.str();
        }

    private:

    };

    template<class REAL_T, class EXPR>
    inline const atl::Exp<REAL_T, EXPR> exp(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Exp<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * function used to protect overflow in exp calculations.
     *
     * Author: Dave Fournier.
     * Original implementation in ADMB.
     *
     * Source: http://admb-project.org/documentation/api/mfexp_8cpp_source.html
     *
     * @param expr
     */
    template<class REAL_T, class EXPR>
    inline const atl::Variable<REAL_T> mfexp(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        REAL_T b = REAL_T(60.0);
        if (expr.GetValue() <= b && expr.GetValue() >= REAL_T(-1) * b) {
            return atl::exp(expr);
        } else if (expr.GetValue() > b) {
            return /*std::exp(b)*/EXP_OF_B * (REAL_T(1.) + REAL_T(2.) * (expr - b)) / (REAL_T(1.) + expr - b);
        } else {
            return std::exp(REAL_T(-1) * b)*(REAL_T(1.) - expr - b) / (REAL_T(1.) + REAL_T(2.) * (REAL_T(-1) * expr - b));
        }
    }
}
namespace std {

    /**
     * Override for the exp function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Exp<REAL_T, EXPR> exp(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Exp<REAL_T, EXPR > (expr.Cast());
    }

    /**
     * Override for the exp function in namespace std.
     *
     * @param expr
     * @return
     */
    template<class REAL_T, class EXPR>
    inline const atl::Variable<REAL_T> mfexp(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::mfexp(expr);
    }
}
#endif /* EXP_HPP */

