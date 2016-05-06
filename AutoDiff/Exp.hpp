/* 
 * File:   Exp.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:37 PM
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

        Exp(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(Compute(expr.GetValue())) {
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
            //            return std::exp(expr_m.GetValue())*(expr_m.EvaluateDerivative(a))*(expr_m.EvaluateDerivative(b)) + std::exp(expr_m.GetValue())*(expr_m.EvaluateDerivative(a, b));

            return std::exp(expr_m.GetValue()) * std::pow(std::log(M_E), 2.0)*(expr_m.EvaluateDerivative(a))*(expr_m.EvaluateDerivative(b)) + std::exp(expr_m.GetValue()) * std::log(M_E)*(expr_m.EvaluateDerivative(a, b));

        }
       
        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            //            return myexp(expr_m.GetValue())*(expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))*(expr_m.EvaluateDerivative(z)) + myexp(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, y))*(expr_m.EvaluateDerivative(z)) +
            //                    myexp(expr_m.GetValue())*(expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y, z)) + myexp(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, z))*(expr_m.EvaluateDerivative(y)) + myexp(expr_m.GetValue())*
            //                    (expr_m.EvaluateDerivative(x, y, z));

            return std::exp( expr_m.GetValue()) * std::pow(std::log(M_E), 3.0) * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))*(expr_m.EvaluateDerivative(z)) + std::exp( expr_m.GetValue()) * std::pow(std::log(M_E), 2.0) * (expr_m.EvaluateDerivative(x, y))*
                    (expr_m.EvaluateDerivative(z)) + std::exp( expr_m.GetValue()) * std::pow(std::log(M_E), 2.0) * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y, z)) + std::exp( expr_m.GetValue()) * std::pow(std::log(M_E), 2.0)* (expr_m.EvaluateDerivative(x, z))*(expr_m.EvaluateDerivative(y))
                    + std::exp( expr_m.GetValue()) * std::log(M_E)*(expr_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicExp<REAL_T>(expr_m.GetDynamicExpession());
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
}
#endif /* EXP_HPP */

