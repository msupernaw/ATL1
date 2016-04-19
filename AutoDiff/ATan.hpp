/* 
 * File:   ATan.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:27 PM
 */

#ifndef ET4AD_ATAN_HPP
#define ET4AD_ATAN_HPP

#include <cmath>
#include "Expression.hpp"
namespace atl {

    /**
     * Expression template for computing the inverse tangent of an expression 
     * template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class ATan;
}

namespace std {

    /**
     * Override for the atan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::ATan<REAL_T, EXPR> atan(const atl::ExpressionBase<REAL_T, EXPR>& a);
}
namespace atl {

    /**
     * Expression template for computing the inverse tangent of an expression 
     * template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class ATan : public ExpressionBase<REAL_T, ATan<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        ATan(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(expr.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan(value_m);
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
            return true;
        }

        bool IsNonlinear()const {
            return true;
        }

        inline void MakeNLInteractions(bool b = false)const {
            //            Cast().MakeNLInteractions();
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            //            expr_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            REAL_T fx = expr_m.GetValue();
            return expr_m.EvaluateDerivative(id) / (fx * fx + 1.0);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fx = expr_m.GetValue();
            return (expr_m.EvaluateDerivative(a, b) / (fx * fx + 1.0)) -
                    (2.0 * fx * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)) / ((fx * fx + 1.0)*(fx * fx + 1.0));

        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return -(2.0 * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))
                    *(expr_m.EvaluateDerivative(z))) /
                    std::pow((std::pow(expr_m.GetValue(), 2.0) + 1), 2.0)+
                    (8.0 * std::pow(expr_m.GetValue(), 2.0)*
                    (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))*
                    (expr_m.EvaluateDerivative(z))) /
                    std::pow((std::pow(expr_m.GetValue(), 2.0) + 1), 3.0)-
                    (2.0 * expr_m.GetValue()*(expr_m.EvaluateDerivative(x, y))*
                    (expr_m.EvaluateDerivative(z))) /
                    std::pow((std::pow(expr_m.GetValue(), 2.0) + 1), 2.0)-
                    (2.0 * expr_m.GetValue()*(expr_m.EvaluateDerivative(x))*
                    (expr_m.EvaluateDerivative(y, z))) /
                    std::pow((std::pow(expr_m.GetValue(), 2.0) + 1), 2.0)-
                    (2.0 * expr_m.GetValue()*(expr_m.EvaluateDerivative(x, z))*
                    (expr_m.EvaluateDerivative(y))) /
                    std::pow((std::pow(expr_m.GetValue(), 2.0) + 1), 2.0) +
                    expr_m.EvaluateDerivative(x, y, z) /
                    (std::pow(expr_m.GetValue(), 2.0) + 1);
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicATan<REAL_T>(expr_m.GetDynamicExpession());
        }

    private:
        const EXPR& expr_m;
        REAL_T value_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::ATan<REAL_T, EXPR> atan(const atl::ExpressionBase<REAL_T, EXPR>& a) {
        return atl::ATan<REAL_T, EXPR > (a.Cast());
    }

}

namespace std {

    /**
     * Override for the atan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::ATan<REAL_T, EXPR> atan(const atl::ExpressionBase<REAL_T, EXPR>& a) {
        return atl::ATan<REAL_T, EXPR > (a.Cast());
    }
}

#endif /* ATAN_HPP */

