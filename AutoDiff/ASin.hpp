/* 
 * File:   ASin.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:23 PM
 */

#ifndef ET4AD_ASIN_HPP
#define ET4AD_ASIN_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {

    /**
     * Expression template for computing the inverse Sine of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class ASin;
}

namespace std {

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::ASin<REAL_T, EXPR> asin(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template for computing the inverse Sine of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class ASin : public ExpressionBase<REAL_T, ASin<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        ASin(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::asin(expr_m.GetValue());
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

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            expr_m.PushAdjoints(adjoints, coefficient *
                    static_cast<REAL_T> (1.0) /
                    std::pow((static_cast<REAL_T> (1.0) -
                    std::pow(expr_m.GetValue(),
                    static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
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
            return expr_m.EvaluateDerivative(id) / std::sqrt(1.0 - fx * fx);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fx = expr_m.GetValue();
            return ((((fx * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)))
                    / std::pow((1.0 - fx * fx), 1.5))
                    +(expr_m.EvaluateDerivative(a, b) / std::sqrt(1.0 - fx * fx)));

        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return ((expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))*
                    (expr_m.EvaluateDerivative(z))) /
                    std::pow((1 - std::pow(expr_m.GetValue(), 2.0)), (3.0 / 2.0))
                    +(3 * std::pow(expr_m.GetValue(), 2.0)*(expr_m.EvaluateDerivative(x))
                    *(expr_m.EvaluateDerivative(y))*(expr_m.EvaluateDerivative(z))) /
                    std::pow((1 - std::pow(expr_m.GetValue(), 2.0)), (5.0 / 2.0))+
                    (expr_m.GetValue()*(expr_m.EvaluateDerivative(x, y))*
                    (expr_m.EvaluateDerivative(z))) /
                    std::pow((1 - std::pow(expr_m.GetValue(), 2.0)), (3.0 / 2.0))
                    +(expr_m.GetValue()*(expr_m.EvaluateDerivative(x))*
                    (expr_m.EvaluateDerivative(y, z))) /
                    std::pow((1 - std::pow(expr_m.GetValue(), 2.0)), (3.0 / 2.0))+
                    (expr_m.GetValue()*(expr_m.EvaluateDerivative(x, z))*
                    (expr_m.EvaluateDerivative(y))) /
                    std::pow((1 - std::pow(expr_m.GetValue(), 2.0)), (3.0 / 2.0)) +
                    expr_m.EvaluateDerivative(x, y, z) /
                    std::sqrt(1 - std::pow(expr_m.GetValue(), 2.0));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicASin<REAL_T>(expr_m.GetDynamicExpession());
        }



    private:
        const EXPR& expr_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::ASin<REAL_T, EXPR> asin(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ASin<REAL_T, EXPR > (expr.Cast());
    }

}

namespace std {

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::ASin<REAL_T, EXPR> asin(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ASin<REAL_T, EXPR > (expr.Cast());
    }

}
#endif /* ASIN_HPP */

