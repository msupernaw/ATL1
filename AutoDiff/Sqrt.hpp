/* 
 * File:   Sqrt.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:31 PM
 */

#ifndef ET4AD_SQRT_HPP
#define ET4AD_SQRT_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {

    /**
     * Expression template for computing the square root of an expression
     * template.
     * @param a
     */
    template <class REAL_T, class EXPR>
    class Sqrt;

}

namespace std {

    /**
     * Override for the sqrt function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Sqrt<REAL_T, EXPR> sqrt(const atl::ExpressionBase<REAL_T, EXPR>& expr);

}

namespace atl {

    /**
     * Expression template for computing the square root of an expression
     * template.
     * @param a
     */
    template <class REAL_T, class EXPR>
    class Sqrt : public ExpressionBase<REAL_T, Sqrt<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Sqrt(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::sqrt(expr_m.GetValue());
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
            expr_m.MakeNLInteractions(b);
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) / (2.0 * this->GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return (expr_m.EvaluateDerivative(a, b) / (2.0 * this->GetValue())) -
                    (expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)) / (4.0 * std::pow(expr_m.GetValue(), 1.5));
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return (3.0 * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))
                    *(expr_m.EvaluateDerivative(z))) / (8.0 * std::pow(expr_m.GetValue(), 5.0 / 2.0))
                    -((expr_m.EvaluateDerivative(x, y))*(expr_m.EvaluateDerivative(z)))
                    / (4.0 * std::pow(expr_m.GetValue(), 3.0 / 2.0))-((expr_m.EvaluateDerivative(x))
                    *(expr_m.EvaluateDerivative(y, z))) / (4 * std::pow(expr_m.GetValue(), 3.0 / 2.0))
                    -((expr_m.EvaluateDerivative(x, z))*(expr_m.EvaluateDerivative(y)))
                    / (4.0 * std::pow(expr_m.GetValue(), 3.0 / 2.0)) + expr_m.EvaluateDerivative(x, y, z)
                    / (2.0 * std::sqrt(expr_m.GetValue()));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicSqrt<REAL_T>(expr_m.GetDynamicExpession());
        }


    private:
        const EXPR& expr_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Sqrt<REAL_T, EXPR> sqrt(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Sqrt<REAL_T, EXPR > (expr.Cast());
    }

}

namespace std {

    /**
     * Override for the sqrt function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Sqrt<REAL_T, EXPR> sqrt(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Sqrt<REAL_T, EXPR > (expr.Cast());
    }
}

#endif /* SQRT_HPP */

