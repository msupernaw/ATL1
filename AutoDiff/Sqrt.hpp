/* 
 * File:   Sqrt.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:31 PM
 */

#ifndef ET4AD_SQRT_HPP
#define	ET4AD_SQRT_HPP

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

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) / (2.0 * this->GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return (expr_m.EvaluateDerivative(a, b) / (2.0 * this->GetValue())) -
                    (expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)) / (4.0 * std::pow(expr_m.GetValue(), 1.5));
        }


    private:
        const EXPR& expr_m;
    } __attribute__((packed));

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

#endif	/* SQRT_HPP */

