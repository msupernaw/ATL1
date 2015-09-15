/* 
 * File:   Sinh.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:48 PM
 */

#ifndef ET4AD_SINH_HPP
#define	ET4AD_SINH_HPP

#include <cmath>
#include "Expression.hpp"
namespace atl {

    /**
     * Expression template for computing the hyperbolic sine of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Sinh;
}

namespace std {

    /**
     * Override for the sinh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Sinh<REAL_T, EXPR> sinh(const atl::ExpressionBase<REAL_T, EXPR>& expr);

}

namespace atl {

    /**
     * Expression template for computing the hyperbolic sine of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Sinh : public ExpressionBase<REAL_T, Sinh<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Sinh(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::sinh(expr_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * std::cosh(expr_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return ((std::sinh(expr_m.GetValue()) * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b))
                    + (std::cosh(expr_m.GetValue()) * expr_m.EvaluateDerivative(a, b)));
        }


    private:
        const EXPR& expr_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Sinh<REAL_T, EXPR> sinh(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Sinh<REAL_T, EXPR > (expr.Cast());
    }
}

namespace std {

    /**
     * Override for the sinh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Sinh<REAL_T, EXPR> sinh(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Sinh<REAL_T, EXPR > (expr.Cast());
    }
}

#endif	/* SINH_HPP */

