/* 
 * File:   ASin.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:23 PM
 */

#ifndef ET4AD_ASIN_HPP
#define	ET4AD_ASIN_HPP

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

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
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




    private:
        const EXPR& expr_m;
    } __attribute__((packed));

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
#endif	/* ASIN_HPP */

