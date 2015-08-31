/* 
 * File:   ACos.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:26 PM
 */

#ifndef ET4AD_ACOS_HPP
#define	ET4AD_ACOS_HPP

#include <cmath>
#include "Expression.hpp"
namespace atl {

    /**
     * Expression template for computing the inverse cosine of an expression
     * template.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    class ACos;
}


namespace std {

    /**
     * Override for the asin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::ACos<REAL_T, EXPR> acos(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template for computing the inverse cosine of an expression
     * template.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    class ACos : public ExpressionBase<REAL_T, ACos<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        ACos(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::acos(expr_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            REAL_T fx = expr_m.GetValue();
            return -1.0 * expr_m.EvaluateDerivative(id) / std::sqrt(1.0 - fx * fx);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fx = expr_m.GetValue();
            return (((-1.0 * (fx * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)))
                    / std::pow((1.0 - fx * fx), 1.5))
                    - (expr_m.EvaluateDerivative(a, b) / std::sqrt(1.0 - fx * fx)));

        }



    private:
        const EXPR& expr_m;
    } __attribute__((packed));

    template<class REAL_T, class EXPR>
    inline const atl::ACos<REAL_T, EXPR> acos(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ACos<REAL_T, EXPR > (expr.Cast());
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
    inline const atl::ACos<REAL_T, EXPR> acos(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ACos<REAL_T, EXPR > (expr.Cast());
    }

}
#endif	/* ACOS_HPP */

