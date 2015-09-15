/* 
 * File:   Exp.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:37 PM
 */

#ifndef ET4AD_EXP_HPP
#define	ET4AD_EXP_HPP

#include <cmath>
#include "Expression.hpp"

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

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * value_m;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fx = value_m;
            return ((fx * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b))
                    + (fx * expr_m.EvaluateDerivative(a, b)));
        }

    private:

    };

    template<class REAL_T, class EXPR>
    inline const atl::Exp<REAL_T, EXPR> exp(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Exp<REAL_T, EXPR > (expr.Cast());
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
#endif	/* EXP_HPP */

