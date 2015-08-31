/* 
 * File:   Tan.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:17 PM
 */

#ifndef ET4AD_TAN_HPP
#define	ET4AD_TAN_HPP

#include <cmath>
#include "Expression.hpp"


namespace atl {

    /**
     * Expression template for computing the tangent of a expression template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Tan;
}


namespace std {

    /**
     * Override for the tan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Tan<REAL_T, EXPR> tan(const atl::ExpressionBase<REAL_T, EXPR>& expr);

}



namespace atl {

    /**
     * Expression template for computing the tangent of a expression template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Tan : public ExpressionBase<REAL_T, Tan<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Tan(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::tan(expr_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * (1.0 / std::cos(expr_m.GetValue()))*(1.0 / std::cos(expr_m.GetValue()));
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T sec2 = (1.0 / std::cos(expr_m.GetValue())) * (1.0 / std::cos(expr_m.GetValue()));
            return 2.0 * sec2 * this->GetValue() * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b) +
                    sec2 * expr_m.EvaluateDerivative(a, b);
        }


    private:
        const EXPR& expr_m;
    } __attribute__((packed));

    template<class REAL_T, class EXPR>
    inline const atl::Tan<REAL_T, EXPR> tan(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Tan<REAL_T, EXPR > (expr.Cast());
    }

}

namespace std {

    /**
     * Override for the tan function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Tan<REAL_T, EXPR> tan(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Tan<REAL_T, EXPR > (expr.Cast());
    }

}
#endif	/* TAN_HPP */

