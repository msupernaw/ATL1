/* 
 * File:   MFExp.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:44 PM
 */

#ifndef ET4AD_MFEXP_HPP
#define	ET4AD_MFEXP_HPP

#include <cmath>
#include "Expression.hpp"

#define EXP_OF_B REAL_T(114200738981568423454048256.0)


namespace atl {

    /**
     * Expression template used to protect overflow in exp calculations. 
     * 
     * Author: Dave Fournier.
     * Original implementation in ADMB.
     * 
     * Source: http://admb-project.org/documentation/api/mfexp_8cpp_source.html
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class MFExp;
}

namespace std {

    template<class REAL_T, class EXPR>
    inline const atl::MFExp<REAL_T, EXPR> mfexp(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template used to protect overflow in exp calculations. 
     * 
     * Author: Dave Fournier.
     * Original implementation in ADMB.
     * 
     * Source: http://admb-project.org/documentation/api/mfexp_8cpp_source.html
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class MFExp : public ExpressionBase<REAL_T, MFExp<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        MFExp(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(Compute(expr.GetValue())) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline const REAL_T Compute(const REAL_T & value) {
            //            REAL_T x = value;
            REAL_T b = REAL_T(60);
            if (value <= b && value >= REAL_T(-1) * b) {
                return std::exp(value);
            } else if (value > b) {
                return /*std::exp(b)*/EXP_OF_B * (REAL_T(1.) + REAL_T(2.) * (value - b)) / (REAL_T(1.) + value - b);
            } else {
                return std::exp(REAL_T(-1) * b)*(REAL_T(1.) - value - b) / (REAL_T(1.) + REAL_T(2.) * (REAL_T(-1) * value - b));
            }
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
        const EXPR& expr_m;
        REAL_T value_m;


    };

}

namespace std {

    template<class REAL_T, class EXPR>
    inline const atl::MFExp<REAL_T, EXPR> mfexp(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::MFExp<REAL_T, EXPR > (expr.Cast());
    }
}

#endif	/* MFEXP_HPP */

