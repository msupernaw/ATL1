/* 
 * File:   Sin.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:04 PM
 */

#ifndef ET4AD_SIN_HPP
#define	ET4AD_SIN_HPP

#include <cmath>

#include "Expression.hpp"


namespace atl {

    /**
     * Expression template for taking the sine of an expression.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    class Sin;
}


namespace std {

    /**
     * Override for the sin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Sin<REAL_T, EXPR> sin(const atl::ExpressionBase<REAL_T, EXPR>& expr);

}

namespace atl {

    /**
     * Expression template for taking the sine of an expression.
     * 
     * @param a
     */
    template <class REAL_T, class EXPR>
    class Sin : public ExpressionBase<REAL_T, Sin<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Sin(const ExpressionBase<REAL_T, EXPR>& a)
        : expr_m(a.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::sin(expr_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * std::cos(expr_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return (std::cos(expr_m.GetValue()) * expr_m.EvaluateDerivative(a, b))-
                    std::sin(expr_m.GetValue()) * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b);
        }



    private:
        const EXPR& expr_m;
    } __attribute__((packed));

    template<class REAL_T, class EXPR>
    inline const atl::Sin<REAL_T, EXPR> sin(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Sin<REAL_T, EXPR > (expr.Cast());
    }

}
namespace std {

    /**
     * Override for the sin function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Sin<REAL_T, EXPR> sin(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Sin<REAL_T, EXPR > (expr.Cast());
    }
}
#endif	/* SIN_HPP */

