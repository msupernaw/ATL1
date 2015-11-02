/* 
 * File:   Cosh.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:50 PM
 */

#ifndef ET4AD_COSH_HPP
#define	ET4AD_COSH_HPP


#include <cmath>
#include "Expression.hpp"


namespace atl {

    /**
     * Expression template for computing the hyperbolic cosine of an 
     * expression template.
     */
    template <class REAL_T, class EXPR>
    class Cosh;

}

namespace std {
    /**
     * Override for the cosh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Cosh<REAL_T, EXPR> cosh(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template for computing the hyperbolic cosine of an 
     * expression template.
     */
    template <class REAL_T, class EXPR>
    class Cosh : public ExpressionBase<REAL_T, Cosh<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Cosh(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(expr.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return std::cosh(value_m);
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent = true)const {
            expr_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * std::sinh(expr_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return ((std::cosh(expr_m.GetValue()) * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b))
                    + (std::sinh(expr_m.GetValue()) * expr_m.EvaluateDerivative(a, b)));
        }


    private:
        const EXPR& expr_m;
        REAL_T value_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Cosh<REAL_T, EXPR> cosh(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Cosh<REAL_T, EXPR > (expr.Cast());
    }

}

namespace std {

    /**
     * Override for the cosh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Cosh<REAL_T, EXPR> cosh(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Cosh<REAL_T, EXPR > (expr.Cast());
    }
}



#endif	/* COSH_HPP */

