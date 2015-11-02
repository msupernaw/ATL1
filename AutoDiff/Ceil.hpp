/* 
 * File:   Ceil.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:56 PM
 */

#ifndef ET4AD_CEIL_HPP
#define	ET4AD_CEIL_HPP

#include <cmath>
#include "Expression.hpp"


namespace atl {

    /**
     * Expression template for handling the ceiling of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Ceil;
}

namespace std {

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Ceil<REAL_T, EXPR> ceil(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template for handling the ceiling of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Ceil : public ExpressionBase<REAL_T, Ceil<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Ceil(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::ceil(expr_m.GetValue());
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
            return 0.0; //expr_m.EvaluateDerivative(id) * GetValue();
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return 0.0; //expr_m.EvaluateDerivative(a, b) * GetValue();

        }


    private:
        const EXPR& expr_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Ceil<REAL_T, EXPR> ceil(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Ceil<REAL_T, EXPR > (expr.Cast());
    }


}

namespace std {

    /**
     * Override for the ceil function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Ceil<REAL_T, EXPR> ceil(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Ceil<REAL_T, EXPR > (expr.Cast());
    }
}


#endif	/* CEIL_HPP */

