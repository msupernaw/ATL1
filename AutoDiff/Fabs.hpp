/* 
 * File:   Fabs.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:53 PM
 */

#ifndef ET4AD_FABS_HPP
#define	ET4AD_FABS_HPP

#include <cmath>
#include "Expression.hpp"
namespace atl {

    /**
     * Expression template for handling the absolute value of a expression 
     * template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Fabs;

}

namespace std {

    /**
     * Override for the fabs function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Fabs<REAL_T, EXPR> fabs(const atl::ExpressionBase<REAL_T, EXPR>& expr);

}

namespace atl {

    /**
     * Expression template for handling the absolute value of a expression 
     * template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Fabs : public ExpressionBase<REAL_T, Fabs<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Fabs(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::fabs(expr_m.GetValue());
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (expr_m.EvaluateDerivative(id) * expr_m.GetValue()) / this->GetValue();
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return (expr_m.EvaluateDerivative(a, b) * expr_m.GetValue()) / this->GetValue();
        }

    private:
        const EXPR& expr_m;
    } __attribute__((packed));

    template<class REAL_T, class EXPR>
    inline const atl::Fabs<REAL_T, EXPR> fabs(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Fabs<REAL_T, EXPR > (expr.Cast());
    }

}

namespace std {

    /**
     * Override for the fabs function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Fabs<REAL_T, EXPR> fabs(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Fabs<REAL_T, EXPR > (expr.Cast());
    }

}

#endif	/* FABS_HPP */

