/* 
 * File:   Log.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:34 PM
 */

#ifndef ET4AD_LOG_HPP
#define	ET4AD_LOG_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {

    /**
     * Expression template to compute the log of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Log;

}

namespace std {

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Log<REAL_T, EXPR> log(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template to compute the log of an expression template.
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Log : public ExpressionBase<REAL_T, Log<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Log(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(expr_m.GetValue()), value2_m(REAL_T(1.0) / value_m) {
        }

        inline const REAL_T GetValue() const {
            return std::log(value_m);
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }


        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent = true)const {
            expr_m.PushIds(ids,include_dependent);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) /expr_m.GetValue(); // expr_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            //            REAL_T fx =expr_m.GetValue();
            return (expr_m.EvaluateDerivative(a, b) /expr_m.GetValue()) - (expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)) / (expr_m.GetValue()*expr_m.GetValue());
        }

    private:
        const EXPR& expr_m;
        const REAL_T value_m;
        const REAL_T value2_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Log<REAL_T, EXPR> log(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Log<REAL_T, EXPR > (expr.Cast());
    }
}

namespace std {

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Log<REAL_T, EXPR> log(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Log<REAL_T, EXPR > (expr.Cast());
    }

}

#endif	/* LOG_HPP */
