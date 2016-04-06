/* 
 * File:   Log10.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:35 PM
 */

#ifndef ET4AD_LOG10_HPP
#define ET4AD_LOG10_HPP

#include <cmath>
#include "Expression.hpp"

#define AD_LOG10 2.30258509299404590109361379290930926799774169921875

namespace atl {

    /**
     * Expression template to compute the log base 10 of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Log10;
}

namespace std {

    /**
     * Override for the log function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Log10<REAL_T, EXPR> log10(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}


namespace atl {

    /**
     * Expression template to compute the log base 10 of an expression 
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Log10 : public ExpressionBase<REAL_T, Log10<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Log10(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::log10(expr_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            expr_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr_m.PushIds(ids);
        }

        bool IsNonlinear()const {
            return true;
        }

        inline void MakeNLInteractions(bool b = false)const {
            //            expr_m.MakeNLInteractions();
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            //            expr_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (expr_m.EvaluateDerivative(id) / (AD_LOG10 * expr_m.GetValue()));
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fx = expr_m.GetValue();
            return (expr_m.EvaluateDerivative(a, b) / (AD_LOG10 * fx)) -
                    ((expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b)) / (AD_LOG10 * (fx * fx)));
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return (2.0 * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))
                    *(expr_m.EvaluateDerivative(z))) / (AD_LOG10 * std::pow(expr_m.GetValue(), 3.0))
                    -((expr_m.EvaluateDerivative(x, y))*(expr_m.EvaluateDerivative(z)))
                    / (AD_LOG10 * std::pow(expr_m.GetValue(), 2.0))-((expr_m.EvaluateDerivative(x))
                    *(expr_m.EvaluateDerivative(y, z))) / (AD_LOG10 * std::pow(expr_m.GetValue(), 2.0))
                    -((expr_m.EvaluateDerivative(x, z))*(expr_m.EvaluateDerivative(y)))
                    / (AD_LOG10 * std::pow(expr_m.GetValue(), 2.0)) +
                    expr_m.EvaluateDerivative(x, y, z) / (AD_LOG10 * expr_m.GetValue());
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicLog10<REAL_T>(expr_m.GetDynamicExpession());
        }


    private:
        const EXPR& expr_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Log10<REAL_T, EXPR> log10(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Log10<REAL_T, EXPR > (expr.Cast());
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
    inline const atl::Log10<REAL_T, EXPR> log10(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Log10<REAL_T, EXPR > (expr.Cast());
    }
}

#endif /* LOG10_HPP */

