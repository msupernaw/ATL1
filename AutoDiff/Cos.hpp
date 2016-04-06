/* 
 * File:   Cos.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:15 PM
 */

#ifndef ET4AD_COS_HPP
#define ET4AD_COS_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {

    /**
     * Expression template for taking the cosine of an expression.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Cos;

}


namespace std {

    /**
     * Override for the cos function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Cos<REAL_T, EXPR> cos(const atl::ExpressionBase<REAL_T, EXPR>& expr);

}


namespace atl {

    /**
     * Expression template for taking the cosine of an expression.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Cos : public ExpressionBase<REAL_T, Cos<REAL_T, EXPR> > {
    public:
        typedef REAL_T BASE_TYPE;

        Cos(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(expr.GetValue()) {

        }

        inline const REAL_T GetValue() const {
            return std::cos(value_m);
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
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
//            expr_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return -1.0 * expr_m.EvaluateDerivative(id) * std::sin(expr_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return ((-1.0 * std::cos(expr_m.GetValue()) * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b))
                    - (1.0 * std::sin(expr_m.GetValue()) * expr_m.EvaluateDerivative(a, b)));
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return std::sin(expr_m.GetValue())*(expr_m.EvaluateDerivative(x))*
                    (expr_m.EvaluateDerivative(y))*(expr_m.EvaluateDerivative(z))
                    - std::cos(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, y))
                    *(expr_m.EvaluateDerivative(z)) - std::cos(expr_m.GetValue())
                    *(expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y, z))
                    - std::cos(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, z))
                    *(expr_m.EvaluateDerivative(y)) -
                    std::sin(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicCos<REAL_T>(expr_m.GetDynamicExpession());
        }



    private:
        const EXPR& expr_m;
        REAL_T value_m;
    };

    template<class REAL_T, class EXPR>
    inline const atl::Cos<REAL_T, EXPR> cos(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Cos<REAL_T, EXPR > (expr.Cast());
    }



}

namespace std {

    /**
     * Override for the cos function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Cos<REAL_T, EXPR> cos(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::Cos<REAL_T, EXPR > (expr.Cast());
    }
}

#endif /* COS_HPP */

