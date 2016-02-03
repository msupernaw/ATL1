/* 
 * File:   Tanh.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:51 PM
 */

#ifndef ET4AD_TANH_HPP
#define ET4AD_TANH_HPP

#include <cmath>
#include "Expression.hpp"


namespace atl {

    /**
     * Expression template for computing the hyperbolic tangent of an expresison
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Tanh;

}

namespace std {

    /**
     * Override for the tanh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Tanh<REAL_T, EXPR> tanh(const atl::ExpressionBase<REAL_T, EXPR>& expr);
}

namespace atl {

    /**
     * Expression template for computing the hyperbolic tangent of an expresison
     * template.
     * 
     * @param expr
     */
    template <class REAL_T, class EXPR>
    class Tanh : public ExpressionBase<REAL_T, Tanh<REAL_T, EXPR> > {
        const REAL_T sech;
        REAL_T value_m;
    public:
        typedef REAL_T BASE_TYPE;

        Tanh(const ExpressionBase<REAL_T, EXPR>& expr)
        : expr_m(expr.Cast()), value_m(std::tanh(expr.GetValue())), sech(1.0 / std::cosh(expr.GetValue())) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
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

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr_m.EvaluateDerivative(id) * sech*sech;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T sech2 = sech*sech;
            return sech2 * expr_m.EvaluateDerivative(a, b) - 2.0 * sech2 * this->GetValue() * expr_m.EvaluateDerivative(a) * expr_m.EvaluateDerivative(b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return 4.0 * std::pow(sech,2.0) * std::pow(std::tanh(expr_m.GetValue()),2.0)
                    * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))
                    *(expr_m.EvaluateDerivative(z)) - 2.0 * std::pow(sech,4.0) 
                    * (expr_m.EvaluateDerivative(x))*(expr_m.EvaluateDerivative(y))
                    *(expr_m.EvaluateDerivative(z)) - 2.0 * std::pow(sech,2.0)
                    * std::tanh(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, y))*
                    (expr_m.EvaluateDerivative(z)) - 2.0 * std::pow(sech,2.0) 
                    * std::tanh(expr_m.GetValue())*(expr_m.EvaluateDerivative(x))
                    *(expr_m.EvaluateDerivative(y, z)) - 2.0 * std::pow(sech,2.0) 
                    * std::tanh(expr_m.GetValue())*(expr_m.EvaluateDerivative(x, z))
                    *(expr_m.EvaluateDerivative(y)) + std::pow(sech,2.0)
                    * (expr_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicTanh<REAL_T>(expr_m.GetDynamicExpession());
        }

    private:
        const EXPR& expr_m;
    };

    /**
     * Override for the tanh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Tanh<REAL_T, EXPR> tanh(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Tanh<REAL_T, EXPR > (expr.Cast());
    }
}

namespace std {

    /**
     * Override for the tanh function in namespace std.
     * 
     * @param expr
     * @return 
     */
    template<class REAL_T, class EXPR>
    inline const atl::Tanh<REAL_T, EXPR> tanh(const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::Tanh<REAL_T, EXPR > (expr.Cast());
    }
}


#endif /* TANH_HPP */

