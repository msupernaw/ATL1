/* 
 * File:   Pow.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:32 PM
 */

#ifndef ET4AD_POW_HPP
#define	ET4AD_POW_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {
    template <class REAL_T, class EXPR1, class EXPR2>
    struct Pow;

    template <class REAL_T, class EXPR1>
    struct PowConstant;

    template <class REAL_T, class EXPR2>
    struct ConstantPow;
    ;
}

namespace std {

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    atl::Pow<REAL_T, EXPR1, EXPR2> pow(const atl::ExpressionBase<REAL_T, EXPR1>& expr1,
            const atl::ExpressionBase<REAL_T, EXPR2>& expr2);

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::PowConstant< REAL_T, EXPR > pow(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val);

    /**
     * Override for the pow function in namespace std.
     * 
     * @param val
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ConstantPow<REAL_T, EXPR> pow(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr);

}

namespace atl {

    /**
     *Expression template for computing the power of a expression template,
     * where both arguments are expression templates. 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    struct Pow : public ExpressionBase<REAL_T, Pow<REAL_T, EXPR1, EXPR2> > {
        typedef REAL_T BASE_TYPE;

        Pow(const ExpressionBase<REAL_T, EXPR1>& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr1_m.GetValue(), expr2_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr1_m.VariableCount(count);
            expr2_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr1_m.PushIds(ids);
            expr2_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            REAL_T g = expr2_m.GetValue();
            REAL_T f = expr1_m.GetValue();
            REAL_T fx = expr1_m.EvaluateDerivative(id);
            REAL_T gx = expr2_m.EvaluateDerivative(id);
            return std::pow(f, g)*(std::log(f) * gx + g * fx / f);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fxy = expr1_m.EvaluateDerivative(a, b);
            REAL_T g = expr2_m.GetValue();
            REAL_T f = expr1_m.GetValue();
            REAL_T fx = expr1_m.EvaluateDerivative(a);
            REAL_T gy = expr2_m.EvaluateDerivative(b);
            REAL_T fy = expr1_m.EvaluateDerivative(b);
            REAL_T gx = expr2_m.EvaluateDerivative(a);
            REAL_T gxy = expr2_m.EvaluateDerivative(a, b);
            return std::pow(f, g)*(((fx * gy) / f) + std::log(f) * gxy + (fy * gx / f) -
                    (g * fx * fy) / (f * f) + g * fxy / f) + std::pow(f, g)*(std::log(f) * gx +
                    g * fx / f)*(std::log(f) * gy + g * fy / f);

        }



        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    /**
     *Expression template for computing the power of a expression template,
     * where first argument is an expression templates, the second a constant. 
     */
    template <class REAL_T, class EXPR1>
    struct PowConstant : public ExpressionBase<REAL_T, PowConstant<REAL_T, EXPR1> > {
        typedef REAL_T BASE_TYPE;

        PowConstant(const ExpressionBase<REAL_T, EXPR1>& expr1, const REAL_T & expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr1_m.GetValue(), expr2_m);
        }

        inline void VariableCount(uint32_t& count) const {
            expr1_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr1_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
//            REAL_T g = expr2_m;
//            REAL_T f = expr1_m.GetValue();
//            REAL_T fx = expr1_m.EvaluateDerivative(id);
//            REAL_T gx = 0.0;;
//            return std::pow(f, g)*(std::log(f) * gx + g * fx / f);
            return expr2_m*expr1_m.GetValue()*expr1_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
                        REAL_T fxy = expr1_m.EvaluateDerivative(a, b);
                        REAL_T g = expr2_m;
                        REAL_T f = expr1_m.GetValue();
                        REAL_T fx = expr1_m.EvaluateDerivative(a);
                        REAL_T gy = 0.0;
                        REAL_T fy = expr1_m.EvaluateDerivative(b);
                        REAL_T gx = 0.0;
                        REAL_T gxy = 0.0;
                        return std::pow(f, g)*(((fx * gy) / f) + std::log(f) * gxy + (fy * gx / f) -
                                (g * fx * fy) / (f * f) + g * fxy / f) + std::pow(f, g)*(std::log(f) * gx +
                                g * fx / f)*(std::log(f) * gy + g * fy / f);
//            return 0.0;
        }



        const EXPR1& expr1_m;
        const REAL_T& expr2_m;
    };

    /**
     *Expression template for computing the power of a expression template,
     * where both the first argument is constant, the second is an expression template. 
     */
    template <class REAL_T, class EXPR2>
    struct ConstantPow : public ExpressionBase<REAL_T, ConstantPow<REAL_T, EXPR2> > {
        typedef REAL_T BASE_TYPE;

        ConstantPow(const REAL_T& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr1_m, expr2_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr2_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr2_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            REAL_T g = expr2_m.GetValue();
            REAL_T f = expr1_m;
            REAL_T fx = 0.0;
            REAL_T gx = expr2_m.EvaluateDerivative(id);
            return std::pow(f, g)*(std::log(f) * gx + g * fx / f);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T fxy = 0.0;
            REAL_T g = expr2_m.GetValue();
            REAL_T f = expr1_m;
            REAL_T fx = 0.0;;
            REAL_T gy = expr2_m.EvaluateDerivative(b);
            REAL_T fy = 0.0;
            REAL_T gx = expr2_m.EvaluateDerivative(a);
            REAL_T gxy = expr2_m.EvaluateDerivative(a, b);
            return std::pow(f, g)*(((fx * gy) / f) + std::log(f) * gxy + (fy * gx / f) -
                    (g * fx * fy) / (f * f) + g * fxy / f) + std::pow(f, g)*(std::log(f) * gx +
                    g * fx / f)*(std::log(f) * gy + g * fy / f);

        }



        const REAL_T& expr1_m;
        const EXPR2& expr2_m;
    };

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    atl::Pow<REAL_T, EXPR1, EXPR2> pow(const atl::ExpressionBase<REAL_T, EXPR1>& expr1,
            const atl::ExpressionBase<REAL_T, EXPR2>& expr2) {
        return atl::Pow<REAL_T, EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    atl::Pow<REAL_T, EXPR1, EXPR2> operator^(const atl::ExpressionBase<REAL_T, EXPR1>& expr1,
            const atl::ExpressionBase<REAL_T, EXPR2>& expr2) {
        return atl::Pow<REAL_T, EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::PowConstant< REAL_T, EXPR > pow(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::PowConstant< REAL_T, EXPR > (expr.Cast(), val);
    }

    template <class REAL_T, class EXPR>
    inline
    atl::PowConstant< REAL_T, EXPR > operator^(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::PowConstant< REAL_T, EXPR > (expr.Cast(), val);
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param val
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ConstantPow<REAL_T, EXPR> pow(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ConstantPow<REAL_T, EXPR > (val, expr.Cast());
    }

    template <class REAL_T, class EXPR>
    inline
    atl::ConstantPow<REAL_T, EXPR> operator^(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::ConstantPow<REAL_T, EXPR > (val, expr.Cast());
    }

}
namespace std {

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    atl::Pow<REAL_T, EXPR1, EXPR2> pow(const atl::ExpressionBase<REAL_T, EXPR1>& expr1,
            const atl::ExpressionBase<REAL_T, EXPR2>& expr2) {
        return atl::Pow<REAL_T, EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::PowConstant< REAL_T, EXPR > pow(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::PowConstant< REAL_T, EXPR > (expr.Cast(), val);
    }

    /**
     * Override for the pow function in namespace std.
     * 
     * @param val
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ConstantPow<REAL_T, EXPR> pow(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::ConstantPow<REAL_T, EXPR > (val, expr.Cast());
    }

}
#endif	/* POW_HPP */

