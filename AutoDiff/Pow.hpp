/* 
 * File:   Pow.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:32 PM
 */

#ifndef ET4AD_POW_HPP
#define ET4AD_POW_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {
    template <class REAL_T, class EXPR1, class EXPR2>
    struct Pow;

    template <class REAL_T, class EXPR1>
    struct PowScalar;

    template <class REAL_T, class EXPR2>
    struct ScalarPow;
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
    atl::PowScalar< REAL_T, EXPR > pow(const atl::ExpressionBase<REAL_T, EXPR>& expr,
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
    atl::ScalarPow<REAL_T, EXPR> pow(const REAL_T& val,
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

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            expr1_m.PushIds(ids, include_dependent);
            expr2_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr1_m.PushIds(ids);
            expr2_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr1_m.PushIds(ids);
            expr2_m.PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            expr1_m.PushAdjoints(adjoints, coefficient * expr2_m.GetValue() * pow(expr1_m.GetValue(), expr2_m.GetValue() - static_cast<REAL_T> (1.0)));
            expr2_m.PushAdjoints(adjoints, coefficient * this->GetValue() * std::log(expr1_m.GetValue()));
        }

        bool IsNonlinear()const {
            return true;
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr1_m.PushNLInteractions(ids);
            expr2_m.PushNLInteractions(ids);
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

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return std::pow(expr1_m.GetValue(), expr2_m.GetValue())*(-((expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(y))*(expr2_m.EvaluateDerivative(z))) / std::pow(expr1_m.GetValue(), 2.0)+((expr1_m.EvaluateDerivative(x, y))*(expr2_m.EvaluateDerivative(z))) / expr1_m.GetValue()+((expr1_m.EvaluateDerivative(x))*(expr2_m.EvaluateDerivative(y, z))) / expr1_m.GetValue()-((expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(z))*(expr2_m.EvaluateDerivative(y))) / std::pow(expr1_m.GetValue(), 2.0)+((expr1_m.EvaluateDerivative(x, z))*(expr2_m.EvaluateDerivative(y))) / expr1_m.GetValue()+
                    ((expr1_m.EvaluateDerivative(y))*(expr2_m.EvaluateDerivative(x, z))) / expr1_m.GetValue() + std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(x, y, z))+((expr1_m.EvaluateDerivative(z))*(expr2_m.EvaluateDerivative(x, y))) / expr1_m.GetValue()-((expr1_m.EvaluateDerivative(y))*(expr1_m.EvaluateDerivative(z))*(expr2_m.EvaluateDerivative(x))) / std::pow(expr1_m.GetValue(), 2.0)+((expr1_m.EvaluateDerivative(y, z))*(expr2_m.EvaluateDerivative(x))) / expr1_m.GetValue()+
                    (2.0 * expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(y))*(expr1_m.EvaluateDerivative(z))) / std::pow(expr1_m.GetValue(), 3.0) - (expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x, y))*(expr1_m.EvaluateDerivative(z))) / std::pow(expr1_m.GetValue(), 2.0)-(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(y, z))) / std::pow(expr1_m.GetValue(), 2.0)-(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x, z))*(expr1_m.EvaluateDerivative(y))) / std::pow(expr1_m.GetValue(), 2.0)+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x, y, z))) / expr1_m.GetValue()) +
                    std::pow(expr1_m.GetValue(), expr2_m.GetValue())*(std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(x))+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x))) / expr1_m.GetValue())*
                    (((expr1_m.EvaluateDerivative(y))*(expr2_m.EvaluateDerivative(z))) / expr1_m.GetValue() + std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(y, z))+((expr1_m.EvaluateDerivative(z))*(expr2_m.EvaluateDerivative(y))) / expr1_m.GetValue()-(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(y))*(expr1_m.EvaluateDerivative(z))) / std::pow(expr1_m.GetValue(), 2.0)+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(y, z))) / expr1_m.GetValue()) + std::pow(expr1_m.GetValue(), expr2_m.GetValue())*
                    (std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(y))+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(y))) / expr1_m.GetValue())*
                    (((expr1_m.EvaluateDerivative(x))*(expr2_m.EvaluateDerivative(z))) / expr1_m.GetValue() + std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(x, z))+((expr1_m.EvaluateDerivative(z))*(expr2_m.EvaluateDerivative(x))) / expr1_m.GetValue()-(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(z))) / std::pow(expr1_m.GetValue(), 2.0)+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x, z))) / expr1_m.GetValue()) + std::pow(expr1_m.GetValue(), expr2_m.GetValue())*
                    (((expr1_m.EvaluateDerivative(x))*(expr2_m.EvaluateDerivative(y))) / expr1_m.GetValue() + std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(x, y))+((expr1_m.EvaluateDerivative(y))*(expr2_m.EvaluateDerivative(x))) / expr1_m.GetValue()-(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(y))) / std::pow(expr1_m.GetValue(), 2.0)+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x, y))) / expr1_m.GetValue())*
                    (std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(z))+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(z))) / expr1_m.GetValue()) + std::pow(expr1_m.GetValue(), expr2_m.GetValue())*(std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(x))+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(x))) / expr1_m.GetValue())*(std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(y))+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(y))) / expr1_m.GetValue())*
                    (std::log(expr1_m.GetValue())*(expr2_m.EvaluateDerivative(z))+(expr2_m.GetValue()*(expr1_m.EvaluateDerivative(z))) / expr1_m.GetValue());
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicPow<REAL_T>(expr1_m.GetDynamicExpession(), expr2_m.GetDynamicExpession());
        }



        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
    };

    /**
     *Expression template for computing the power of a expression template,
     * where first argument is an expression templates, the second a constant. 
     */
    template <class REAL_T, class EXPR1>
    struct PowScalar : public ExpressionBase<REAL_T, PowScalar<REAL_T, EXPR1> > {
        typedef REAL_T BASE_TYPE;

        PowScalar(const ExpressionBase<REAL_T, EXPR1>& expr1, const REAL_T & expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr1_m.GetValue(), expr2_m);
        }

        inline void VariableCount(uint32_t& count) const {
            expr1_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            expr1_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr1_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr1_m.PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            expr1_m.PushAdjoints(adjoints, coefficient * expr2_m * pow(expr1_m.GetValue(), expr2_m - static_cast<REAL_T> (1.0)));
        }

        bool IsNonlinear()const {
            return true;
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr1_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 1.0) * expr1_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return ((expr2_m - 1) * expr2_m)*std::pow(expr1_m.GetValue(), expr2_m - 2.0) * expr1_m.EvaluateDerivative(a) * expr1_m.EvaluateDerivative(b) +
                    expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 1.0) * expr1_m.EvaluateDerivative(a, b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return (expr2_m - 2.0)*(expr2_m - 1.0) * expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 3.0)
                    *(expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(y))
                    *(expr1_m.EvaluateDerivative(z))+(expr2_m - 1.0)
                    * expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 2.0)*(expr1_m.EvaluateDerivative(x, y))
                    *(expr1_m.EvaluateDerivative(z))+(expr2_m - 1.0) * expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 2.0)*
                    (expr1_m.EvaluateDerivative(x))*(expr1_m.EvaluateDerivative(y, z))+(expr2_m - 1)
                    * expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 2.0) * (expr1_m.EvaluateDerivative(x, z))
                    *(expr1_m.EvaluateDerivative(y)) + expr2_m * std::pow(expr1_m.GetValue(), expr2_m - 1.0)
                    * (expr1_m.EvaluateDerivative(x, y, z));
        }



        const EXPR1& expr1_m;
        const REAL_T& expr2_m;
    };

    /**
     *Expression template for computing the power of a expression template,
     * where both the first argument is constant, the second is an expression template. 
     */
    template <class REAL_T, class EXPR2>
    struct ScalarPow : public ExpressionBase<REAL_T, ScalarPow<REAL_T, EXPR2> > {
        typedef REAL_T BASE_TYPE;

        ScalarPow(const REAL_T& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::pow(expr1_m, expr2_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr2_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            expr2_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            expr2_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            expr2_m.PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            expr2_m.PushAdjoints(adjoints, coefficient * this->GetValue() * std::log(expr1_m.value()));
        }

        bool IsNonlinear()const {
            return true;
        }
        
         inline void MakeNLInteractions(bool b = false)const {
         
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
//            expr2_m.PushNLInteractions(ids);
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
            REAL_T fx = 0.0;
            ;
            REAL_T gy = expr2_m.EvaluateDerivative(b);
            REAL_T fy = 0.0;
            REAL_T gx = expr2_m.EvaluateDerivative(a);
            REAL_T gxy = expr2_m.EvaluateDerivative(a, b);
            return std::pow(f, g)*(((fx * gy) / f) + std::log(f) * gxy + (fy * gx / f) -
                    (g * fx * fy) / (f * f) + g * fxy / f) + std::pow(f, g)*(std::log(f) * gx +
                    g * fx / f)*(std::log(f) * gy + g * fy / f);

        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return std::pow(expr1_m, expr2_m.GetValue()) * std::pow(std::log(expr1_m), 3.0)
                    * (expr2_m.EvaluateDerivative(x))*(expr2_m.EvaluateDerivative(y))
                    *(expr2_m.EvaluateDerivative(z)) + std::pow(expr1_m, expr2_m.GetValue())
                    * std::pow(std::log(expr1_m), 2.0) * (expr2_m.EvaluateDerivative(x, y))
                    *(expr2_m.EvaluateDerivative(z)) + std::pow(expr1_m, expr2_m.GetValue())
                    * std::pow(std::log(expr1_m), 2.0) * (expr2_m.EvaluateDerivative(x))
                    *(expr2_m.EvaluateDerivative(y, z)) + std::pow(expr1_m, expr2_m.GetValue())
                    * std::pow(std::log(expr1_m), 2.0) * (expr2_m.EvaluateDerivative(x, z))
                    *(expr2_m.EvaluateDerivative(y)) + std::pow(expr1_m, expr2_m.GetValue())
                    * std::log(expr1_m)*(expr2_m.EvaluateDerivative(x, y, z));

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
    atl::PowScalar< REAL_T, EXPR > pow(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::PowScalar< REAL_T, EXPR > (expr.Cast(), val);
    }
    
    template <class REAL_T, class EXPR>
    inline
    atl::PowScalar< REAL_T, EXPR > square(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::PowScalar< REAL_T, EXPR > (expr.Cast(), static_cast<REAL_T>(2.0));
    }
    
    template <class REAL_T, class EXPR>
    inline
    atl::PowScalar< REAL_T, EXPR > cube(const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::PowScalar< REAL_T, EXPR > (expr.Cast(), static_cast<REAL_T>(3.0));
    }

    template <class REAL_T, class EXPR>
    inline
    atl::PowScalar< REAL_T, EXPR > operator^(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::PowScalar< REAL_T, EXPR > (expr.Cast(), val);
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
    atl::ScalarPow<REAL_T, EXPR> pow(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ScalarPow<REAL_T, EXPR > (val, expr.Cast());
    }

    template <class REAL_T, class EXPR>
    inline
    atl::ScalarPow<REAL_T, EXPR> operator^(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::ScalarPow<REAL_T, EXPR > (val, expr.Cast());
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
    atl::PowScalar< REAL_T, EXPR > pow(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::PowScalar< REAL_T, EXPR > (expr.Cast(), val);
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
    atl::ScalarPow<REAL_T, EXPR> pow(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {

        return atl::ScalarPow<REAL_T, EXPR > (val, expr.Cast());
    }

}
#endif /* POW_HPP */

