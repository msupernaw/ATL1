/* 
 * File:   ATan2.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:29 PM
 */

#ifndef ET4AD_ATAN2_HPP
#define	ET4AD_ATAN2_HPP

#include <cmath>
#include "Expression.hpp"

namespace atl {

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both arguments are expression templates. 
     * 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    struct ATan2;

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both the first argument is an expression template, the second
     * is constant. 
     * 
     */
    template <class REAL_T, class EXPR1>
    struct ATan2Constant;

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both arguments are expression templates. 
     * 
     */
    template <class REAL_T, class EXPR2>
    struct ConstantATan2;

}



namespace std {

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    atl::ATan2<REAL_T, EXPR1, EXPR2> atan2(const atl::ExpressionBase<REAL_T, EXPR1>& expr1,
            const atl::ExpressionBase<REAL_T, EXPR2>& expr2);

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ATan2Constant<REAL_T, EXPR > atan2(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val);

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ConstantATan2<REAL_T, EXPR> atan2(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr);



}

namespace atl {

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both arguments are expression templates. 
     * 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    struct ATan2 : public ExpressionBase<REAL_T, ATan2<REAL_T, EXPR1, EXPR2> > {
        typedef REAL_T BASE_TYPE;

        ATan2(const ExpressionBase<REAL_T, EXPR1>& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2.Cast()), value1_m(expr1.GetValue()), value2_m(expr2.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan2(value1_m, value2_m);
        }

        inline void VariableCount(uint32_t& count) const {
            expr1_m.VariableCount(count);
            expr2_m.VariableCount(count);
        }

        //        operator REAL_T() {
        //            return this->GetValue();
        //        }
        //
        //        operator REAL_T()const {
        //            return this->GetValue();
        //        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            return 0.0;
        }

      

        const EXPR1& expr1_m;
        const EXPR2& expr2_m;
        REAL_T value1_m;
        REAL_T value2_m;
    };

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both the first argument is an expression template, the second
     * is constant. 
     * 
     */
    template <class REAL_T, class EXPR1>
    struct ATan2Constant : public ExpressionBase<REAL_T, ATan2Constant<REAL_T, EXPR1> > {
        typedef REAL_T BASE_TYPE;

        ATan2Constant(const ExpressionBase<REAL_T, EXPR1>& expr1, const REAL_T & expr2)
        : expr1_m(expr1.Cast()), expr2_m(expr2), value1_m(expr1_m.GetValue()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan2(value1_m, expr2_m);
        }

        inline void VariableCount(uint32_t& count) const {
            expr1_m.VariableCount(count);
        }

        //        operator REAL_T() {
        //            return this->GetValue();
        //        }
        //
        //        operator REAL_T()const {
        //            return this->GetValue();
        //        }

        inline void PushStackEntry(Entry& entry, REAL_T coefficient = 1.0) const {
            std::cout << __func__ << "not yet implemented!";
            exit(0);
        }

    

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            return 0.0;
        }



        const EXPR1& expr1_m;
        const REAL_T& expr2_m;
        REAL_T value1_m;
    } ;

    /**
     * Expression template for computing the two argument inverse tangent,
     * where both arguments are expression templates. 
     * 
     */
    template <class REAL_T, class EXPR2>
    struct ConstantATan2 : public ExpressionBase<REAL_T, ConstantATan2<REAL_T, EXPR2> > {
        typedef REAL_T BASE_TYPE;

        ConstantATan2(const REAL_T& expr1, const ExpressionBase<REAL_T, EXPR2>& expr2)
        : expr1_m(expr1), expr2_m(expr2.Cast()) {
        }

        inline const REAL_T GetValue() const {
            return std::atan2(expr1_m, expr2_m.GetValue());
        }

        inline void VariableCount(uint32_t& count) const {
            expr2_m.VariableCount(count);
        }

        //        operator REAL_T() {
        //            return this->GetValue();
        //        }
        //
        //        operator REAL_T()const {
        //            return this->GetValue();
        //        }

        inline void PushStackEntry(Entry& entry, REAL_T coefficient = 1.0) const {
            std::cout << __func__ << "not yet implemented!";
            exit(0);
        }

      

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            return 0.0;
        }


        const REAL_T& expr1_m;
        const EXPR2& expr2_m;
    };

}

namespace std {

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR1, class EXPR2>
    inline
    atl::ATan2<REAL_T, EXPR1, EXPR2> atan2(const atl::ExpressionBase<REAL_T, EXPR1>& expr1,
            const atl::ExpressionBase<REAL_T, EXPR2>& expr2) {
        return atl::ATan2<REAL_T, EXPR1, EXPR2 > (expr1.Cast(), expr2.Cast());
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param val
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ATan2Constant<REAL_T, EXPR > atan2(const atl::ExpressionBase<REAL_T, EXPR>& expr,
            const REAL_T& val) {
        return atl::ATan2Constant<REAL_T, EXPR> (expr.Cast(), val);
    }

    /**
     * Override for the atan2 function in namespace std.
     * 
     * @param expr1
     * @param expr2
     * @return 
     */
    template <class REAL_T, class EXPR>
    inline
    atl::ConstantATan2<REAL_T, EXPR> atan2(const REAL_T& val,
            const atl::ExpressionBase<REAL_T, EXPR>& expr) {
        return atl::ConstantATan2<REAL_T, EXPR > (val, expr.Cast());
    }

}
#endif	/* ATAN2_HPP */

