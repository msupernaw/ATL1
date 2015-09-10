/* 
 * File:   Multiply.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:00 PM
 */

#ifndef ET4AD_MULTIPLY_HPP
#define	ET4AD_MULTIPLY_HPP

#include "Expression.hpp"
#include "Constant.hpp"
#include "Log.hpp"



namespace atl {

    template <class REAL_T, class LHS, class RHS>
    struct Multiply : public ExpressionBase<REAL_T, Multiply<REAL_T, LHS, RHS> > {
        typedef REAL_T BASE_TYPE;
        const LHS& lhs_m;
        const RHS& rhs_m;
        const REAL_T lhs_value_m;
        const REAL_T rhs_value_m;
        mutable REAL_T coeff;

        Multiply(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), lhs_value_m(lhs_m.GetValue()), rhs_value_m(rhs_m.GetValue())/*, value_m(lhs_value_m * rhs_value_m) */ {

        }

        inline const REAL_T GetValue() const {
            return lhs_value_m * rhs_value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_value_m * lhs_m.EvaluateDerivative(id) +
                    rhs_m.EvaluateDerivative(id) * lhs_value_m);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            //
            //            REAL_T fxy = lhs_m.EvaluateDerivative(a, b);
            //            REAL_T g = rhs_value_m;
            //            REAL_T f = lhs_value_m;
            //            REAL_T fx = lhs_m.EvaluateDerivative(a);
            //            REAL_T gy = rhs_m.EvaluateDerivative(b);
            //            REAL_T fy = lhs_m.EvaluateDerivative(b);
            //            REAL_T gx = rhs_m.EvaluateDerivative(a);
            //            REAL_T gxy = rhs_m.EvaluateDerivative(a, b);
            //
            //            return fxy * g + fx * gy + fy * gx + gxy*f;

            return lhs_m.EvaluateDerivative(a, b) * rhs_value_m +
                    lhs_m.EvaluateDerivative(a) * rhs_m.EvaluateDerivative(b) +
                    lhs_m.EvaluateDerivative(b) * rhs_m.EvaluateDerivative(a) +
                    rhs_m.EvaluateDerivative(a, b) * lhs_value_m;
        }




    } __attribute__((packed));

    /**
     * Operator for addition of two expression templates.
     * @param a
     * @param b
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline Multiply<REAL_T, LHS, RHS> operator*(const ExpressionBase<REAL_T, LHS>& a,
            const ExpressionBase<REAL_T, RHS>& b) {
        return Multiply<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class REAL_T, class LHS>
    struct MultiplyConstant : public ExpressionBase<REAL_T, MultiplyConstant<REAL_T, LHS> > {
        typedef REAL_T BASE_TYPE;

        MultiplyConstant(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {

        }

        inline const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            lhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_m * lhs_m.EvaluateDerivative(id));
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {

            REAL_T fxy = lhs_m.EvaluateDerivative(a, b);
            REAL_T g = rhs_m;
//            REAL_T f = lhs_m.GetValue();
//            REAL_T fx = lhs_m.EvaluateDerivative(a);
//            REAL_T gy = 0;
//            REAL_T fy = lhs_m.EvaluateDerivative(b);
//            REAL_T gx = 0;
//            REAL_T gxy = 0;

            return fxy * g;
        }

        const LHS& lhs_m;
        const REAL_T rhs_m;
    } __attribute__((packed));

    /**
     * Operator for adding a expression templates to a constant .
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class LHS, class REAL_T>
    inline const MultiplyConstant<REAL_T, LHS > operator*(const ExpressionBase<REAL_T, LHS>& lhs,
            const REAL_T& rhs) {

        return MultiplyConstant<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    struct ConstantMultiply : public ExpressionBase<REAL_T, ConstantMultiply<REAL_T, RHS> > {
        typedef REAL_T BASE_TYPE;

        ConstantMultiply(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m * rhs_m.GetValue()) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            rhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return rhs_m.EvaluateDerivative(id) * lhs_m;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return (rhs_m.EvaluateDerivative(a, b) * lhs_m);
        }


        const REAL_T lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    }
    __attribute__((packed));

    /**
     * Operator for adding a constant value to a expression templates.
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS >
    inline const ConstantMultiply< REAL_T, RHS > operator*(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs) {
        return ConstantMultiply<REAL_T, RHS > (lhs, rhs.Cast());
    }

}


#endif	/* MULTIPLY_HPP */

