/* 
 * File:   Divide.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:54 AM
 */

#ifndef ET4AD_DIVIDE_HPP
#define	ET4AD_DIVIDE_HPP

#include "Expression.hpp"


namespace atl {

    template <class REAL_T, class LHS, class RHS>
    class Divide : public ExpressionBase<REAL_T, Divide<REAL_T, LHS, RHS> > {
    public:
        typedef REAL_T BASE_TYPE;

        typedef Subtract<REAL_T,
        Divide<REAL_T, typename LHS::DIFF_EXPRESSION, RHS>,
        Multiply<REAL_T, Divide<REAL_T, typename RHS::DIFF_EXPRESSION, Multiply<REAL_T, RHS, RHS> >,
        LHS> > DIFF_EXPRESSION;


        const LHS& lhs_m;
        const RHS& rhs_m;
        REAL_T lhs_value_m;
        REAL_T rhs_value_m;
        REAL_T mult;
        REAL_T y2;
        REAL_T y3;

        Divide(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), lhs_value_m(lhs_m.GetValue()),
        rhs_value_m(rhs_m.GetValue()), mult(1.0 / (rhs_value_m)), y2(rhs_value_m*rhs_value_m),
        y3(rhs_value_m*rhs_value_m*rhs_value_m) {


        }

        inline const REAL_T GetValue() const {
            return lhs_value_m *mult; /// rhs_value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent = true)const {
            lhs_m.PushIds(ids, include_dependent);
            rhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_value_m * lhs_m.EvaluateDerivative(id) -
                    rhs_m.EvaluateDerivative(id) * lhs_value_m) / y2; // (mult * mult);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            //            REAL_T y = rhs_value_m;
            //            REAL_T x = lhs_value_m;
            //            REAL_T fxy = lhs_m.EvaluateDerivative(a, b);
            //            REAL_T fx = lhs_m.EvaluateDerivative(a);
            REAL_T gy = rhs_m.EvaluateDerivative(b);
            //            REAL_T fy = lhs_m.EvaluateDerivative(b);
            REAL_T gx = rhs_m.EvaluateDerivative(a);
            //            REAL_T gxy = rhs_m.EvaluateDerivative(a, b);
            //            return fxy * (1.0 / y) + (fx * gy + fy * gx + gxy * x)*(-1.0 / (y2)) + (x * gx * gy) * (2.0 / (y3));
            return lhs_m.EvaluateDerivative(a, b) * (1.0 / rhs_value_m) +
                    (lhs_m.EvaluateDerivative(a) * gy + lhs_m.EvaluateDerivative(b) *
                    gx + rhs_m.EvaluateDerivative(a, b) * lhs_value_m)*(-1.0 / (y2)) +
                    (lhs_value_m * gx * gy) * (2.0 / (y3));



        }

    private:

    };

    /**
     * Operator for addition of two expression templates.
     * @param a
     * @param b
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline Divide<REAL_T, LHS, RHS> operator/(const ExpressionBase<REAL_T, LHS>& a,
            const ExpressionBase<REAL_T, RHS>& b) {
        return Divide<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class REAL_T, class LHS>
    class DivideConstant : public ExpressionBase<REAL_T, DivideConstant<REAL_T, LHS> > {
    public:
        typedef REAL_T BASE_TYPE;
        typedef MultiplyConstant<REAL_T, typename LHS::DIFF_EXPRESSION> DIFF_EXPRESSION;

        DivideConstant(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs_m.GetValue() / rhs_m) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent = true)const {
            lhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            lhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_m * lhs_m.EvaluateDerivative(id) -
                    0 * lhs_m.GetValue()) / (rhs_m * rhs_m);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T y = rhs_m;
            REAL_T x = lhs_m.GetValue();
            REAL_T fxy = lhs_m.EvaluateDerivative(a, b);
            REAL_T fx = lhs_m.EvaluateDerivative(a);
            REAL_T gy = 0;
            REAL_T fy = lhs_m.EvaluateDerivative(b);
            REAL_T gx = 0;
            REAL_T gxy = 0;
            return fxy * (1.0 / y) + (fx * gy + fy * gx + gxy * x)*(-1.0 / (y * y)) + (x * gx * gy) * (2.0 / (y * y * y));

        }



    private:

        const LHS& lhs_m;
        const REAL_T rhs_m;
        const REAL_T value_m;
    };


    /**
     * Operator for adding a expression templates to a constant .
     * @param lhs
     * @param rhs
     * @return 
     */

    /**
     * Operator for adding a expression templates to a constant .
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class LHS, class REAL_T>
    inline const DivideConstant<REAL_T, LHS > operator/(const ExpressionBase<REAL_T, LHS>& lhs,
            const REAL_T& rhs) {
        return DivideConstant<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    class ConstantDivide : public ExpressionBase<REAL_T, ConstantDivide<REAL_T, RHS> > {
    public:
        typedef REAL_T BASE_TYPE;

        ConstantDivide(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m / rhs_m.GetValue()) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent = true)const {
            rhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            rhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_m.GetValue() * 0 -
                    rhs_m.EvaluateDerivative(id) * lhs_m) / (rhs_m.GetValue() * rhs_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            REAL_T y = rhs_m.GetValue();
            REAL_T x = lhs_m;
            REAL_T fxy = 0;
            REAL_T fx = 0;
            REAL_T gy = rhs_m.EvaluateDerivative(b);
            REAL_T fy = 0;
            REAL_T gx = rhs_m.EvaluateDerivative(a);
            REAL_T gxy = rhs_m.EvaluateDerivative(a, b);
            return fxy * (1.0 / y) + (fx * gy + fy * gx + gxy * x)*(-1.0 / (y * y)) + (x * gx * gy) * (2.0 / (y * y * y));

        }




    private:

        const REAL_T lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    //    /**
    //     * Operator for adding a constant value to a expression templates.
    //     * @param lhs
    //     * @param rhs
    //     * @return 
    //     */

    template <class REAL_T, class RHS >
    inline const ConstantDivide< REAL_T, RHS > operator/(const REAL_T& lhs,
            const ExpressionBase<REAL_T, RHS>& rhs) {
        return ConstantDivide<REAL_T, RHS > (lhs, rhs.Cast());
    }

}

#endif	/* DIVIDE_HPP */

