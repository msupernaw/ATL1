/* 
 * File:   Subtraction.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:52 AM
 */

#ifndef ET4AD_SUBTRACT_HPP
#define	ET4AD_SUBTRACT_HPP

#include "Expression.hpp"
//#include "Constant.hpp"
namespace atl {

    template <class REAL_T, class LHS, class RHS>
    struct Subtract : public ExpressionBase<REAL_T, Subtract<REAL_T, LHS, RHS> > {
        typedef REAL_T BASE_TYPE;

        Subtract(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), value_m(lhs_m.GetValue() - rhs_m.GetValue()) {

        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<double>* >& ids)const {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return lhs_m.EvaluateDerivative(id) - rhs_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return lhs_m.EvaluateDerivative(a, b) - rhs_m.EvaluateDerivative(a, b);
        }



        const LHS& lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for subtracting an expression templates from another.
     * @param a
     * @param b
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline Subtract<REAL_T, LHS, RHS> operator-(const ExpressionBase<REAL_T, LHS>& a,
            const ExpressionBase<REAL_T, RHS>& b) {
        return Subtract<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class REAL_T, class LHS>
    struct SubtractConstant : public ExpressionBase<REAL_T, SubtractConstant<REAL_T, LHS> > {
        typedef REAL_T BASE_TYPE;
        const LHS& lhs_m;
        const REAL_T rhs_m;
        const REAL_T value_m;

        SubtractConstant(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs.GetValue() - rhs) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            lhs_m.PushIds(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return lhs_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return lhs_m.EvaluateDerivative(a, b);
        }



    };

    /**
     * Operator for subtracting a constant from a expression templates .
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class LHS, class REAL_T>
    inline const SubtractConstant<REAL_T, LHS> operator-(const ExpressionBase<REAL_T, LHS>& lhs,
            const REAL_T& rhs) {
        return SubtractConstant<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    struct ConstantSubtract : public ExpressionBase<REAL_T, ConstantSubtract<REAL_T, RHS> > {
        typedef REAL_T BASE_TYPE;

        ConstantSubtract(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m - rhs_m.GetValue()) {


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
            return -1.0 * rhs_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return -1.0 * rhs_m.EvaluateDerivative(a, b);
        }

        const REAL_T lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    } __attribute__((packed));

    /**
     * Operator for subtracting  a expression templates from a constant.
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS >
    inline const ConstantSubtract< REAL_T, RHS > operator-(const REAL_T& lhs,
            const ExpressionBase<REAL_T, RHS> &rhs) {
        return ConstantSubtract<REAL_T, RHS > (lhs, rhs.Cast());
    }

}


#endif	/* SUBTRACTION_HPP */

