/* 
 * File:   Multiply.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 12:00 PM
 */

#ifndef ET4AD_MULTIPLY_HPP
#define ET4AD_MULTIPLY_HPP

#include "Expression.hpp"
#include "Scalar.hpp"
#include "Log.hpp"



namespace atl {

    template <class REAL_T, class LHS, class RHS>
    struct Multiply : public ExpressionBase<REAL_T, Multiply<REAL_T, LHS, RHS> > {
        const LHS& lhs_m;
        const RHS& rhs_m;

        Multiply(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {

        }

        inline const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m.GetValue();
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            lhs_m.PushIds(ids, include_dependent);
            rhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            lhs_m.PushAdjoints(adjoints, coefficient * rhs_m.GetValue());
            rhs_m.PushAdjoints(adjoints, coefficient * lhs_m.GetValue());
        }

        inline const REAL_T EvaluateDerivative(uint32_t id) const {
            return (lhs_m.GetValue() * rhs_m.EvaluateDerivative(id) +
                    lhs_m.EvaluateDerivative(id) * rhs_m.GetValue());
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return lhs_m.EvaluateDerivative(a) * rhs_m.EvaluateDerivative(b) + lhs_m.GetValue() * rhs_m.EvaluateDerivative(a, b) +
                    lhs_m.EvaluateDerivative(b) * rhs_m.EvaluateDerivative(a) + rhs_m.GetValue() * lhs_m.EvaluateDerivative(a, b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return (lhs_m.EvaluateDerivative(x, y))*(rhs_m.EvaluateDerivative(z))+
                    (lhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y, z))
                    +(lhs_m.EvaluateDerivative(x, z))*(rhs_m.EvaluateDerivative(y))
                    +(lhs_m.EvaluateDerivative(y))*(rhs_m.EvaluateDerivative(x, z))
                    + lhs_m.GetValue()*(rhs_m.EvaluateDerivative(x, y, z))+
                    (lhs_m.EvaluateDerivative(z))*(rhs_m.EvaluateDerivative(x, y))+
                    (lhs_m.EvaluateDerivative(y, z))*(rhs_m.EvaluateDerivative(x))
                    + rhs_m.GetValue()*(lhs_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicMultiply<REAL_T>(lhs_m.GetDynamicExpession(), rhs_m.GetDynamicExpession());
        }

    };

    /**
     * Operator for addition of two expression templates.
     * @param a
     * @param b
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline const Multiply<REAL_T, LHS, RHS> operator*(const ExpressionBase<REAL_T, LHS>& a,
            const ExpressionBase<REAL_T, RHS>& b) {
        return Multiply<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class REAL_T, class LHS>
    inline const Multiply<REAL_T, LHS, LHS> square(const ExpressionBase<REAL_T, LHS>& a) {
        return Multiply<REAL_T, LHS, LHS > (a.Cast(), a.Cast());
    }

    template <class REAL_T, class LHS>
    struct MultiplyScalar : public ExpressionBase<REAL_T, MultiplyScalar<REAL_T, LHS> > {
        typedef REAL_T BASE_TYPE;

        MultiplyScalar(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {

        }

        inline const REAL_T GetValue() const {
            return lhs_m.GetValue() * rhs_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            lhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            lhs_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            lhs_m.PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            lhs_m.PushAdjoints(adjoints, coefficient * rhs_m);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_m * lhs_m.EvaluateDerivative(id));
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return lhs_m.EvaluateDerivative(a, b) * rhs_m;
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return rhs_m * (lhs_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicMultiply<REAL_T>(lhs_m.GetDynamicExpession(), atl::DynamicScalar<REAL_T>(rhs_m));
        }


        const LHS& lhs_m;
        const REAL_T rhs_m;
    };

    /**
     * Operator for adding a expression templates to a constant .
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class LHS, class REAL_T>
    inline const MultiplyScalar<REAL_T, LHS > operator*(const ExpressionBase<REAL_T, LHS>& lhs,
            const REAL_T& rhs) {

        return MultiplyScalar<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    struct ScalarMultiply : public ExpressionBase<REAL_T, ScalarMultiply<REAL_T, RHS> > {
        typedef REAL_T BASE_TYPE;

        ScalarMultiply(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m * rhs_m.GetValue()) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            rhs_m.VariableCount(count);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            rhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            rhs_m.PushIds(ids);
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            rhs_m.PushIds(ids);
        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {
            rhs_m.PushAdjoints(adjoints, coefficient * lhs_m);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return rhs_m.EvaluateDerivative(id) * lhs_m;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return (rhs_m.EvaluateDerivative(a, b) * lhs_m);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return lhs_m * (rhs_m.EvaluateDerivative(x, y, z));
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicMultiply<REAL_T>(new atl::DynamicScalar<REAL_T>(lhs_m), rhs_m.GetDynamicExpession());
        }



        const REAL_T lhs_m;
        const RHS& rhs_m;
        const REAL_T value_m;
    };

    /**
     * Operator for adding a constant value to a expression templates.
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS >
    inline const ScalarMultiply< REAL_T, RHS > operator*(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs) {
        return ScalarMultiply<REAL_T, RHS > (lhs, rhs.Cast());
    }

}


#endif /* MULTIPLY_HPP */

