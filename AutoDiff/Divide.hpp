/* 
 * File:   Divide.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:54 AM
 */

#ifndef ET4AD_DIVIDE_HPP
#define ET4AD_DIVIDE_HPP

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

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {
            lhs_m.PushIds(ids, include_dependent);
            rhs_m.PushIds(ids, include_dependent);
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            //            if (lhs_m.IsNonlinear() && !rhs_m.IsNonlinear()) {
            //                lhs_m.PushIds(ids);
            //                rhs_m.PushIds(ids, true);
            //            } else if (!lhs_m.IsNonlinear() && rhs_m.IsNonlinear()) {
            //                lhs_m.PushIds(ids, true);
            //                rhs_m.PushIds(ids);
            //            } else {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
            //            }
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            lhs_m.PushIds(ids);
            rhs_m.PushIds(ids);
        }

        bool IsNonFunction() const {
            if (lhs_m.IsNonFunction() || rhs_m.IsNonFunction()) {
                return true;
            } else {
                return false;
            }
        }

        bool IsNonlinear()const {
            if (lhs_m.IsNonlinear() || rhs_m.IsNonlinear()) {
                return true;
            } else {
                return false;
            }
        }

        inline void MakeNLInteractions(bool b = false)const {
            if (lhs_m.IsNonlinear() && !rhs_m.IsNonlinear()) {
                rhs_m.MakeNLInteractions(true);
            }
            if (!lhs_m.IsNonlinear() && rhs_m.IsNonlinear()) {
                lhs_m.MakeNLInteractions(true);
            }
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            if (lhs_m.IsNonlinear() && !rhs_m.IsNonlinear()) {
                rhs_m.PushIds(ids, true);
            }
            if (!lhs_m.IsNonlinear() && rhs_m.IsNonlinear()) {
                lhs_m.PushIds(ids, true);
            }
            //            lhs_m.PushNLInteractions(ids);
            //            rhs_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return (rhs_value_m * lhs_m.EvaluateDerivative(id) -
                    rhs_m.EvaluateDerivative(id) * lhs_value_m) / y2; // (mult * mult);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {

            REAL_T aa = ((2.0 * lhs_m.GetValue() * rhs_m.EvaluateDerivative(a) * rhs_m.EvaluateDerivative(b)) / (rhs_m.GetValue() * rhs_m.GetValue() * rhs_m.GetValue())); //(2*f(a,b)*('diff(g(a,b),a,1))*('diff(g(a,b),b,1)))/g(a,b)^3
            REAL_T bb = ((lhs_m.EvaluateDerivative(a) * rhs_m.EvaluateDerivative(b)) / (rhs_m.GetValue() * rhs_m.GetValue())); //(('diff(f(a,b),a,1))*('diff(g(a,b),b,1)))/g(a,b)^2
            REAL_T cc = ((lhs_m.GetValue() * rhs_m.EvaluateDerivative(a, b)) / (rhs_m.GetValue() * rhs_m.GetValue())); //(f(a,b)*('diff(g(a,b),a,1,b,1)))/g(a,b)^2
            REAL_T dd = ((lhs_m.EvaluateDerivative(b) * rhs_m.EvaluateDerivative(a)) / (rhs_m.GetValue() * rhs_m.GetValue())); //(('diff(f(a,b),b,1))*('diff(g(a,b),a,1)))/g(a,b)^2
            REAL_T ee = lhs_m.EvaluateDerivative(a, b) / rhs_m.GetValue();
            return aa - bb - cc - dd + ee;

        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return -(6 * lhs_m.GetValue()*(rhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y))*(rhs_m.EvaluateDerivative(z))) / std::pow(rhs_m.GetValue(), 4.0)+(2 * (lhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y))*(rhs_m.EvaluateDerivative(z))) / std::pow(rhs_m.GetValue(), 3.0)+
                    (2 * lhs_m.GetValue()*(rhs_m.EvaluateDerivative(x, y))*(rhs_m.EvaluateDerivative(z))) / std::pow(rhs_m.GetValue(), 3.0)+(2 * (lhs_m.EvaluateDerivative(y))*(rhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(z))) / std::pow(rhs_m.GetValue(), 3.0)-
                    ((lhs_m.EvaluateDerivative(x, y))*(rhs_m.EvaluateDerivative(z))) / std::pow(rhs_m.GetValue(), 2.0)+(2 * lhs_m.GetValue()*(rhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y, z))) / std::pow(rhs_m.GetValue(), 3.0)-((lhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y, z))) / std::pow(rhs_m.GetValue(), 2.0)+
                    (2 * lhs_m.GetValue()*(rhs_m.EvaluateDerivative(x, z))*(rhs_m.EvaluateDerivative(y))) / std::pow(rhs_m.GetValue(), 3.0)+(2 * (lhs_m.EvaluateDerivative(z))*(rhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y))) / std::pow(rhs_m.GetValue(), 3.0)-
                    ((lhs_m.EvaluateDerivative(x, z))*(rhs_m.EvaluateDerivative(y))) / std::pow(rhs_m.GetValue(), 2.0)-((lhs_m.EvaluateDerivative(y))*(rhs_m.EvaluateDerivative(x, z))) / std::pow(rhs_m.GetValue(), 2.0)-(lhs_m.GetValue()*(rhs_m.EvaluateDerivative(x, y, z))) / std::pow(rhs_m.GetValue(), 2.0)-
                    ((lhs_m.EvaluateDerivative(z))*(rhs_m.EvaluateDerivative(x, y))) / std::pow(rhs_m.GetValue(), 2.0)-((lhs_m.EvaluateDerivative(y, z))*(rhs_m.EvaluateDerivative(x))) / std::pow(rhs_m.GetValue(), 2.0) + lhs_m.EvaluateDerivative(x, y, z) / rhs_m.GetValue()
                    ;

        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicDivide<REAL_T>(lhs_m.GetDynamicExpession(), rhs_m.GetDynamicExpession());
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
    class DivideScalar : public ExpressionBase<REAL_T, DivideScalar<REAL_T, LHS> > {
    public:
        typedef REAL_T BASE_TYPE;
        //        typedef MultiplyScalar<REAL_T, typename LHS::DIFF_EXPRESSION> DIFF_EXPRESSION;

        DivideScalar(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs), value_m(lhs_m.GetValue() / rhs_m) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
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

        bool IsNonFunction() const {
            if (lhs_m.IsNonFunction()) {
                return true;
            } else {
                return false;
            }
        }

        bool IsNonlinear()const {
            if (lhs_m.IsNonlinear()) {
                return true;
            } else {
                return false;
            }
        }

        inline void MakeNLInteractions(bool b = false)const {

            if (!lhs_m.IsNonlinear()) {
                lhs_m.MakeNLInteractions(b);
            }
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            if (!lhs_m.IsNonlinear())
                lhs_m.PushNLInteractions(ids);
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

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return lhs_m.EvaluateDerivative(x, y, z)
                    / rhs_m;
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicDivide<REAL_T>(lhs_m.GetDynamicExpession(), new atl::DynamicScalar<REAL_T>(rhs_m));
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
    inline const DivideScalar<REAL_T, LHS > operator/(const ExpressionBase<REAL_T, LHS>& lhs,
            const REAL_T& rhs) {
        return DivideScalar<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    class ScalarDivide : public ExpressionBase<REAL_T, ScalarDivide<REAL_T, RHS> > {
    public:
        typedef REAL_T BASE_TYPE;

        ScalarDivide(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m / rhs_m.GetValue()) {


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
            rhs_m.PushAdjoints(adjoints, -1.0 * coefficient * this->GetValue() / rhs_m);
        }

        bool IsNonFunction() const {
            if (rhs_m.IsNonFunction()) {
                return true;
            } else {
                return false;
            }
        }

        bool IsNonlinear()const {
            if (rhs_m.IsNonlinear()) {
                return true;
            } else {
                return false;
            }
        }

        inline void MakeNLInteractions(bool b = false)const {
            if (!rhs_m.IsNonlinear()) {
                rhs_m.MakeNLInteractions(b);
            }
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            if (!rhs_m.IsNonlinear())
                rhs_m.PushNLInteractions(ids);
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

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return -(6 * lhs_m * (rhs_m.EvaluateDerivative(x))*
                    (rhs_m.EvaluateDerivative(y))*(rhs_m.EvaluateDerivative(z)))
                    / std::pow(rhs_m.GetValue(), 4.0)+(2 * lhs_m * (rhs_m.EvaluateDerivative(x, y))
                    *(rhs_m.EvaluateDerivative(z))) / std::pow(rhs_m.GetValue(), 3.0)+
                    (2 * lhs_m * (rhs_m.EvaluateDerivative(x))*(rhs_m.EvaluateDerivative(y, z)))
                    / std::pow(rhs_m.GetValue(), 3.0)+(2 * lhs_m * (rhs_m.EvaluateDerivative(x, z))
                    *(rhs_m.EvaluateDerivative(y))) / std::pow(rhs_m.GetValue(), 3.0)-
                    (lhs_m * (rhs_m.EvaluateDerivative(x, y, z))) / std::pow(rhs_m.GetValue(), 2.0);
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicDivide<REAL_T>(new atl::DynamicScalar<REAL_T>(lhs_m), rhs_m.GetDynamicExpession());
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
    inline const ScalarDivide< REAL_T, RHS > operator/(const REAL_T& lhs,
            const ExpressionBase<REAL_T, RHS>& rhs) {
        return ScalarDivide<REAL_T, RHS > (lhs, rhs.Cast());
    }

}

#endif /* DIVIDE_HPP */

