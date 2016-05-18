/* 
 * File:   Addition.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:36 AM
 */

/**
 *
 * @author  Matthew R. Supernaw
 *
 * Public Domain Notice
 * National Oceanic And Atmospheric Administration
 *
 * This software is a "United States Government Work" under the terms of the
 * United States Copyright Act.  It was written as part of the author's official
 * duties as a United States Government employee and thus cannot be copyrighted.
 * This software is freely available to the public for use. The National Oceanic
 * And Atmospheric Administration and the U.S. Government have not placed any
 * restriction on its use or reproduction.  Although all reasonable efforts have
 * been taken to ensure the accuracy and reliability of the software and data,
 * the National Oceanic And Atmospheric Administration and the U.S. Government
 * do not and cannot warrant the performance warrant the performance or results
 * that may be obtained by using this  software or data. The National Oceanic
 * And Atmospheric Administration and the U.S. Government disclaim all
 * warranties, express or implied, including warranties of performance,
 * merchantability or fitness for any particular purpose.
 *
 * Please cite the author(s) in any work or product based on this material.
 *
 */


#ifndef ET4AD_ADD_HPP
#define ET4AD_ADD_HPP

#include "Expression.hpp"


namespace atl {

    template <class REAL_T, class LHS, class RHS>
    struct Add : public ExpressionBase<REAL_T, Add<REAL_T, LHS, RHS> > {
        typedef REAL_T BASE_TYPE;

        Add(const ExpressionBase<REAL_T, LHS>& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()), value_m(lhs_m.GetValue() + rhs_m.GetValue()) {

            //            std::cout<<lhs.GetId()<<" - "<<rhs.GetId()<<"\n";
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            lhs_m.VariableCount(count);
            rhs_m.VariableCount(count);
        }

        operator REAL_T() {
            return this->GetValue();
        }

        operator REAL_T()const {
            return this->GetValue();
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

        bool IsNonFunction()const {
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
            //            if (lhs_m.IsNonlinear() && !rhs_m.IsNonlinear()) {
            //                rhs_m.MakeNLInteractions(true);
            //            }
            //            if (!lhs_m.IsNonlinear() && rhs_m.IsNonlinear()) {
            //                lhs_m.MakeNLInteractions(true);
            //            }
            lhs_m.MakeNLInteractions(b);
            rhs_m.MakeNLInteractions(b);
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            //            if (lhs_m.IsNonlinear() && !rhs_m.IsNonlinear()) {
            //                rhs_m.PushIds(ids, true);
            //            }
            //            if (!lhs_m.IsNonlinear() && rhs_m.IsNonlinear()) {
            //                lhs_m.PushIds(ids, true);
            //            }
            //            if(!lhs_m.IsNonlinear())
            lhs_m.PushNLInteractions(ids);
            rhs_m.PushNLInteractions(ids);
        }

        inline const REAL_T EvaluateDerivative(uint32_t id) const {
            return lhs_m.EvaluateDerivative(id) + rhs_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return lhs_m.EvaluateDerivative(a, b) + rhs_m.EvaluateDerivative(a, b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return lhs_m.EvaluateDerivative(x, y, z) + rhs_m.EvaluateDerivative(x, y, z);
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicAdd<REAL_T>(lhs_m.GetDynamicExpession(), rhs_m.GetDynamicExpession());
        }


        const LHS& lhs_m;
        const RHS& rhs_m;
        REAL_T value_m;
    };

    /**
     * Operator for addition of two expression templates.
     * @param a
     * @param b
     * @return 
     */
    template <class REAL_T, class LHS, class RHS>
    inline Add<REAL_T, LHS, RHS> operator+(const ExpressionBase<REAL_T, LHS>& a,
            const ExpressionBase<REAL_T, RHS>& b) {
        return Add<REAL_T, LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class REAL_T, class LHS>
    struct AddScalar : public ExpressionBase<REAL_T, AddScalar<REAL_T, LHS> > {
        typedef REAL_T BASE_TYPE;

        AddScalar(const ExpressionBase<REAL_T, LHS>& lhs, const REAL_T & rhs)
        : lhs_m(lhs.Cast()), rhs_m(rhs) {


        }

        inline const REAL_T GetValue() const {
            return lhs_m.GetValue() + rhs_m;
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
            lhs_m.PushAdjoints(adjoints, coefficient);
        }

        bool IsNonFunction()const {
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

            lhs_m.MakeNLInteractions(b);

        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            lhs_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return lhs_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return lhs_m.EvaluateDerivative(a, b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return lhs_m.EvaluateDerivative(x, y, z);
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicAdd<REAL_T>(lhs_m.GetDynamicExpession(), new atl::DynamicScalar<REAL_T>(rhs_m));
        }

        const LHS& lhs_m;
        REAL_T rhs_m;
        //        const REAL_T value_m;
    };

    /**
     * Operator for adding a expression templates to a constant .
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class LHS, class REAL_T>
    inline const AddScalar<REAL_T, LHS> operator+(const ExpressionBase<REAL_T, LHS>& lhs,
            const REAL_T& rhs) {
        return AddScalar<REAL_T, LHS > (lhs.Cast(), rhs);
    }

    template <class REAL_T, class RHS>
    struct ScalarAdd : public ExpressionBase<REAL_T, ScalarAdd<REAL_T, RHS> > {
        typedef REAL_T BASE_TYPE;

        ScalarAdd(const REAL_T& lhs, const ExpressionBase<REAL_T, RHS>& rhs)
        : lhs_m(lhs), rhs_m(rhs.Cast()), value_m(lhs_m + rhs_m.GetValue()) {


        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
            rhs_m.VariableCount(count);
        }

        inline void PushStackEntry(Entry& entry, REAL_T coefficient = 1.00000e+0) const {
            rhs_m.PushStackEntry(entry, coefficient);
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
            rhs_m.PushAdjoints(adjoints, coefficient);
        }

        bool IsNonFunction()const {
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
            rhs_m.MakeNLInteractions(b);

        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
            rhs_m.PushNLInteractions(ids);
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return rhs_m.EvaluateDerivative(id);
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return rhs_m.EvaluateDerivative(a, b);
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return rhs_m.EvaluateDerivative(x, y, z);
        }

        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicAdd<REAL_T>(new atl::DynamicScalar<REAL_T>(lhs_m), rhs_m.GetDynamicExpession());
        }

        REAL_T lhs_m;
        const RHS& rhs_m;
        REAL_T value_m;
    };

    /**
     * Operator for adding a constant value to a expression templates.
     * @param lhs
     * @param rhs
     * @return 
     */
    template <class REAL_T, class RHS >
    inline const ScalarAdd< REAL_T, RHS > operator+(const REAL_T& lhs,
            const ExpressionBase<REAL_T, RHS> &rhs) {
        return ScalarAdd<REAL_T, RHS > (lhs, rhs.Cast());
    }

}




#endif /* ADDITION_HPP */

