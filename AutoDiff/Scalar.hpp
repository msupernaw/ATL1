/* 
 * File:   Constant.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:18 AM
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


#ifndef CONSTANT_HPP
#define CONSTANT_HPP

#include "Expression.hpp"


namespace atl {

    template<class REAL_T>
    class Scalar : public atl::ExpressionBase<REAL_T, Scalar<REAL_T> > {
    public:
        typedef REAL_T BASE_TYPE;

        Scalar(REAL_T value) : value_m(value) {

        }

        Scalar(const Scalar<REAL_T>& other) :
        value_m(other.value_m) {
        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {

        }

        operator REAL_T() {
            return this->GetValue();
        }

        operator REAL_T()const {
            return this->GetValue();
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent)const {

        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {

        }

        inline void PushIds(IDSet<uint32_t >& ids)const {

        }

        inline void PushAdjoints(std::vector<std::pair<atl::VariableInfo<REAL_T>*, REAL_T> >& adjoints, REAL_T coefficient = 1.0) const {

        }

              bool IsNonFunction()const {
            return false;
        }
        
        bool IsNonlinear()const {
            return false;
        }
        
         inline void MakeNLInteractions(bool b = false)const {
         
        }

        inline void PushNLInteractions(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
   
        }

        inline const REAL_T EvaluateDerivative(uint32_t id) const {
            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return 0.0;
        }
        
        inline atl::DynamicExpression<REAL_T>* GetDynamicExpession() const {
            return new atl::DynamicScalar<REAL_T>(value_m);
        }


    private:
        REAL_T value_m;
    };



}


#endif /* CONSTANT_HPP */

