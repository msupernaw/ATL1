/* 
 * File:   Constant.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:18 AM
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

        inline const REAL_T EvaluateDerivative(uint32_t id) const {
            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t x, uint32_t y, uint32_t z) const {
            return 0.0;
        }


    private:
        REAL_T value_m;
    };



}


#endif /* CONSTANT_HPP */

