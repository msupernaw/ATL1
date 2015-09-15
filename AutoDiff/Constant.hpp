/* 
 * File:   Constant.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:18 AM
 */

#ifndef CONSTANT_HPP
#define	CONSTANT_HPP

#include "Expression.hpp"


namespace atl {

    template<class REAL_T>
    class Constant : public atl::ExpressionBase<REAL_T, Constant<REAL_T> > {
    public:
        typedef REAL_T BASE_TYPE;

        Constant() : value_m(static_cast<REAL_T> (0.0)) {

        }

        Constant(const REAL_T & value) : value_m(value) {

        }

        inline const REAL_T GetValue() const {
            return value_m;
        }

        inline void VariableCount(uint32_t& count) const {
        }

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
        }

        inline REAL_T EvaluateDerivative(uint32_t id) const {
            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return 0.0;

        }

    private:
        REAL_T value_m;
    };



}


#endif	/* CONSTANT_HPP */

