/* 
 * File:   ArrayMultiply.hpp
 * Author: matthewsupernaw
 *
 * Created on August 25, 2014, 1:29 PM
 */

#ifndef ARRAYMULTIPLY_HPP
#define	ARRAYMULTIPLY_HPP
#include "ArrayExpressionBase.hpp"
#include "ArrayTraits.hpp"

namespace atl {

    template< class LHS, class RHS>
    struct ArrayMultiply : ArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, ArrayMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit ArrayMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            for (int i = 0; i < 7; i++) {
                assert(lhs_m.Size(i) == rhs_m.Size(i));
            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return lhs_m(i) * rhs_m(i);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) * rhs_m(i, j);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return lhs_m(i, j, k) * rhs_m(i, j, k);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return lhs_m(i, j, k, l) * rhs_m(i, j, k, l);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return lhs_m(i, j, k, j, m) * rhs_m(i, j, k, l, m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return lhs_m(i, j, k, j, m, n) * rhs_m(i, j, k, l, m, n);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return lhs_m(i, j, k, j, m, n, o) * rhs_m(i, j, k, l, m, n, o);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }


    };
    //

    template< class LHS, class T >
    struct ArrayMultiplyScalar : ArrayExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, ArrayMultiplyScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;

        inline explicit ArrayMultiplyScalar(const ArrayExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return lhs_m(i) * rhs_m;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) * rhs_m;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return lhs_m(i, j, k) * rhs_m;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return lhs_m(i, j, k, l) * rhs_m;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return lhs_m(i, j, k, j, m) * rhs_m;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return lhs_m(i, j, k, j, m, n) * rhs_m;
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return lhs_m(i, j, k, j, m, n, o) * rhs_m;
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }


    };

    template<class T, class RHS>
    struct ArrayScalarMultiply : ArrayExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ArrayScalarMultiply<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;


        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ArrayScalarMultiply(const T& lhs, const ArrayExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return lhs_m * rhs_m(i);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m * rhs_m(i, j);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return lhs_m * rhs_m(i, j, k);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return lhs_m * rhs_m(i, j, k, l);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return lhs_m * rhs_m(i, j, k, j, m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return lhs_m * rhs_m(i, j, k, j, m, n);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return lhs_m * rhs_m(i, j, k, j, m, n, o);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            rhs_m.ExpressionLength(length);
        }

    };
    //
    //typename et4ad::promote_trait<typename LHS::RET_TYPE , typename RHS::RET_TYPE >::return_type

    template <class LHS, class RHS>
    inline const ArrayMultiply< LHS, RHS> operator*(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return ArrayMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_ARRAY_MULTIPLY_SCALAR(TYPE) \
    template< class LHS>      \
        inline const ArrayMultiplyScalar<LHS,TYPE> \
    operator*(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return ArrayMultiplyScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_ARRAY_MULTIPLY_SCALAR(short)
    ATL_ARRAY_MULTIPLY_SCALAR(unsigned short)
    ATL_ARRAY_MULTIPLY_SCALAR(int)
    ATL_ARRAY_MULTIPLY_SCALAR(unsigned int)
    ATL_ARRAY_MULTIPLY_SCALAR(long)
    ATL_ARRAY_MULTIPLY_SCALAR(unsigned long)
    ATL_ARRAY_MULTIPLY_SCALAR(float)
    ATL_ARRAY_MULTIPLY_SCALAR(double)
    ATL_ARRAY_MULTIPLY_SCALAR(long double)
    ATL_ARRAY_MULTIPLY_SCALAR(atl::Variable<float>)
    ATL_ARRAY_MULTIPLY_SCALAR(atl::Variable<double>)
    ATL_ARRAY_MULTIPLY_SCALAR(atl::Variable<long double>)


#define ATL_MULTIPLY_SCALAR_ARRAY(TYPE) \
    template< class RHS>      \
        inline const ArrayScalarMultiply<TYPE,RHS> \
    operator*(const TYPE& a, const ArrayExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ArrayScalarMultiply<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_MULTIPLY_SCALAR_ARRAY(short)
    ATL_MULTIPLY_SCALAR_ARRAY(unsigned short)
    ATL_MULTIPLY_SCALAR_ARRAY(int)
    ATL_MULTIPLY_SCALAR_ARRAY(unsigned int)
    ATL_MULTIPLY_SCALAR_ARRAY(long)
    ATL_MULTIPLY_SCALAR_ARRAY(unsigned long)
    ATL_MULTIPLY_SCALAR_ARRAY(float)
    ATL_MULTIPLY_SCALAR_ARRAY(double)
    ATL_MULTIPLY_SCALAR_ARRAY(long double)
    ATL_MULTIPLY_SCALAR_ARRAY(atl::Variable<float>)
    ATL_MULTIPLY_SCALAR_ARRAY(atl::Variable<double>)
    ATL_MULTIPLY_SCALAR_ARRAY(atl::Variable<long double>)


}






#endif	/* ARRAYMULTIPLY_HPP */

