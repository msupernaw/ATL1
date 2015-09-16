/* 
 * File:   MatrixMultiply.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 3:24 PM
 */

#ifndef MATRIXMULTIPLY_HPP
#define	MATRIXMULTIPLY_HPP
#include "ContainerDefs.hpp"

#ifdef ATL_HAS_SSE
#include <emmintrin.h>
#endif
#include <iostream>
#include <cassert>
#include "MatrixExpressionBase.hpp"
#include "ArrayTraits.hpp"
#include "../AutoDiff/Variable.hpp"
namespace atl {

    template< class LHS, class RHS>
    struct MatrixMultiply : MatrixExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;
        typedef typename atl::intrinsic_trait<typename LHS::RET_TYPE> INTR_TRAITL;
        typedef typename atl::intrinsic_trait<typename RHS::RET_TYPE> INTR_TRAITR;

        const LHS& lhs_m;
        const RHS& rhs_m;
        size_t lrows;
        size_t lcols;
        size_t end;

        inline explicit MatrixMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
            lrows = lhs_m.Size(0);
            lcols = lhs_m.Size(1);
            end = (((lcols - 1UL) & size_t(-2)) + 1UL);
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            //                        for (int i = 0; i < 2; i++) {
            //                            assert(lhs_m.Size(i) == rhs_m.Size(i));
            //                        }
            assert(lhs_m.Size(1) == rhs_m.Size(0));
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {

            if (lhs_m.Size(1) == rhs_m.Size(0)) {

                switch (dimension) {
                    case 0:
                        return lhs_m.Size(0);
                    case 1:
                        return rhs_m.Size(1);
                    default:
                        return 0;
                }
            }
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
//#warning matrix operation needs review.
            RET_TYPE ret = 0.0;// = lhs_m(i, 0) * rhs_m(0, j);


            for (size_t k = lhs_m.IndexMin(0); k <= lhs_m.IndexMax(0); k++) {
                ret += lhs_m(i, k) * rhs_m(k, j);
            }

            return ret;
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return std::max(lhs_m.IndexMin(d), rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return std::min(lhs_m.IndexMax(d), rhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            RET_TYPE ret = lhs_m(i, 0) * rhs_m(0, j);


            for (size_t k = 0; k < lhs_m.Size(0); k += 2UL) {
                ret += lhs_m(i, k) * rhs_m(k, j) + lhs_m(i, k + 1UL) * rhs_m(k + 1UL, j);
                //                        ret += lhs_m(i, k + 1UL) * rhs_m(k + 1UL, j);
            }
            if (end < lcols) {
                ret += lhs_m(i, end) * rhs_m(end, j);
            }
            return ret;
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            lhs_m.IsAliased(aliased, ptr);
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }

    };

    template< class LHS, class T >
    struct MatrixMultiplyScalar : MatrixExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, MatrixMultiplyScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit MatrixMultiplyScalar(const MatrixExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) * rhs_m;
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (lhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (lhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) * rhs_m;
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            lhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }



    };

    template<class T, class RHS>
    struct MatrixScalarMultiply : MatrixExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, MatrixScalarMultiply<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixScalarMultiply(const T& lhs, const MatrixExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m * rhs_m(i, j);
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (rhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
              return lhs_m * rhs_m(i, j);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            rhs_m.ExpressionLength(length);
        }


    };
    //
    //typename et4ad::promote_trait<typename LHS::RET_TYPE , typename RHS::RET_TYPE >::return_type

    template <class LHS, class RHS>
    inline const MatrixMultiply< LHS, RHS> operator*(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        //        std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented!!!\n" << std::flush;
        //        exit(0);
        return MatrixMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_MATRIX_MULTIPLY_SCALAR(TYPE) \
    template< class LHS>      \
        inline const MatrixMultiplyScalar<LHS,TYPE> \
    operator*(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return MatrixMultiplyScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_MATRIX_MULTIPLY_SCALAR(short)
    ATL_MATRIX_MULTIPLY_SCALAR(unsigned short)
    ATL_MATRIX_MULTIPLY_SCALAR(int)
    ATL_MATRIX_MULTIPLY_SCALAR(unsigned int)
    ATL_MATRIX_MULTIPLY_SCALAR(long)
    ATL_MATRIX_MULTIPLY_SCALAR(unsigned long)
    ATL_MATRIX_MULTIPLY_SCALAR(float)
    ATL_MATRIX_MULTIPLY_SCALAR(double)
    ATL_MATRIX_MULTIPLY_SCALAR(long double)
    ATL_MATRIX_MULTIPLY_SCALAR(atl::Variable<float>)
    ATL_MATRIX_MULTIPLY_SCALAR(atl::Variable<double>)
    ATL_MATRIX_MULTIPLY_SCALAR(atl::Variable<long double>)


#define ATL_MULTIPLY_SCALAR_MATRIX(TYPE) \
    template< class RHS>      \
        inline const MatrixScalarMultiply<TYPE,RHS> \
    operator*(const TYPE& a, const MatrixExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return MatrixScalarMultiply<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_MULTIPLY_SCALAR_MATRIX(short)
    ATL_MULTIPLY_SCALAR_MATRIX(unsigned short)
    ATL_MULTIPLY_SCALAR_MATRIX(int)
    ATL_MULTIPLY_SCALAR_MATRIX(unsigned int)
    ATL_MULTIPLY_SCALAR_MATRIX(long)
    ATL_MULTIPLY_SCALAR_MATRIX(unsigned long)
    ATL_MULTIPLY_SCALAR_MATRIX(float)
    ATL_MULTIPLY_SCALAR_MATRIX(double)
    ATL_MULTIPLY_SCALAR_MATRIX(long double)
    ATL_MULTIPLY_SCALAR_MATRIX(atl::Variable<float>)
    ATL_MULTIPLY_SCALAR_MATRIX(atl::Variable<double>)
    ATL_MULTIPLY_SCALAR_MATRIX(atl::Variable<long double>)



}




#endif	/* MATRIXMULTIPLY_HPP */

