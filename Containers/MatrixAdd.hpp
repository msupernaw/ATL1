/* 
 * File:   MatrixAdd.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 2:28 PM
 */

#ifndef MATRIXADD_HPP
#define	MATRIXADD_HPP

#include "MatrixExpressionBase.hpp"
#include "ArrayTraits.hpp"
#include "../AutoDiff/Variable.hpp"

namespace atl {

    template< class LHS, class RHS>
    struct MatrixAdd : MatrixExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixAdd<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            for (int i = 0; i < 7; i++) {
                assert(lhs_m.Size(i) == rhs_m.Size(i));
            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) + rhs_m(i, j);
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
            return lhs_m.AtRaw(i, j) + rhs_m.AtRaw(i, j);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            lhs_m.IsAliased(aliased, ptr);
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }


    };
    //

    template< class LHS, class T >
    struct MatrixAddScalar : MatrixExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, MatrixAddScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit MatrixAddScalar(const MatrixExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) + rhs_m;
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
            return lhs_m.AtRaw(i, j) + rhs_m;
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
    struct MatrixScalarAdd : MatrixExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, MatrixScalarAdd<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixScalarAdd(const T& lhs, const MatrixExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m + rhs_m(i, j);
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
            return lhs_m + rhs_m.AtRaw(i, j);
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
    inline const MatrixAdd< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixAdd< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_MATRIX_ADD_SCALAR(TYPE) \
    template< class LHS>      \
        inline const MatrixAddScalar<LHS,TYPE> \
    operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return MatrixAddScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_MATRIX_ADD_SCALAR(short)
    ATL_MATRIX_ADD_SCALAR(unsigned short)
    ATL_MATRIX_ADD_SCALAR(int)
    ATL_MATRIX_ADD_SCALAR(unsigned int)
    ATL_MATRIX_ADD_SCALAR(long)
    ATL_MATRIX_ADD_SCALAR(unsigned long)
    ATL_MATRIX_ADD_SCALAR(float)
    ATL_MATRIX_ADD_SCALAR(double)
    ATL_MATRIX_ADD_SCALAR(long double)
    ATL_MATRIX_ADD_SCALAR(atl::Variable<float>)
    ATL_MATRIX_ADD_SCALAR(atl::Variable<double>)
    ATL_MATRIX_ADD_SCALAR(atl::Variable<long double>)


#define ATL_ADD_SCALAR_MATRIX(TYPE) \
    template< class RHS>      \
        inline const MatrixScalarAdd<TYPE,RHS> \
    operator+(const TYPE& a, const MatrixExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return MatrixScalarAdd<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_ADD_SCALAR_MATRIX(short)
    ATL_ADD_SCALAR_MATRIX(unsigned short)
    ATL_ADD_SCALAR_MATRIX(int)
    ATL_ADD_SCALAR_MATRIX(unsigned int)
    ATL_ADD_SCALAR_MATRIX(long)
    ATL_ADD_SCALAR_MATRIX(unsigned long)
    ATL_ADD_SCALAR_MATRIX(float)
    ATL_ADD_SCALAR_MATRIX(double)
    ATL_ADD_SCALAR_MATRIX(long double)
    ATL_ADD_SCALAR_MATRIX(atl::Variable<float>)
    ATL_ADD_SCALAR_MATRIX(atl::Variable<double>)
    ATL_ADD_SCALAR_MATRIX(atl::Variable<long double>)




}




#endif	/* MATRIXADD_HPP */

