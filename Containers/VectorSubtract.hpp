/* 
 * File:   VectorSubtract.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 1:17 PM
 */

#ifndef VECTORSUBTRACT_HPP
#define	VECTORSUBTRACT_HPP
#include "VectorExpressionBase.hpp"
#include "ArrayTraits.hpp"
#include "../AutoDiff/Variable.hpp"

namespace atl {

    template< class LHS, class RHS>
    struct VectorSubtract : VectorExpression<typename atl::PromoteBinaryOpReturnType<typename LHS::RET_TYPE, typename RHS::RET_TYPE, atl::SUBTRACT>::return_type, VectorSubtract<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

//        typedef typename atl::PromoteBinaryOpReturnType<RET_TYPEL, RET_TYPER, atl::SUBTRACT>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;
typedef typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE>::return_type RET_TYPE;
        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorSubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ENABLE_BOUNDS_CHECKING
            for (int i = 0; i < 7; i++) {
                assert(lhs_m.Size(i) == rhs_m.Size(i));
            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        //        inline const size_t Dimensions() const {
        //            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        //        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return lhs_m(i) - rhs_m(i);
        }

        /*
         *
         *Returns the first valid index for this vector. 
         */
        inline const size_t IndexMin() const {
            return std::max(lhs_m.IndexMin(), rhs_m.IndexMin());
        }

        /*
         *
         *Returns the last valid index for this vector. 
         */
        inline const size_t IndexMax() const {
            return std::min(lhs_m.IndexMax(), rhs_m.IndexMax());
        }

        /**
         * Get a value based on the raw index for the underlying 
         * data. valid index is 0 - (length -1).
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return lhs_m(i) - rhs_m(i);
        }

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
    //

    template< class LHS, class T >
    struct VectorSubtractScalar : VectorExpression<typename atl::PromoteBinaryOpReturnType<typename LHS::RET_TYPE, T, atl::SUBTRACT>::return_type, VectorSubtractScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        //        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteBinaryOpReturnType<typename LHS::RET_TYPE, T, atl::SUBTRACT>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit VectorSubtractScalar(const VectorExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        //        inline const size_t Dimensions() const {
        //            return lhs_m.Dimensions();
        //        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return lhs_m(i) - rhs_m;
        }

        /*
         *
         *Returns the first valid index for this vector. 
         */
        inline const size_t IndexMin() const {
            return (lhs_m.IndexMin());
        }

        /*
         *
         *Returns the last valid index for this vector. 
         */
        inline const size_t IndexMax() const {
            return (lhs_m.IndexMax());
        }

        /**
         * Get a value based on the raw index for the underlying 
         * data. valid index is 0 - (length -1).
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return lhs_m(i) - rhs_m;
        }

         inline void IsAliased(bool& aliased, void* ptr) const{
            lhs_m.IsAliased(aliased, ptr);
        }
        
        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }



    };

    template<class T, class RHS>
    struct VectorScalarSubtract : VectorExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, VectorScalarSubtract<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorScalarSubtract(const T& lhs, const VectorExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        //        inline const size_t Dimensions() const {
        //            return rhs_m.Dimensions();
        //        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return lhs_m - rhs_m(i);
        }

        /*
         *
         *Returns the first valid index for this vector. 
         */
        inline const size_t IndexMin() const {
            return (rhs_m.IndexMin());
        }

        /*
         *
         *Returns the last valid index for this vector. 
         */
        inline const size_t IndexMax() const {
            return (rhs_m.IndexMax());
        }

        /**
         * Get a value based on the raw index for the underlying 
         * data. valid index is 0 - (length -1).
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return lhs_m - rhs_m(i);
        }

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
    inline const VectorSubtract< LHS, RHS> operator-(const VectorExpression<typename LHS::RET_TYPE, LHS>& a,
            const VectorExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorSubtract< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_VECTOR_SUBTRACT_SCALAR(TYPE) \
    template< class LHS>      \
        inline const VectorSubtractScalar<LHS,TYPE> \
    operator-(const VectorExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return VectorSubtractScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_VECTOR_SUBTRACT_SCALAR(short)
    ATL_VECTOR_SUBTRACT_SCALAR(unsigned short)
    ATL_VECTOR_SUBTRACT_SCALAR(int)
    ATL_VECTOR_SUBTRACT_SCALAR(unsigned int)
    ATL_VECTOR_SUBTRACT_SCALAR(long)
    ATL_VECTOR_SUBTRACT_SCALAR(unsigned long)
    ATL_VECTOR_SUBTRACT_SCALAR(float)
    ATL_VECTOR_SUBTRACT_SCALAR(double)
    ATL_VECTOR_SUBTRACT_SCALAR(long double)
    ATL_VECTOR_SUBTRACT_SCALAR(atl::Variable<float>)
    ATL_VECTOR_SUBTRACT_SCALAR(atl::Variable<double>)
    ATL_VECTOR_SUBTRACT_SCALAR(atl::Variable<long double>)


#define ATL_SUBTRACT_SCALAR_VECTOR(TYPE) \
    template< class RHS>      \
        inline const VectorScalarSubtract<TYPE,RHS> \
    operator-(const TYPE& a, const VectorExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return VectorScalarSubtract<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_SUBTRACT_SCALAR_VECTOR(short)
    ATL_SUBTRACT_SCALAR_VECTOR(unsigned short)
    ATL_SUBTRACT_SCALAR_VECTOR(int)
    ATL_SUBTRACT_SCALAR_VECTOR(unsigned int)
    ATL_SUBTRACT_SCALAR_VECTOR(long)
    ATL_SUBTRACT_SCALAR_VECTOR(unsigned long)
    ATL_SUBTRACT_SCALAR_VECTOR(float)
    ATL_SUBTRACT_SCALAR_VECTOR(double)
    ATL_SUBTRACT_SCALAR_VECTOR(long double)
    ATL_SUBTRACT_SCALAR_VECTOR(atl::Variable<float>)
    ATL_SUBTRACT_SCALAR_VECTOR(atl::Variable<double>)
    ATL_SUBTRACT_SCALAR_VECTOR(atl::Variable<long double>)


}





#endif	/* VECTORSUBTRACT_HPP */

