/* 
 * File:   VectorExpressionBase.hpp
 * Author: matthewsupernaw
 *
 * Created on October 21, 2014, 10:04 AM
 */

#ifndef VECTOREXPRESSIONBASE_HPP
#define	VECTOREXPRESSIONBASE_HPP
#include "../Traits/Promote.hpp"
#include "../AutoDiff/AutoDiff.hpp"
namespace atl {

    template< class T, class A>
    struct VectorExpression {
        typedef T RET_TYPE;
        typedef T BASE_TYPE;

        const A & Cast() const {
            return static_cast<const A&> (*this);
        }

        inline void ExpressionLength(uint32_t& length)const {
            Cast().ExpressionLength(length);
        }

        inline const size_t Size(const int32_t & dimension) const {
            return Cast().Size(dimension);
        }
        

        inline const size_t Dimensions() const {
            return 1;
        }
        //

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return Cast().operator()(i);
        }

        /*
         *
         *Returns the first valid index for this vector. 
         */
        inline const size_t IndexMin() const {
            return Cast().IndexMin();
        }

        /*
         *
         *Returns the last valid index for this vector. 
         */
        inline const size_t IndexMax() const {
            return Cast().IndexMax();
        }

        /**
         * Get a value based on the raw index for the underlying 
         * data. valid index is 0 - (length -1).
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return Cast().AtRaw(i);
        }
        
        inline void IsAliased(bool& aliased, void* ptr)const{
            Cast().IsAliased(aliased, ptr);
        }

        VectorExpression& operator=(const VectorExpression & exp) const {
            return *this;
        }

    };
    
    
    template<class T2, class A>
    std::ostream& operator<<(std::ostream& out, atl::VectorExpression< T2, A> &expr) {

        int dims = expr.Dimensions();
        std::stringstream ss;

        ss << "[0]:";
        out << ss.str();
        for (int i = 0; i < expr.Size(0); i++) {
            if ((i % 15) == 0 && i != 0) {
                std::cout << "\n" << ss.str();
            }
            out << expr.AtRaw(i) << " ";
        }
        out << "\n";

        return out;
    }

}

#endif	/* VECTOREXPRESSIONBASE_HPP */

