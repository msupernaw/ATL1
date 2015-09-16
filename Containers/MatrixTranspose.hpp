/* 
 * File:   MatrixTranspose.hpp
 * Author: matthewsupernaw
 *
 * Created on November 10, 2014, 11:29 AM
 */

#ifndef MATRIXTRANSPOSE_HPP
#define	MATRIXTRANSPOSE_HPP

#include "MatrixExpressionBase.hpp"
namespace atl {

    template<class M>
    struct MatrixTranspose : atl::MatrixExpression<typename M::RET_TYPE, MatrixTranspose<M> > {
        const M& m_m;
        typedef typename M::RET_TYPE RET_TYPE;
        typedef typename M::BASE_TYPE BASE_TYPE;

        inline explicit MatrixTranspose(const M& m) : m_m(m) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return m_m.Size(1);
                case 1:
                    return m_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(i < m_m.Size(1));
            assert(j < m_m.Size(0));
#endif   

            return m_m(j, i);

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
        inline void IsAliased(bool& aliased, void* ptr) const{
            lhs_m.IsAliased(aliased, ptr);
            rhs_m.IsAliased(aliased, ptr);
        }


        inline void ExpressionLength(uint32_t& length)const {
            length++;
            m_m.ExpressionLength(length);
        }

    };

    template<class M>
    inline const MatrixTranspose<M> Transpose(const M& m) {
        return MatrixTranspose<M>(m);
    }
}

#endif	/* MATRIXTRANSPOSE_HPP */

