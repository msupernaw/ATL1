/* 
 * File:   MatrixRow.hpp
 * Author: matthewsupernaw
 *
 * Created on December 13, 2014, 8:11 AM
 */

#ifndef MATRIXROW_HPP
#define	MATRIXROW_HPP

#include "VectorExpressionBase.hpp"
#include "MatrixExpressionBase.hpp"

namespace atl {

    template<class MAT>
    struct MatrixRow :
    VectorExpression<MAT::RET_TYPE, MatrixRow<MAT> > {
        MAT& matrix;
        uint32_t row;
        typedef typename MAT::BASE_TYPE BASE_TYPE;

        MatrixRow(const MAT& mat, uint32_t row) : matrix(mat.Cast()), row(row) {

        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            matrix.ExpressionLength(length);
        }

    };



}


#endif	/* MATRIXROW_HPP */

