/* 
 * File:   MatrixColumn.hpp
 * Author: matthewsupernaw
 *
 * Created on December 15, 2014, 9:00 AM
 */

#ifndef MATRIXCOLUMN_HPP
#define	MATRIXCOLUMN_HPP

namespace atl {

    template<class MAT>
    struct MatrixColumn :
    VectorExpression<MAT::RET_TYPE, MatrixColumn<MAT> > {
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

#endif	/* MATRIXCOLUMN_HPP */

