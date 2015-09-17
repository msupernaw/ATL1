
/* 
 * File:   Matrix.hpp
 * Author: matthewsupernaw
 *
 * Created on October 14, 2014, 11:14 AM
 */

#ifndef MATRIX_HPP
#define	MATRIX_HPP
#include "ContainerDefs.hpp"
#if defined(ATL_CONCURRENCY_ENABLED)
#include  "../third_party/tbb42_20140601oss/include/tbb/concurrent_vector.h"
#else
#include <vector>
#endif

#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "MatrixExpressionBase.hpp"
#include "Array.hpp"
#include "VectorExpressionBase.hpp"
#include "ArrayExpressionBase.hpp"
#include "MatrixArrayOperators.hpp"
#include "VectorExpressionBase.hpp"
#include "VectorArrayOperators.hpp"
#include "MatrixVectorOperators.hpp"
#include "ArrayTraits.hpp"

namespace atl {

    template<class T>
    class Matrix : public MatrixExpression<T, Matrix<T> > {
        //Friends
        //=
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpEquals(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr);
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpEqualsThread(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr, int start, int end);
        //+=
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpPlusEquals(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr);
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpPlusEqualsThread(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr, int start, int end);

        //-=
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpMinuEquals(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr);
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpMinuEqualsThread(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr, int start, int end);
        //*=
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpTimesEquals(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr);
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpTimesEqualsThread(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr, int start, int end);

        // /=
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpDivideEquals(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr);
        template<class TT, class TT2, class RHS>
        friend void ConcurrentOpDivideEqualsThread(Matrix<TT>& v, const MatrixExpression<TT2, RHS>& expr, int start, int end);

    protected:
        size_t isize;
        size_t jsize;
#if defined(ATL_CONCURRENCY_ENABLED)
        tbb::concurrent_vector<T> data_m;
#else
        std::vector<T> data_m;
#endif

        void CheckBounds(const uint32_t& i, const uint32_t& j) const {

            if (i >= isize) {
                std::cout << "-->" << i << " >= " << isize << std::endl;
            }
            assert(i < isize);
            assert(j < jsize);
        }

        inline void Set(size_t i, size_t j, const T& value) {
            data_m[(i * jsize) + j] = value;
        }

        inline T& Get(size_t i, size_t j) {
            return data_m[(i * jsize) + j];
        }

        inline const T& Get(size_t i, size_t j) const {
            return data_m[(i * jsize) + j];
        }

    public:


        typedef typename std::vector<T>::iterator iterator;
        typedef typename std::vector<T>::reverse_iterator reverse_iterator;
        typedef typename std::vector<T>::const_iterator const_iterator;
        typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

        typedef typename IntrinsicBaseType<T>::TYPE INTRINSIC_BASE;
        typedef T RET_TYPE;
        typedef T BASE_TYPE;

        Matrix()
        :
        isize(0), jsize(0) {
        }

        Matrix(const std::initializer_list<std::initializer_list<T> >& l) {
            isize = l.size();
            jsize = l.begin()->size();
            typename std::initializer_list<std::initializer_list<T> >::iterator it;
            typename std::initializer_list<T>::iterator jt;
            data_m.resize(isize * jsize);
            int iindex = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                int jindex = 0;
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    T v = (*jt);
                    //                    data_m[index] = v;
                    Set(iindex, jindex, v);
                    jindex++;

                }
                iindex++;
            }
        }

        Matrix(const size_t& i, const size_t& j)
        :
        isize(i), jsize(j) {
            data_m.resize(i * j);
        }

        Matrix(const Matrix &orig)
        : isize(orig.isize), jsize(orig.jsize) {
            data_m.resize(orig.data_m.size());
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = orig.data_m[i];
            }
            //            data_m.insert(data_m.begin(), orig.data_m.begin(), orig.data_m.end());

        }

        template<class T2, class A>
        Matrix(const MatrixExpression<T2, A> &expr)
        : isize(0), jsize(0) {

            isize = expr.Size(0);
            jsize = expr.Size(1);
            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }
        }

        /**
         * Returns an n x n identity Matrix or a Matrix
         * with value v on the diagonal.
         * 
         * @param rows
         * @param cols
         */
        static const atl::Matrix<T> Identity(size_t n, T v = 1.0) {
            Matrix<T> id(n, n);
            for (int i = 0; i < n; i++) {
                id(i, i) = v;
            }
            return id;
        }

        void Resize(size_t rows, size_t cols) {
            this->isize = rows;
            this->jsize = cols;
            this->data_m.resize(rows * cols);
        }

        void SetBounds(INTRINSIC_BASE minb, INTRINSIC_BASE maxb) {
            std::cout << "warning atl::Matrix<>::" << __func__ << "not implemented for primitive types";
        }

        Matrix& operator=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = val;
            }
            return *this;
        }

        inline Matrix operator-() {
            return static_cast<T>(-1.0) * (*this);
        }

        Matrix& operator=(const std::initializer_list<std::initializer_list<T> >& l) {
            isize = l.size();
            jsize = l.begin()->size();
            typename std::initializer_list<std::initializer_list<T> >::iterator it;
            typename std::initializer_list<T>::iterator jt;
            data_m.resize(isize * jsize);
            int iindex = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                int jindex = 0;
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    T v = (*jt);
                    //                    data_m[index] = v;
                    Set(iindex, jindex, v);
                    jindex++;

                }
                iindex++;
            }
            return *this;
        }

        Matrix& operator=(const Matrix &other) {
            //            data_m.resize(other.data_m.size());
            //            data_m.insert(data_m.begin(), other.data_m.begin(), other.data_m.end());
            isize = other.isize;
            jsize = other.jsize;
            data_m.resize(other.data_m.size());
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = other.data_m[i];
            }
            return *this;
        }

        template<class T2, class A>
        Matrix& operator=(const MatrixExpression< T2, A> &expr) {
            isize = expr.Size(0);
            jsize = expr.Size(1);
            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }

            return *this;

        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator=(const VectorArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() < 2);
#endif

            isize = expr.Size(0);
            jsize = expr.Size(1);

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }
            return *this;
        }
#endif
        template<class T2, class A>
        Matrix& operator=(const ArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() < 2);
#endif

            isize = expr.Size(0);
            jsize = expr.Size(1);

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }
            return *this;
        }

        template<class T2, class A>
        Matrix& operator=(const VectorExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() < 2);
#endif

            isize = expr.Size(0);
            jsize = 1;

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }

            return *this;
        }

        template<class T2, class A>
        Matrix& operator=(const MatrixVectorExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() < 2);
#endif

            isize = expr.Size(0);
            jsize = expr.Size(1);

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator=(const MatrixArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() < 2);
#endif

            isize = expr.Size(0);
            jsize = expr.Size(1);

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    //                    data_m[(i * jsize) + j] = expr(i, j);
                    Set(i, j, expr(i, j));
                }
            }
            return *this;
        }
#endif
        Matrix& operator+=(const T& val) {
            *this = *this+val;
            return *this;
        }

        Matrix& operator+=(const Matrix &other) {
            *this = *this+other;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator+=(const MatrixExpression< T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator+=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }
#endif
        template<class T2, class A>
        Matrix& operator+=(const ArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator+=(const VectorExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator+=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator+=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }
#endif
        Matrix& operator-=(const T& val) {
            *this = *this-val;
            return *this;
        }

        Matrix& operator-=(const Matrix &other) {
            *this = *this-other;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator-=(const MatrixExpression< T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator-=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }
#endif
        template<class T2, class A>
        Matrix& operator-=(const ArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator-=(const VectorExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator-=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator-=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }
#endif
        Matrix& operator*=(const T& val) {
            *this = *this*val;
            return *this;
        }

        Matrix& operator*=(const Matrix &other) {
            *this = *this*other;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator*=(const MatrixExpression< T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator*=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }
#endif
        template<class T2, class A>
        Matrix& operator*=(const ArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator*=(const VectorExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator*=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator*=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }
#endif
        Matrix& operator/=(const T& val) {
            *this = *this / val;
            return *this;
        }

        Matrix& operator/=(const Matrix &other) {
            *this = *this / other;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator/=(const MatrixExpression< T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator/=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }
#endif
        template<class T2, class A>
        Matrix& operator/=(const ArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator/=(const VectorExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Matrix& operator/=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }
#ifdef ATL_CONTAINER_EXPERIMENTAL
        template<class T2, class A>
        Matrix& operator/=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }
#endif
        inline Matrix& operator++() {
            *this = *this+static_cast<T> (1.0);
            return *this;
        }

        inline const Matrix operator++(int i) {
            Matrix temp = *this;
            *this = static_cast<T> (1.0)+ (*this);
            return temp;
        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return isize;
                case 1:
                    return jsize;
                default:
                    return 0;

            }
        }

        inline const T& operator()(size_t i, size_t j) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j);
#endif

            return Get(i, j); //data_m[(i * jsize) + j];
        }

        inline T& operator()(size_t i, size_t j) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j);
#endif

            return Get(i, j); //data_m[(i * jsize) + j];
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return 0;
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {

            switch (d) {
                case 0:
                    return this->isize - 1;
                case 1:
                    return this->jsize - 1;
                default:
                    return 0;
            }
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return data_m[(i * jsize) + j];
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            if ((void*) & this->data_m == ptr) {
                aliased = true;
            }
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
        }

        iterator begin() {
            return this->data_m.begin();
        }

        iterator end() {
            return this->data_m.end();
        }

        const_iterator begin() const {
            return this->data_m.begin();
        }

        const_iterator end()const {
            return this->data_m.end();
        }

        iterator rbegin() {
            return this->data_m.rbegin();
        }

        iterator rend() {
            return this->data_m.rend();
        }

        const_iterator rbegin() const {
            return this->data_m.rbegin();
        }

        const_iterator rend()const {
            return this->data_m.rend();
        }

    };


#define MAKE_VARIABLE_MATRIX_TYPE(TYPE)\
    template<>\
    void Matrix<atl::Variable<TYPE> >::SetBounds(TYPE minb, TYPE maxb) {\
        for (int i = 0; i < this->data_m.size(); i++) {\
            this->data_m[i].SetBounds(minb, maxb);\
        }\
    }\


    MAKE_VARIABLE_MATRIX_TYPE(float)
    MAKE_VARIABLE_MATRIX_TYPE(double)
    MAKE_VARIABLE_MATRIX_TYPE(long double)

    template<class T, class A>
    T Det(const atl::MatrixExpression< T, A>& m) {


        T d = T(1);
        T ratio;
        size_t rows = m.Size(0);
        size_t cols = m.Size(1);
        if (rows != cols) {
            std::cout << "Determinant Error, matrix is not square.";
            return d;
        }



        if (rows == 1 && cols == 1) {
            return m(0, 0);
        }

        if (rows == 2 && cols == 2) {
            d = (m(0, 0) * m(1, 1) -
                    m(0, 1) * m(1, 0));
            return d;
        }


        //BELOW Here IS EXTREMELY SLOW!!
        size_t n = m.Size(0);

        Matrix<T> matrix(rows, cols);
        size_t i, j, k;

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (j > i) {
                    ratio = m(j, i) / m(i, i);
                    for (k = 0; k < n; k++) {
                        T temp = matrix(j, k) - ratio * matrix(i, k);
                        //matrix[j][k] -= ratio * matrix[i][k];


                        //  matrix.SetValue(j, k, temp);
                        matrix(i, k) = temp;
                    }
                }
            }
        }

        for (i = 0; i < n; i++) {
            d *= matrix(i, i); //[i][i];
        }
        return d;

    }

    template<class T, class A>
    T Det(int n, const atl::MatrixExpression< T, A>& mat) {
        T d = 0;
        //        int n = mat.Size(0);
        if (mat.Size(0) != mat.Size(1)) {
            std::cout << __func__ << " Error: matrix must be square";
            exit(0);
        }
        int c, subi, i, j, subj;
        atl::Matrix<T> submat(mat.Size(0), mat.Size(1));
        if (n == 2) {
            return ( (mat(0, 0) * mat(1, 1)) - (mat(1, 0) * mat(0, 1)));
        } else {
            for (c = 0; c < n; c++) {
                subi = 0;
                for (i = 1; i < n; i++) {
                    subj = 0;
                    for (j = 0; j < n; j++) {
                        if (j == c) {
                            continue;
                        }
                        submat(subi, subj) = mat(i, j);
                        subj++;
                    }
                    subi++;
                }
                d = d + (pow(-1.0, (T) c) * mat(0, c) * Det(n - 1, submat));
            }
        }
        return d;
    }

    template<class T2, class A>
    std::ostream& operator<<(std::ostream& out, const atl::MatrixExpression< T2, A> &expr) {

        int dims = expr.Dimensions();
        std::stringstream ss;

        out.precision(4);
        //        out << std::scientific;
        for (int i = 0; i < expr.Size(0); i++) {
            out << "[ ";
            for (int j = 0; j < expr.Size(1); j++) {

                out << std::setw(7) << expr(i, j) << " ";
            }
            out << "]\n";
        }


        return out;
    }

    /**
     * \ingroup Matrix
     * Minor of a Matrix.
     * @param m
     * @param row
     * @param column
     * @return 
     */
    template <class T, class A>
    Matrix<T> Minor(const MatrixExpression<T, A> &m, unsigned int row, unsigned int column) {
        Matrix<T> ret(m.Size(0) - 1, m.Size(1) - 1);
        if (row <= m.Size(0) && column <= m.Size(1)) {
            for (unsigned int i = 0; i < m.Size(0); i++) {
                for (unsigned int j = 0; j < m.Size(1); j++) {

                    if (i != row && j != column) {
                        ret((i - (i > row)), (j - (j > column))) = m(i, j);
                    }

                }
            }
        } else {
            std::cout << "Matrix exception. call to Minor, Index for minor out of range";
            exit(0);
        }
        return ret;
    }

    /**
     * \ingroup Matrix
     * Determine if this Matrix is positive definite.
     * @param m
     * @return 
     */
    template<class T, class A>
    bool IsPositveDefinite(const MatrixExpression<T, A> &m) {

        for (unsigned int i = 0; i < m.Size(0); i++) {

            if (Det(Minor(m, i, i)) <= T(0)) {

                return false;
            }
        }
        return true;
    }

    /**
     * \ingroup Matrix
     * Determine if this Matrix is negative definite.
     * @param m
     * @return 
     */
    template<class T, class A>
    bool IsNegativeDefinite(const MatrixExpression<T, A> &m) {
        for (unsigned int i = 0; i < m.Size(0); i++) {

            if (Det(Minor(m, i, i)) >= T(0)) {

                return false;
            }
        }
        return true;
    }

    /**
     * \ingroup Matrix
     * Determine if this Matrix is semi-positive definite.
     * @param m
     * @return 
     */
    template<class T, class A>
    bool IsPositveSemidefinite(const MatrixExpression<T, A> &m) {
        for (unsigned int i = 0; i < m.Size(0); i++) {

            if (Det(Minor(m, i, i)) < T(0)) {

                return false;
            }
        }
        return true;

    }

    /**
     * \ingroup Matrix
     * Determine if this Matrix is semi-negative definite.
     * @param m
     * @return 
     */
    template<class T, class A>
    bool IsNegativeSemidefinite(const MatrixExpression<T, A> &m) {
        for (unsigned int i = 0; i < m.Size(0); i++) {

            if (Det(Minor(m, i, i)) > T(0)) {

                return false;
            }
        }
        return true;
    }

}
#endif	/* MATRIX_HPP */

