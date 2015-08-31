/* 
 * File:   Vector.hpp
 * Author: matthewsupernaw
 *
 * Created on October 14, 2014, 10:12 AM
 */

#ifndef VECTOR_HPP
#define	VECTOR_HPP

#if defined(ATL_CONCURRENCY_ENABLED)
#include  "tbb42_20140601oss/include/tbb/concurrent_vector.h"
#else
#include <vector>
#endif

#include "Array.hpp"
#include "Matrix.hpp"
#include "VectorExpressionBase.hpp"
#include "ArrayExpressionBase.hpp"
#include "MatrixExpressionBase.hpp"
#include "MatrixArrayOperators.hpp"
#include "VectorExpressionBase.hpp"
#include "VectorArrayOperators.hpp"
#include "MatrixVectorOperators.hpp"
#include "ArrayTraits.hpp"

#include <iostream>
#include <sstream>

namespace atl {

    template<class T>
    class Vector : public VectorExpression<T, Vector<T> > {
    protected:
        size_t isize;

#if defined(ATL_CONCURRENCY_ENABLED)
        tbb::concurrent_vector<T> data_m;
#else
        std::vector<T> data_m;
#endif

        void CheckBounds(const uint32_t& i) const {

            assert(i < isize);

        }


    public:

        typedef T RET_TYPE;
        typedef T BASE_TYPE;

        /**
         * Default constructor.
         * Constructs a 1d Vector.
         * @param i
         */
        Vector(size_t i = 1)
        :
        isize(i) {
            data_m.resize(i);
        }

        Vector(const Vector &orig)
        : isize(orig.isize) {
            this->isize = orig.isize;
            data_m.resize(orig.data_m.size());
#if defined(ATL_CONCURRENCY_ENABLED)
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = orig.data_m[i];
            }
#else
            data_m.insert(data_m.begin(), orig.data_m.begin(), orig.data_m.end());
#endif
        }

        template<class T2, class A>
        Vector(const VectorExpression<T2, A> &expr)
        : isize(0) {

            isize = expr.Size(0);

            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }

        }

        Vector& operator=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = val;
            }
            return *this;
        }

        Vector& operator=(const Vector &other) {
            this->isize = other.isize;
            data_m.resize(other.data_m.size());
#if defined(ATL_CONCURRENCY_ENABLED)
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = other.data_m[i];
            }
#else
            data_m.insert(data_m.begin(), other.data_m.begin(), other.data_m.end());
#endif
            return *this;
        }

        template<class T2, class A>
        Vector& operator=(const VectorExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() == 1);
#endif
            isize = expr.Size(0);


            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
            return *this;
        }

        template<class T2, class A>
        Vector& operator=(const VectorArrayExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() == 1);
#endif
            isize = expr.Size(0);


            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
            return *this;
        }

        template<class T2, class A>
        Vector& operator=(const ArrayExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(expr.Dimensions() == 1);
#endif
            isize = expr.Size(0);


            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
            return *this;
        }

        template<class T2, class A>
        Vector& operator=(const MatrixExpression<T2, A> &expr) {

            isize = expr.Size(0);

#ifdef ENABLE_BOUNDS_CHECKING
            assert(expr.Size(1) == 1);
#endif


            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i, 0);
            }

            return *this;
        }

        template<class T2, class A>
        Vector& operator=(const MatrixVectorExpression<T2, A> &expr) {

            isize = expr.Size(0);

#ifdef ENABLE_BOUNDS_CHECKING
            assert(expr.Size(1) == 1);
#endif


            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i, 0);
            }

            return *this;
        }

        template<class T2, class A>
        Vector& operator=(const MatrixArrayExpression<T2, A> &expr) {

            isize = expr.Size(0);

#ifdef ENABLE_BOUNDS_CHECKING
            assert(expr.Size(1) == 1);
#endif


            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i, 0);
            }

            return *this;
        }

        inline Vector& operator+=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] += val;
            }
            return *this;
        }

        inline Vector& operator+=(const Vector &other) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == other.isize);
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] += other.data_m[i];
            }
            return *this;
        }

        template<class T2, class A>
        inline Vector& operator+=(const VectorExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == expr.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] += expr(i);
            }

            return *this;

        }

        template<class T2, class A>
        Vector& operator+=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator+=(const ArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator+=(const MatrixExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator+=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator+=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        inline Vector& operator-=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] -= val;
            }
            return *this;
        }

        inline Vector& operator-=(const Vector &other) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == other.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] -= other.data_m[i];
            }

            return *this;

        }

        template<class T2, class A>
        inline Vector& operator-=(const VectorExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == expr.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] -= expr(i);
            }

            return *this;

        }

        template<class T2, class A>
        Vector& operator-=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator-=(const ArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator-=(const MatrixExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator-=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator-=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        inline Vector& operator*=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] *= val;
            }
            return *this;
        }

        inline Vector& operator*=(const Vector &other) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == other.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] *= other.data_m[i];
            }

            return *this;

        }

        template<class T2, class A>
        inline Vector& operator*=(const VectorExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == expr.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] *= expr(i);
            }

            return *this;

        }

        template<class T2, class A>
        Vector& operator*=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator*=(const ArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator*=(const MatrixExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator*=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator*=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        inline Vector& operator/=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] /= val;
            }
            return *this;
        }

        inline Vector& operator/=(const Vector &other) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == other.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] /= other.data_m[i];
            }

            return *this;

        }

        template<class T2, class A>
        inline Vector& operator/=(const VectorExpression<T2, A> &expr) {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(isize == expr.Size(0));
#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] /= expr(i);
            }

            return *this;

        }

        template<class T2, class A>
        Vector& operator/=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator/=(const ArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator/=(const MatrixExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator/=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Vector& operator/=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return isize;
                default:
                    return 0;

            }
        }

        inline const size_t Dimensions() const {
            return 1;
        }

        inline const T& operator()(const uint32_t & i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            CheckBounds(i);
#endif

            return data_m[i];
        }

        inline T& operator()(const uint32_t & i) {
#ifdef ENABLE_BOUNDS_CHECKING
            CheckBounds(i);
#endif

            return data_m[i];
        }

        /*
         *
         *Returns the first valid index for this vector. 
         */
        inline const size_t IndexMin() const {
            return 0;
        }

        /*
         *
         *Returns the last valid index for this vector. 
         */
        inline const size_t IndexMax() const {
            return this->isize-1;
        }
        
        inline void IsAliased(bool& aliased, void* ptr) const{
            if((void*)&this->data_m == ptr){
                aliased = true;
            }
        }

        /**
         * Get a value based on the raw index for the underlying 
         * data. valid index is 0 - (length -1).
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return  data_m[i];
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
        }

    };
    
    template<typename REAL_T>
    std::ifstream& operator >> (std::ifstream& in, const atl::Vector<REAL_T>& v){
        for(int i = v.IndexMin(); i<=v.IndexMax(); i++){
            in>>v(i);
        }
        return in;
    }
    

}

#endif	/* VECTOR_HPP */

