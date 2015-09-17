/* 
 * File:   Array.hpp
 * Author: matthewsupernaw
 *
 * Created on August 25, 2014, 9:22 AM
 */

#ifndef ARRAY_HPP
#define	ARRAY_HPP

#include <assert.h>
#include <iostream>
#include <sstream>
#include <initializer_list>

#if defined(ATL_CONCURRENCY_ENABLED)
#include  "tbb42_20140601oss/include/tbb/concurrent_vector.h"
#else
#include <vector>
#endif


#include "ArrayExpressionBase.hpp"
#include "MatrixExpressionBase.hpp"
#include "MatrixArrayOperators.hpp"
#include "VectorExpressionBase.hpp"
#include "VectorArrayOperators.hpp"
#include "MatrixVectorOperators.hpp"
#include "ArrayTraits.hpp"


//#define ATL_ENABLE_BOUNDS_CHECKING


namespace atl {

    template<class T, int group = 0 >
    class Array : public ArrayExpression< T, Array<T, group> > {
    protected:
        size_t dims_m;
        size_t isize;
        size_t jsize;
        size_t ksize;
        size_t lsize;
        size_t msize;
        size_t nsize;
        size_t osize;

#if defined(ATL_CONCURRENCY_ENABLED)
        tbb::concurrent_vector<T> data_m;
#else
        std::vector<T> data_m;
#endif

        void CheckBounds(const uint32_t& i,
                const uint32_t& j = 0,
                const uint32_t& k = 0,
                const uint32_t& l = 0,
                const uint32_t& m = 0,
                const uint32_t& n = 0,
                const uint32_t& o = 0) const {

            assert(i < isize);
            assert(j < jsize);
            assert(k < ksize);
            assert(l < lsize);
            assert(m < msize);
            assert(n < nsize);
            assert(o < osize);
        }

        //        template<class A>
        //        inline void SetData(const ArrayExpression<A> &expr) const {
        //
        //
        //
        //        }

    public:

        typedef typename std::vector<T>::iterator iterator;
        typedef typename std::vector<T>::reverse_iterator reverse_iterator;
        typedef typename std::vector<T>::const_iterator const_iterator;
        typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

        typedef typename IntrinsicBaseType<T>::TYPE INTRINSIC_BASE;
        typedef T RET_TYPE;
        typedef T BASE_TYPE;

        /**
         * Default constructor.
         * Constructs a 1d Array.
         * @param i
         */
        Array(size_t i = 1)
        : dims_m(1),
        isize(i),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            data_m.resize(i, T(0));
        }

        /**
         *  Constructs a 2d Array.
         * @param i
         * @param j
         */
        Array(size_t i, size_t j)
        : dims_m(2),
        isize(i),
        jsize(j),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {

            data_m.resize(i * j);

        }

        /**
         *  Constructs a 3d Array.
         * @param i
         * @param j
         * @param k
         */
        Array(size_t i, size_t j, size_t k)
        : dims_m(3),
        isize(i),
        jsize(j),
        ksize(k),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            data_m.resize((i * j * k), T(0));

        }

        /**
         *  Constructs a 4d Array.
         * @param i
         * @param j
         * @param k
         * @param l
         */
        Array(size_t i, size_t j, size_t k, size_t l)
        : dims_m(4),
        isize(i),
        jsize(j),
        ksize(k),
        lsize(l),
        msize(1),
        nsize(1),
        osize(1) {
            data_m.resize((i * j * k * l), T(0));

        }

        /**
         *  Constructs a 5d Array.
         * @param i
         * @param j
         * @param k
         * @param l
         * @param m
         */
        Array(size_t i, size_t j, size_t k, size_t l, size_t m)
        : dims_m(5),
        isize(i),
        jsize(j),
        ksize(k),
        lsize(l),
        msize(m),
        nsize(1),
        osize(1) {
            data_m.resize((i * j * k * l * m), T(0));

        }

        /**
         *  Constructs a 6d Array.
         * @param i
         * @param j
         * @param k
         * @param l
         * @param m
         * @param n
         * @param o
         */
        Array(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n)
        : dims_m(6),
        isize(i),
        jsize(j),
        ksize(k),
        lsize(l),
        msize(m),
        nsize(n),
        osize(1) {
            data_m.resize((i * j * k * l * m * n), T(0));

        }

        /**
         *  Constructs a 7d Array.
         * @param i
         * @param j
         * @param k
         * @param l
         * @param m
         * @param n
         * @param o
         */
        Array(size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o)
        : dims_m(7),
        isize(i),
        jsize(j),
        ksize(k),
        lsize(l),
        msize(m),
        nsize(n),
        osize(o) {
            data_m.resize((i * j * k * l * m * n * o), T(0));

        }

        Array(const Array &orig)
        : dims_m(orig.dims_m),
        isize(orig.isize),
        jsize(orig.jsize),
        ksize(orig.ksize),
        lsize(orig.lsize),
        msize(orig.msize),
        nsize(orig.nsize),
        osize(orig.osize) {
            data_m.resize(orig.data_m.size());
            data_m.insert(data_m.begin(), orig.data_m.begin(), orig.data_m.end());

        }

        Array(const std::initializer_list<T>& l)
        : dims_m(1),
        isize(l.size()),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            typename std::initializer_list<T>::iterator it;
            data_m.resize(isize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                T v = (*it);
                data_m[index] = v;
                index++;
            }
        }

        Array(const std::initializer_list<std::initializer_list<T> >& l)
        : dims_m(2),
        isize(l.size()),
        jsize(l.begin()->size()),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            typename std::initializer_list<std::initializer_list<T> >::iterator it;
            typename std::initializer_list<T>::iterator jt;
            data_m.resize(isize * jsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    T v = (*jt);
                    data_m[index] = v;
                    index++;
                }
            }
        }

        Array(const std::initializer_list<std::initializer_list<std::initializer_list<T> > >& l)
        : dims_m(3),
        isize(l.size()),
        jsize(l.begin()->size()),
        ksize(l.begin()->begin()->size()),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator it;
            typename std::initializer_list<std::initializer_list<T> >::iterator jt;
            typename std::initializer_list<T>::iterator kt;
            data_m.resize(isize * jsize * ksize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        T v = (*kt);
                        data_m[index] = v;
                        index++;
                    }
                }
            }
        }

        Array(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >& l)
        : dims_m(4),
        isize(l.size()),
        jsize(l.begin()->size()),
        ksize(l.begin()->begin()->size()),
        lsize(l.begin()->begin()->begin()->size()),
        msize(1),
        nsize(1),
        osize(1) {
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator jt;
            typename std::initializer_list<std::initializer_list<T> >::iterator kt;
            typename std::initializer_list<T>::iterator lt;
            data_m.resize(isize * jsize * ksize * lsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            T v = (*lt);
                            data_m[index] = v;
                            index++;
                        }
                    }
                }
            }
        }

        Array(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >& l)
        : dims_m(5),
        isize(l.size()),
        jsize(l.begin()->size()),
        ksize(l.begin()->begin()->size()),
        lsize(l.begin()->begin()->begin()->size()),
        msize(l.begin()->begin()->begin()->begin()->size()),
        nsize(1),
        osize(1) {
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator jt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator kt;
            typename std::initializer_list<std::initializer_list<T> >::iterator lt;
            typename std::initializer_list<T>::iterator mt;
            ;
            data_m.resize(isize * jsize * ksize * lsize * msize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            for (mt = lt->begin(); mt != lt->end(); ++mt) {
                                T v = (*mt);
                                data_m[index] = v;
                                index++;
                            }
                        }
                    }
                }
            }
        }

        Array(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > >& l)
        : dims_m(6),
        isize(l.size()),
        jsize(l.begin()->size()),
        ksize(l.begin()->begin()->size()),
        lsize(l.begin()->begin()->begin()->size()),
        msize(l.begin()->begin()->begin()->begin()->size()),
        nsize(l.begin()->begin()->begin()->begin()->begin()->size()),
        osize(1) {
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >::iterator jt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator kt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator lt;
            typename std::initializer_list<std::initializer_list<T> >::iterator mt;
            typename std::initializer_list<T>::iterator nt;

            data_m.resize(isize * jsize * ksize * lsize * msize * nsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            for (mt = lt->begin(); mt != lt->end(); ++mt) {
                                for (nt = mt->begin(); nt != mt->end(); ++nt) {
                                    T v = (*nt);
                                    data_m[index] = v;
                                    index++;
                                }
                            }
                        }
                    }
                }
            }
        }

        Array(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > > >& l)
        : dims_m(7),
        isize(l.size()),
        jsize(l.begin()->size()),
        ksize(l.begin()->begin()->size()),
        lsize(l.begin()->begin()->begin()->size()),
        msize(l.begin()->begin()->begin()->begin()->size()),
        nsize(l.begin()->begin()->begin()->begin()->begin()->size()),
        osize(l.begin()->begin()->begin()->begin()->begin()->begin()->size()) {
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > >::iterator jt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >::iterator kt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator lt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator mt;
            typename std::initializer_list<std::initializer_list<T> >::iterator nt;
            typename std::initializer_list<T>::iterator ot;
            ;
            data_m.resize(isize * jsize * ksize * lsize * msize * nsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            for (mt = lt->begin(); mt != lt->end(); ++mt) {
                                for (nt = mt->begin(); nt != mt->end(); ++nt) {
                                    for (ot = nt->begin(); ot != nt->end(); ++ot) {
                                        T v = (*ot);
                                        data_m[index] = v;
                                        index++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        template<class T2, class A>
        Array(const ArrayExpression<T2, A> &expr)
        : isize(0),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {

            dims_m = expr.Dimensions();

            for (int i = 0; i < dims_m; i++) {
                switch (dims_m) {
                    case 1:
                        isize = expr.Size(i);
                        break;
                    case 2:
                        jsize = expr.Size(i);
                        break;
                    case 3:
                        ksize = expr.Size(i);
                        break;
                    case 4:
                        ksize = expr.Size(i);
                        break;
                    case 5:
                        lsize = expr.Size(i);
                        break;
                    case 6:
                        msize = expr.Size(i);
                        break;
                    case 7:
                        nsize = expr.Size(i);
                        break;
                    case 8:
                        osize = expr.Size(i);
                        break;
                }
            }


            //            this->SetData(expr);
            switch (dims_m) {

                case 1:
                    data_m.resize(isize);

                    for (int i = 0; i < isize; i++) {
                        data_m[i] = expr(i);
                    }
                    break;
                case 2:
                    data_m.resize(isize * jsize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            data_m[(i * jsize) + j] = expr(i, j);
                        }
                    }
                    break;
                case 3:
                    data_m.resize(isize * jsize * ksize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                data_m[(i * jsize) + (j * ksize) + k] = expr(i, j, k);
                            }
                        }
                    }
                    break;
                case 4:
                    data_m.resize(isize * jsize * ksize * lsize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] = expr(i, j, k, l);
                                }
                            }
                        }
                    }
                    break;

                case 5:
                    data_m.resize(isize * jsize * ksize * lsize * msize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] = expr(i, j, k, l, m);
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 6:
                    data_m.resize(isize * jsize * ksize * lsize * msize * nsize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] = expr(i, j, k, l, m, n);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case 7:
                    data_m.resize(isize * jsize * ksize * lsize * msize * nsize, osize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            for (int o = 0; o < osize; o++) {
                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] = expr(i, j, k, l, m, n, o);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;

            }


        }

        template<class T2, class A>
        Array(const VectorExpression<T2, A> &expr)
        : isize(0),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
        }

        template<class T2, class A>
        Array(const VectorArrayExpression<T2, A> &expr)
        : isize(0),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
        }

        template<class T2, class A>
        Array(const MatrixExpression<T2, A> &expr)
        : isize(0),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = expr.Size(1);
            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    data_m[(i * jsize) + j] = expr(i, j);
                }
            }
        }

        template<class T2, class A>
        Array(const MatrixVectorExpression<T2, A> &expr)
        : isize(0),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = expr.Size(1);

            if (jsize == 0) {
                jsize = 1;
            }
            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    data_m[(i * jsize) + j] = expr(i, j);
                }
            }
        }

        template<class T2, class A>
        Array(const MatrixArrayExpression<T2, A> &expr)
        : isize(0),
        jsize(1),
        ksize(1),
        lsize(1),
        msize(1),
        nsize(1),
        osize(1) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = expr.Size(1);
            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    data_m[(i * jsize) + j] = expr(i, j);
                }
            }
        }

        void Resize(size_t i, size_t j =1, size_t k=1, size_t l=1, size_t m=1, size_t n=1, size_t o=1){
            data_m.resize((i * j * k * l * m * n * o), T(0));
        }
        
        
        Array& operator=(const T& val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] = val;
            }
            return *this;
        }

        inline Array operator-() {
            return static_cast<T> (-1.0) * (*this);
        }

        Array& operator=(const Array &other) {
            isize = other.isize;
            jsize = other.jsize;
            ksize = other.ksize;
            lsize = other.lsize;
            msize = other.msize;
            nsize = other.nsize;
            osize = other.osize;
            dims_m = other.Dimensions();
            data_m = other.data_m;
            return *this;
        }

        template<class T2, class A>
        Array& operator=(const ArrayExpression<T2, A> &expr) {

            isize = 1;
            jsize = 1;
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;
            dims_m = expr.Dimensions();

            for (int i = 0; i < dims_m; i++) {
                switch (dims_m) {
                    case 1:
                        isize = expr.Size(i);
                        break;
                    case 2:
                        jsize = expr.Size(i);
                        break;
                    case 3:
                        ksize = expr.Size(i);
                        break;
                    case 4:
                        ksize = expr.Size(i);
                        break;
                    case 5:
                        lsize = expr.Size(i);
                        break;
                    case 6:
                        msize = expr.Size(i);
                        break;
                    case 7:
                        nsize = expr.Size(i);
                        break;
                    case 8:
                        osize = expr.Size(i);
                        break;
                }
            }


            //            this->SetData(expr);
            switch (dims_m) {

                case 1:
                    data_m.resize(isize);

                    for (int i = 0; i < isize; i++) {
                        data_m[i] = expr(i);
                    }
                    break;
                case 2:
                    data_m.resize(isize * jsize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            data_m[(i * jsize) + j] = expr(i, j);
                        }
                    }
                    break;
                case 3:
                    data_m.resize(isize * jsize * ksize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                data_m[(i * jsize) + (j * ksize) + k] = expr(i, j, k);
                            }
                        }
                    }
                    break;
                case 4:
                    data_m.resize(isize * jsize * ksize * lsize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] = expr(i, j, k, l);
                                }
                            }
                        }
                    }
                    break;

                case 5:
                    data_m.resize(isize * jsize * ksize * lsize * msize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] = expr(i, j, k, l, m);
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 6:
                    data_m.resize(isize * jsize * ksize * lsize * msize * nsize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] = expr(i, j, k, l, m, n);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case 7:
                    data_m.resize(isize * jsize * ksize * lsize * msize * nsize, osize);
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            for (int o = 0; o < osize; o++) {
                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] = expr(i, j, k, l, m, n, o);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;

            }

            return *this;

        }

        Array& operator=(const std::initializer_list<T>& l) {
            dims_m = 1;
            isize = l.size();
            jsize = 1;
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;
            typename std::initializer_list<T>::iterator it;
            data_m.resize(isize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                T v = (*it);
                data_m[index] = v;
                index++;
            }
        }

        Array& operator=(const std::initializer_list<std::initializer_list<T> >& l) {
            dims_m = 2;
            isize = l.size();
            jsize = l.begin()->size();
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;
            typename std::initializer_list<std::initializer_list<T> >::iterator it;
            typename std::initializer_list<T>::iterator jt;
            data_m.resize(isize * jsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    T v = (*jt);
                    data_m[index] = v;
                    index++;
                }
            }
        }

        Array& operator=(const std::initializer_list<std::initializer_list<std::initializer_list<T> > >& l) {
            dims_m = 3;
            isize = l.size();
            jsize = l.begin()->size();
            ksize = l.begin()->begin()->size();
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator it;
            typename std::initializer_list<std::initializer_list<T> >::iterator jt;
            typename std::initializer_list<T>::iterator kt;
            data_m.resize(isize * jsize * ksize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        T v = (*kt);
                        data_m[index] = v;
                        index++;
                    }
                }
            }
        }

        Array& operator=(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >& l) {
            dims_m = 4;
            isize = l.size();
            jsize = l.begin()->size();
            ksize = l.begin()->begin()->size();
            lsize = l.begin()->begin()->begin()->size();
            msize = 1;
            nsize = 1;
            osize = 1;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator jt;
            typename std::initializer_list<std::initializer_list<T> >::iterator kt;
            typename std::initializer_list<T>::iterator lt;
            data_m.resize(isize * jsize * ksize * lsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            T v = (*lt);
                            data_m[index] = v;
                            index++;
                        }
                    }
                }
            }
        }

        Array& operator=(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >& l) {
            dims_m = 5;
            isize = l.size();
            jsize = l.begin()->size();
            ksize = l.begin()->begin()->size();
            lsize = l.begin()->begin()->begin()->size();
            msize = l.begin()->begin()->begin()->begin()->size();
            nsize = 1;
            osize = 1;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator jt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator kt;
            typename std::initializer_list<std::initializer_list<T> >::iterator lt;
            typename std::initializer_list<T>::iterator mt;
            ;
            data_m.resize(isize * jsize * ksize * lsize * msize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            for (mt = lt->begin(); mt != lt->end(); ++mt) {
                                T v = (*mt);
                                data_m[index] = v;
                                index++;
                            }
                        }
                    }
                }
            }
        }

        Array& operator=(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > >& l) {
            dims_m = 6;
            isize = l.size();
            jsize = l.begin()->size();
            ksize = l.begin()->begin()->size();
            lsize = l.begin()->begin()->begin()->size();
            msize = l.begin()->begin()->begin()->begin()->size();
            nsize = l.begin()->begin()->begin()->begin()->begin()->size();
            osize = 1;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >::iterator jt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator kt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator lt;
            typename std::initializer_list<std::initializer_list<T> >::iterator mt;
            typename std::initializer_list<T>::iterator nt;

            data_m.resize(isize * jsize * ksize * lsize * msize * nsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            for (mt = lt->begin(); mt != lt->end(); ++mt) {
                                for (nt = mt->begin(); nt != mt->end(); ++nt) {
                                    T v = (*nt);
                                    data_m[index] = v;
                                    index++;
                                }
                            }
                        }
                    }
                }
            }
        }

        Array& operator=(const std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > > >& l) {
            dims_m = 7;
            isize = l.size();
            jsize = l.begin()->size();
            ksize = l.begin()->begin()->size();
            lsize = l.begin()->begin()->begin()->size();
            msize = l.begin()->begin()->begin()->begin()->size();
            nsize = l.begin()->begin()->begin()->begin()->begin()->size();
            osize = l.begin()->begin()->begin()->begin()->begin()->begin()->size();
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > > >::iterator it;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > > >::iterator jt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > > >::iterator kt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<std::initializer_list<T> > > >::iterator lt;
            typename std::initializer_list<std::initializer_list<std::initializer_list<T> > >::iterator mt;
            typename std::initializer_list<std::initializer_list<T> >::iterator nt;
            typename std::initializer_list<T>::iterator ot;
            ;
            data_m.resize(isize * jsize * ksize * lsize * msize * nsize);
            int index = 0;
            for (it = l.begin(); it != l.end(); ++it) {
                for (jt = it->begin(); jt != it->end(); ++jt) {
                    for (kt = jt->begin(); kt != jt->end(); ++kt) {
                        for (lt = kt->begin(); lt != kt->end(); ++lt) {
                            for (mt = lt->begin(); mt != lt->end(); ++mt) {
                                for (nt = mt->begin(); nt != mt->end(); ++nt) {
                                    for (ot = nt->begin(); ot != nt->end(); ++ot) {
                                        T v = (*ot);
                                        data_m[index] = v;
                                        index++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        template<class T2, class A>
        Array& operator=(const VectorExpression<T2, A> &expr) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = 1;
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;

            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
            return *this;
        }

        template<class T2, class A>
        Array& operator=(const VectorArrayExpression<T2, A> &expr) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = 1;
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;

            data_m.resize(isize);

            for (int i = 0; i < isize; i++) {
                data_m[i] = expr(i);
            }
            return *this;
        }

        template<class T2, class A>
        Array& operator=(const MatrixExpression<T2, A> &expr) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = expr.Size(1);
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    data_m[(i * jsize) + j] = expr(i, j);
                }
            }

            return *this;
        }

        template<class T2, class A>
        Array& operator=(const MatrixVectorExpression<T2, A> &expr) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = expr.Size(1);
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;

            if (jsize == 0) {
                jsize = 1;
            }
            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    data_m[(i * jsize) + j] = expr(i, j);
                }
            }
            return *this;
        }

        template<class T2, class A>
        Array& operator=(const MatrixArrayExpression<T2, A> &expr) {
            dims_m = expr.Dimensions();
            isize = expr.Size(0);
            jsize = expr.Size(1);
            ksize = 1;
            lsize = 1;
            msize = 1;
            nsize = 1;
            osize = 1;

            data_m.resize(isize * jsize);

            for (int i = 0; i < isize; i++) {
                for (int j = 0; j < jsize; j++) {
                    data_m[(i * jsize) + j] = expr(i, j);
                }
            }
            return *this;
        }

        inline Array& operator+=(const T & val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] += val;
            }
            return *this;
        }

        inline Array& operator+=(const Array & other) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == other.Dimensions());
            assert(isize == other.isize);
            assert(jsize == other.jsize);
            assert(ksize == other.ksize);
            assert(lsize == other.lsize);
            assert(msize == other.msize);
            assert(nsize == other.nsize);
            assert(osize == other.osize);

#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] += other.data_m[i];
            }
            return *this;
        }

        template< class T2, class A>
        inline Array& operator+=(const ArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == expr.Dimensions());
            assert(isize == expr.Size(0));
            assert(jsize == expr.Size(1));
            assert(ksize == expr.Size(2));
            assert(lsize == expr.Size(3));
            assert(msize == expr.Size(4));
            assert(nsize == expr.Size(5));
            assert(osize == expr.Size(6));

#endif  


            switch (dims_m) {

                case 1:
                    for (int i = 0; i < isize; i++) {
                        data_m[i] += expr(i);
                    }
                    break;
                case 2:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            data_m[(i * jsize) + j] += expr(i, j);
                        }
                    }
                    break;
                case 3:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                data_m[(i * jsize) + (j * ksize) + k] += expr(i, j, k);
                            }
                        }
                    }
                    break;
                case 4:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] += expr(i, j, k, l);
                                }
                            }
                        }
                    }
                    break;

                case 5:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] += expr(i, j, k, l, m);
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 6:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] += expr(i, j, k, l, m, n);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case 7:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            for (int o = 0; o < osize; o++) {
                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] += expr(i, j, k, l, m, n, o);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;

            }

            return *this;

        }

        template<class T2, class A>
        Array& operator+=(const VectorExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator+=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator+=(const MatrixExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator+=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this+Array<T2>(expr);
            return *this;
        }

        template<class T2, class A>
        Array& operator+=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this+expr;
            return *this;
        }

        inline Array& operator-=(const T & val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] -= val;
            }
        }

        inline Array& operator-=(const Array & other) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == other.Dimensions());
            assert(isize == other.isize);
            assert(jsize == other.jsize);
            assert(ksize == other.ksize);
            assert(lsize == other.lsize);
            assert(msize == other.msize);
            assert(nsize == other.nsize);
            assert(osize == other.osize);

#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] -= other.data_m[i];
            }
            return *this;
        }

        template<class T2, class A>
        inline Array& operator-=(const ArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == expr.Dimensions());
            assert(isize == expr.Size(0));
            assert(jsize == expr.Size(1));
            assert(ksize == expr.Size(2));
            assert(lsize == expr.Size(3));
            assert(msize == expr.Size(4));
            assert(nsize == expr.Size(5));
            assert(osize == expr.Size(6));

#endif  


            switch (dims_m) {

                case 1:
                    for (int i = 0; i < isize; i++) {
                        data_m[i] -= expr(i);
                    }
                    break;
                case 2:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            data_m[(i * jsize) + j] -= expr(i, j);
                        }
                    }
                    break;
                case 3:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                data_m[(i * jsize) + (j * ksize) + k] -= expr(i, j, k);
                            }
                        }
                    }
                    break;
                case 4:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] -= expr(i, j, k, l);
                                }
                            }
                        }
                    }
                    break;

                case 5:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] -= expr(i, j, k, l, m);
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 6:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] -= expr(i, j, k, l, m, n);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case 7:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            for (int o = 0; o < osize; o++) {
                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] -= expr(i, j, k, l, m, n, o);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;

            }

            return *this;

        }

        template<class T2, class A>
        Array& operator-=(const VectorExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator-=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator-=(const MatrixExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator-=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this-Array<T2>(expr);
            return *this;
        }

        template<class T2, class A>
        Array& operator-=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this-expr;
            return *this;
        }

        inline Array& operator*=(const T & val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] *= val;
            }
            return *this;
        }

        inline Array& operator*=(const Array & other) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == other.Dimensions());
            assert(isize == other.isize);
            assert(jsize == other.jsize);
            assert(ksize == other.ksize);
            assert(lsize == other.lsize);
            assert(msize == other.msize);
            assert(nsize == other.nsize);
            assert(osize == other.osize);

#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] *= other.data_m[i];
            }
            return *this;
        }

        template<class T2, class A>
        inline Array& operator*=(const ArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == expr.Dimensions());
            assert(isize == expr.Size(0));
            assert(jsize == expr.Size(1));
            assert(ksize == expr.Size(2));
            assert(lsize == expr.Size(3));
            assert(msize == expr.Size(4));
            assert(nsize == expr.Size(5));
            assert(osize == expr.Size(6));

#endif  


            switch (dims_m) {

                case 1:
                    for (int i = 0; i < isize; i++) {
                        data_m[i] *= expr(i);
                    }
                    break;
                case 2:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            data_m[(i * jsize) + j] *= expr(i, j);
                        }
                    }
                    break;
                case 3:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                data_m[(i * jsize) + (j * ksize) + k] *= expr(i, j, k);
                            }
                        }
                    }
                    break;
                case 4:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] *= expr(i, j, k, l);
                                }
                            }
                        }
                    }
                    break;

                case 5:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] *= expr(i, j, k, l, m);
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 6:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] *= expr(i, j, k, l, m, n);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case 7:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            for (int o = 0; o < osize; o++) {
                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] *= expr(i, j, k, l, m, n, o);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;

            }

            return *this;

        }

        template<class T2, class A>
        Array& operator*=(const VectorExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator*=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator*=(const MatrixExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator*=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this*Array<T2>(expr);
            ;
            return *this;
        }

        template<class T2, class A>
        Array& operator*=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this*expr;
            return *this;
        }

        inline Array& operator/=(const T & val) {
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] /= val;
            }
        }

        inline Array& operator/=(const Array & other) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == other.Dimensions());
            assert(isize == other.isize);
            assert(jsize == other.jsize);
            assert(ksize == other.ksize);
            assert(lsize == other.lsize);
            assert(msize == other.msize);
            assert(nsize == other.nsize);
            assert(osize == other.osize);

#endif   
            for (int i = 0; i < data_m.size(); i++) {
                data_m[i] /= other.data_m[i];
            }
            return *this;
        }

        template<class T2, class A>
        inline Array& operator/=(const ArrayExpression<T2, A> &expr) {

#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(dims_m == expr.Dimensions());
            assert(isize == expr.Size(0));
            assert(jsize == expr.Size(1));
            assert(ksize == expr.Size(2));
            assert(lsize == expr.Size(3));
            assert(msize == expr.Size(4));
            assert(nsize == expr.Size(5));
            assert(osize == expr.Size(6));

#endif  


            switch (dims_m) {

                case 1:
                    for (int i = 0; i < isize; i++) {
                        data_m[i] /= expr(i);
                    }
                    break;
                case 2:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            data_m[(i * jsize) + j] /= expr(i, j);
                        }
                    }
                    break;
                case 3:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                data_m[(i * jsize) + (j * ksize) + k] /= expr(i, j, k);
                            }
                        }
                    }
                    break;
                case 4:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] /= expr(i, j, k, l);
                                }
                            }
                        }
                    }
                    break;

                case 5:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] /= expr(i, j, k, l, m);
                                    }
                                }
                            }
                        }
                    }
                    break;
                case 6:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] /= expr(i, j, k, l, m, n);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                case 7:
                    for (int i = 0; i < isize; i++) {
                        for (int j = 0; j < jsize; j++) {
                            for (int k = 0; k < ksize; k++) {
                                for (int l = 0; l < lsize; l++) {
                                    for (int m = 0; m < msize; m++) {
                                        for (int n = 0; n < nsize; n++) {
                                            for (int o = 0; o < osize; o++) {
                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] /= expr(i, j, k, l, m, n, o);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;

            }

            return *this;

        }

        template<class T2, class A>
        Array& operator/=(const VectorExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator/=(const VectorArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator/=(const MatrixExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        template<class T2, class A>
        Array& operator/=(const MatrixVectorExpression<T2, A> &expr) {
            *this = *this / Array<T2>(expr);
            return *this;
        }

        template<class T2, class A>
        Array& operator/=(const MatrixArrayExpression<T2, A> &expr) {
            *this = *this / expr;
            return *this;
        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return isize;
                case 1:
                    return jsize;
                case 2:
                    return ksize;
                case 3:
                    return lsize;
                case 4:
                    return msize;
                case 5:
                    return nsize;
                case 6:
                    return osize;
                default:
                    return 0;

            }
        }

        inline const size_t Dimensions() const {
            return dims_m;
        }

        inline const T& operator()(const uint32_t & i) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i);
#endif

            return data_m[i];
        }

        inline const T& operator()(const uint32_t& i, const uint32_t & j) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j);
#endif

            return data_m[i * jsize + j];
        }

        inline const T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k);
#endif

            return data_m[(i * jsize * ksize)+(j * ksize) + k];
        }

        inline const T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l);
#endif

            return data_m[(i * jsize * ksize * lsize)+(j * ksize * lsize) + (k * lsize) + l];
        }

        inline const T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l, m);
#endif

            return data_m[(i * jsize * ksize * lsize * msize)+(j * ksize * lsize * msize) + (k * lsize * msize) + (l * msize) + m];
        }

        inline const T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l, m, n);
#endif

            return data_m[(i * jsize * ksize * lsize * msize * nsize)+(j * ksize * lsize * msize * nsize)
                    + (k * lsize * msize * nsize) + (l * msize * nsize) + (m * nsize) + n];
        }

        inline const T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l, m, n, o);
#endif

            return data_m[(i * jsize * ksize * lsize * msize * nsize * osize)+
                    (j * ksize * lsize * msize * nsize * osize) + (k * lsize * msize * nsize * osize)
                    + (l * msize * nsize * osize) + (m * nsize * osize) + (n * osize) + o];
        }

        inline T& operator()(const uint32_t & i) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i);
#endif

            return data_m[i];
        }

        inline T& operator()(const uint32_t& i, const uint32_t & j) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j);
#endif

            return data_m[i * jsize + j];
        }

        inline T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k);
#endif

            return data_m[(i * jsize * ksize)+(j * ksize) + k];
        }

        inline T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l);
#endif

            return data_m[(i * jsize * ksize * lsize)+(j * ksize * lsize) + (k * lsize) + l];
        }

        inline T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l, m);
#endif

            return data_m[(i * jsize * ksize * lsize * msize)+(j * ksize * lsize * msize) + (k * lsize * msize) + (l * msize) + m];
        }

        inline T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l, m, n);
#endif

            return data_m[(i * jsize * ksize * lsize * msize * nsize)+(j * ksize * lsize * msize * nsize)
                    + (k * lsize * msize * nsize) + (l * msize * nsize) + (m * nsize) + n];
        }

        inline T& operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            CheckBounds(i, j, k, l, m, n, o);
#endif

            return data_m[(i * jsize * ksize * lsize * msize * nsize * osize)+
                    (j * ksize * lsize * msize * nsize * osize) + (k * lsize * msize * nsize * osize)
                    + (l * msize * nsize * osize) + (m * nsize * osize) + (n * osize) + o];
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
        }

        void SetBounds(INTRINSIC_BASE minb, INTRINSIC_BASE maxb) {
            std::cout << "warning atl::Array<>::" << __func__ << "not implemented for primitive types";
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

#define MAKE_VARIABLE_ARRAY_TYPE(TYPE)\
    template<>\
    void Array<atl::Variable<TYPE> >::SetBounds(TYPE minb, TYPE maxb) {\
        for (int i = 0; i < this->data_m.size(); i++) {\
            this->data_m[i].SetBounds(minb, maxb);\
        }\
    }\


    MAKE_VARIABLE_ARRAY_TYPE(float)
    MAKE_VARIABLE_ARRAY_TYPE(double)
    MAKE_VARIABLE_ARRAY_TYPE(long double)

    template<class T2, class A>
    std::ostream& operator<<(std::ostream& out, const atl::ArrayExpression< T2, A> &expr) {

        int dims = expr.Dimensions();

        std::stringstream ss;
        switch (dims) {

            case 1:
                ss << "[0]:";
                out << ss.str();
                for (int i = 0; i < expr.Size(0); i++) {
                    if ((i % 15) == 0 && i != 0) {
                        std::cout << "\n" << ss.str();
                    }
                    out << expr(i) << " ";
                }
                out << "\n";
                break;
            case 2:
                ss << "[0]:";
                out << ss.str();
                for (int i = 0; i < expr.Size(0); i++) {
                    if ((i % 15) == 0 && i != 0) {
                        std::cout << "\n" << ss.str();
                    }
                    for (int j = 0; j < expr.Size(1); j++) {
                        out << expr(i, j) << " ";
                    }
                }
                out << "\n";
                break;
                //                case 3:
                //                    for (int i = 0; i < isize; i++) {
                //                        for (int j = 0; j < jsize; j++) {
                //                            for (int k = 0; k < ksize; k++) {
                //                                data_m[(i * jsize) + (j * ksize) + k] /= expr(i, j, k);
                //                            }
                //                        }
                //                    }
                //                    break;
                //                case 4:
                //                    for (int i = 0; i < isize; i++) {
                //                        for (int j = 0; j < jsize; j++) {
                //                            for (int k = 0; k < ksize; k++) {
                //                                for (int l = 0; l < lsize; l++) {
                //                                    data_m[(i * jsize) + (j * ksize) + (k * lsize) + l] /= expr(i, j, k, l);
                //                                }
                //                            }
                //                        }
                //                    }
                //                    break;
                //
                //                case 5:
                //                    for (int i = 0; i < isize; i++) {
                //                        for (int j = 0; j < jsize; j++) {
                //                            for (int k = 0; k < ksize; k++) {
                //                                for (int l = 0; l < lsize; l++) {
                //                                    for (int m = 0; m < msize; m++) {
                //                                        data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + m] /= expr(i, j, k, l, m);
                //                                    }
                //                                }
                //                            }
                //                        }
                //                    }
                //                    break;
                //                case 6:
                //                    for (int i = 0; i < isize; i++) {
                //                        for (int j = 0; j < jsize; j++) {
                //                            for (int k = 0; k < ksize; k++) {
                //                                for (int l = 0; l < lsize; l++) {
                //                                    for (int m = 0; m < msize; m++) {
                //                                        for (int n = 0; n < nsize; n++) {
                //                                            data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + n] /= expr(i, j, k, l, m, n);
                //                                        }
                //                                    }
                //                                }
                //                            }
                //                        }
                //                    }
                //                    break;
                //
                //                case 7:
                //                    for (int i = 0; i < isize; i++) {
                //                        for (int j = 0; j < jsize; j++) {
                //                            for (int k = 0; k < ksize; k++) {
                //                                for (int l = 0; l < lsize; l++) {
                //                                    for (int m = 0; m < msize; m++) {
                //                                        for (int n = 0; n < nsize; n++) {
                //                                            for (int o = 0; o < osize; o++) {
                //                                                data_m[(i * jsize) + (j * ksize) + (k * lsize) + (l * msize) + (m * nsize) + (n * osize) + o] /= expr(i, j, k, l, m, n, o);
                //                                            }
                //                                        }
                //                                    }
                //                                }
                //                            }
                //                        }
                //                    }
                //                    break;

            default:
                break;

        }
        return out;

    }


}


#endif	/* ARRAY_HPP */

