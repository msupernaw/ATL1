/* 
 * File:   Containers.hpp
 * Author: matthewsupernaw
 *
 * Created on December 3, 2014, 10:09 AM
 */

#ifndef CONTAINERS_HPP
#define	CONTAINERS_HPP
#include "Array.hpp"
#include "Vector.hpp"
#include "VectorAdd.hpp"
#include "VectorSubtract.hpp"
#include "VectorMultiply.hpp"
#include "VectorDivide.hpp"
#include "Matrix.hpp"
#include "ConcurrentOperators.hpp"
#include <thread>
#include <vector>
#include "../Traits/Type.hpp"

namespace atl {

    template<class T,  class M>
    struct MatrixRow : atl::VectorExpression< T, MatrixRow<T, M> > {
        const M& m_m;
        size_t row_m;

        typedef typename M::RET_TYPE RET_TYPE;
        typedef typename M::BASE_TYPE BASE_TYPE;

        MatrixRow(const atl::MatrixExpression<T, M>& m, size_t row) : m_m(m.Cast()), row_m(row) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return m_m.Size(1);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(i < m_m.Size(1));
#endif   

            return m_m(row_m, i);

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return m_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return m_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return m_m.AtRaw(row_m, i);
        }

        /**
         * Used to determine if the resulting container is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            m_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            m_m.ExpressionLength(length);
        }

    };

    template<class T, class M>
    const MatrixRow<T, M> Row(const atl::MatrixExpression<T, M>& m, size_t row) {
        return MatrixRow<T, M>(m, row);
    }
    
    
    
    template<class T, class M>
    struct MatrixColumn : atl::VectorExpression< T, MatrixColumn<T, M> > {
        const M& m_m;
        size_t col_m;

        typedef typename M::RET_TYPE RET_TYPE;
        typedef typename M::BASE_TYPE BASE_TYPE;

        MatrixColumn(const atl::MatrixExpression<T, M>& m, size_t col) : m_m(m.Cast()), col_m(col) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return m_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(i < m_m.Size(0));
#endif   

            return m_m(i,col_m);

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return m_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return m_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return m_m.AtRaw(i,col_m);
        }

        /**
         * Used to determine if the resulting container is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            m_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            m_m.ExpressionLength(length);
        }

    };

    template<class T, class M>
    const MatrixColumn<T, M> Column(const atl::MatrixExpression<T, M>& m, size_t row) {
        return MatrixColumn<T, M>(m, row);
    }


    template<class T2, class A>
    void SumThread(typename A::BASE_TYPE& ret, int start, int end, const VectorExpression<T2, A> &expr) {
        //        size_t s = end;
        typedef typename A::BASE_TYPE RET_TYPE;
        RET_TYPE r;

        size_t e = (((end - 1UL) & size_t(-2)) + 1UL);
        r = expr(start);
        for (size_t i = start + 1; i < e; i += 2UL) {
            r += expr(i) + expr(i + 1);
        }

        if (e < end) {
            r += expr(end);
        }

        ret = r;
    }

    template<class T2, class A>
    inline const typename A::BASE_TYPE Sum(const VectorExpression<T2, A> &expr, bool concurrent = false) {
        typename A::BASE_TYPE ret; // = TT(0.0);


        if (concurrent) {

            int range = expr.Size(0) / std::thread::hardware_concurrency();
            std::vector<std::thread> pool;
            std::vector<typename A::BASE_TYPE > temp(std::thread::hardware_concurrency());

            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
                if (i < (std::thread::hardware_concurrency() - 1)) {
                    //                    std::cout << i*range << " - " << (i + 1) * range << "\n";
                    pool.push_back(std::thread(SumThread<T2, A>, std::ref(temp[i]), i*range, (i + 1) * range, std::ref(expr)));
                } else {
                    pool.push_back(std::thread(SumThread<T2, A>, std::ref(temp[i]), i*range, expr.Size(0), std::ref(expr)));
                    //                    std::cout << i*range << " - " << expr.Size(0) << "\n";
                }
            }

            for (int i = 0; i < pool.size(); i++) {
                pool[i].join();
            }

            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
                ret += temp[i];
            }

        } else {

            size_t s = expr.Size(0);
            size_t end = (((s - 1UL) & size_t(-2)) + 1UL);
            ret = expr(0);
            for (size_t i = 1UL; i < end; i += 2UL) {
                ret += expr(i) + expr(i + 1);
            }

            if (end < s) {
                ret += expr(end);
            }
        }
        return ret;
    }

    template<class T2, class A>
    void Norm2Thread(typename A::BASE_TYPE& ret, int start, int end, const VectorExpression<T2, A> &expr) {
        //        size_t s = end;
        typedef typename A::BASE_TYPE RET_TYPE;
        //        RET_TYPE r;

#pragma unroll
        for (int i = start; i < end; i++) {
            ret += expr(i) * expr(i);
        }
        //        size_t e = (((end - 1UL) & size_t(-2)) + 1UL);
        //        ret = expr(start)*expr(start);
        //        for (size_t i = start + 1; i < e; i += 2UL) {
        //            ret += expr(i)*expr(i) + expr(i + 1)*expr(i + 1);
        //        }
        //
        //        if (e < end) {
        //            ret += expr(e)*expr(e);
        //        }

        //        ret = r;
    }

    template<class T2, class A>
    inline const typename A::BASE_TYPE Norm2(const VectorExpression<T2, A> &expr, bool concurrent = false) {
        typename A::BASE_TYPE ret; // = TT(0.0);


        if (concurrent) {

            int range = expr.Size(0) / std::thread::hardware_concurrency();
            std::vector<std::thread> pool;
            std::vector<typename A::BASE_TYPE > temp(std::thread::hardware_concurrency());

            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
                if (i < (std::thread::hardware_concurrency() - 1)) {
                    //                    std::cout << i*range << " - " << (i + 1) * range << "\n";
                    pool.push_back(std::thread(Norm2Thread<T2, A>, std::ref(temp[i]), i*range, (i + 1) * range, std::ref(expr)));
                } else {
                    pool.push_back(std::thread(Norm2Thread<T2, A>, std::ref(temp[i]), i*range, expr.Size(0), std::ref(expr)));
                    //                    std::cout << i*range << " - " << expr.Size(0) << "\n";
                }
            }

            for (int i = 0; i < pool.size(); i++) {
                pool[i].join();
            }

            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
                ret += temp[i];
            }

        } else {

            size_t s = expr.Size(0);
            size_t end = (((s - 1UL) & size_t(-2)) + 1UL);
            ret = expr(0) * expr(0);
            for (size_t i = 1UL; i < end; i += 2UL) {
                ret += expr(i) * expr(i) + expr(i + 1) * expr(i + 1);
            }

            if (end < s) {
                ret += expr(end) * expr(end);
            }
        }
        return (ret);
    }

       template<class T2, class A>
    void Norm2ThreadA(typename A::BASE_TYPE& ret, int start, int end, const ArrayExpression<T2, A> &expr) {
        //        size_t s = end;
        typedef typename A::BASE_TYPE RET_TYPE;
        //        RET_TYPE r;

#pragma unroll
        for (int i = start; i < end; i++) {
            ret += expr(i) * expr(i);
        }
        //        size_t e = (((end - 1UL) & size_t(-2)) + 1UL);
        //        ret = expr(start)*expr(start);
        //        for (size_t i = start + 1; i < e; i += 2UL) {
        //            ret += expr(i)*expr(i) + expr(i + 1)*expr(i + 1);
        //        }
        //
        //        if (e < end) {
        //            ret += expr(e)*expr(e);
        //        }

        //        ret = r;
    }

    template<class T2, class A>
    inline const typename A::BASE_TYPE Norm2(const ArrayExpression<T2, A> &expr, bool concurrent = false) {
        typename A::BASE_TYPE ret; // = TT(0.0);


        if (concurrent) {

            int range = expr.Size(0) / std::thread::hardware_concurrency();
            std::vector<std::thread> pool;
            std::vector<typename A::BASE_TYPE > temp(std::thread::hardware_concurrency());

            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
                if (i < (std::thread::hardware_concurrency() - 1)) {
                    //                    std::cout << i*range << " - " << (i + 1) * range << "\n";
                    pool.push_back(std::thread(Norm2ThreadA<T2, A>, std::ref(temp[i]), i*range, (i + 1) * range, std::ref(expr)));
                } else {
                    pool.push_back(std::thread(Norm2ThreadA<T2, A>, std::ref(temp[i]), i*range, expr.Size(0), std::ref(expr)));
                    //                    std::cout << i*range << " - " << expr.Size(0) << "\n";
                }
            }

            for (int i = 0; i < pool.size(); i++) {
                pool[i].join();
            }

            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
                ret += temp[i];
            }

        } else {

            size_t s = expr.Size(0);
            size_t end = (((s - 1UL) & size_t(-2)) + 1UL);
            ret = expr(0) * expr(0);
            for (size_t i = 1UL; i < end; i += 2UL) {
                ret += expr(i) * expr(i) + expr(i + 1) * expr(i + 1);
            }

            if (end < s) {
                ret += expr(end) * expr(end);
            }
        }
        return (ret);
    }

    
    
    template<class T2, class A>
    inline const typename A::BASE_TYPE Norm(const VectorExpression<T2, A> &expr, bool concurrent = false) {
        return std::sqrt(Norm2(expr, concurrent));
    }

//    template<class T2, class A>
//    inline typename A::BASE_TYPE Sum(const VectorExpression<T2, A> &expr){
//        typedef typename A::BASE_TYPE RET_TYPE;
//        RET_TYPE ret;
//        for(int i = 0; i < expr.Size(0); i++){
//            ret+=expr(i);
//        }
//        return ret;
//    }
    
    
    template<class T2, class A>
    inline const atl::Vector<typename A::BASE_TYPE> Square(const VectorExpression<T2, A> &expr){
       
        atl::Vector<typename A::BASE_TYPE> ret(expr.Size(0));
        for(int i = 0; i < expr.Size(0); i++){
            ret(i)=expr(i)*expr(i);
        }
        return ret;
    }
    
    
//    
//    template<class T2, class A>
//    inline void Sum(typename A::BASE_TYPE& ret, const VectorExpression<T2, A> &expr, bool concurrent = false) {
        //        typename A::BASE_TYPE ret; // = TT(0.0);


        //
        //        if (concurrent) {
        //
        //            int range = expr.Size(0) / std::thread::hardware_concurrency();
        //            std::vector<std::thread> pool;
        //            std::vector<typename A::BASE_TYPE > temp(std::thread::hardware_concurrency());
        //
        //            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
        //                if (i < (std::thread::hardware_concurrency() - 1)) {
        //                    std::cout << i*range << " - " << (i + 1) * range << "\n";
        //                    pool.push_back(std::thread(SumThread<T2, A>, std::ref(temp[i]), i*range, (i + 1) * range, std::ref(expr)));
        //                } else {
        //                    pool.push_back(std::thread(SumThread<T2, A>, std::ref(temp[i]), i*range, expr.Size(0), std::ref(expr)));
        //                    std::cout << i*range << " - " << expr.Size(0) << "\n";
        //                }
        //            }
        //
        //            for (int i = 0; i < pool.size(); i++) {
        //                pool[i].join();
        //            }
        //
        //            for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
        //                ret += temp[i];
        //            }
        //
        //        } else {
        //            
        //            size_t s = expr.Size(0);
        //            size_t end = (((s - 1UL) & size_t(-2)) + 1UL);
        //            ret = expr(0);
        //            for (size_t i = 1UL; i < end; i += 2UL) {
        //                ret += expr(i) + expr(i + 1);
        //            }
        //
        //            if (end < s) {
        //                ret += expr(end);
        //            }
        //        }
        //        return ret;
//    }
    ////    
    //     template<atl::Variable<double>, class A>
    //      atl::Variable<double> Sum(const VectorExpression<atl::Variable<double>, A> &expr) {
    //        atl::Variable<double> ret; // = TT(0.0);
    //        size_t s = expr.Size(0);
    //        for (int i = 0; i < s; i++) {
    //            ret += expr(i);
    //        }
    //        return ret;
    //    }
    
    template<class T>
    const atl::Matrix<typename T::RET_TYPE> Transpose(const atl::VectorExpression<typename T::RET_TYPE, T>& v ){
        return atl::Matrix<typename T::RET_TYPE>(v);
    }
}

#endif	/* CONTAINERS_HPP */

