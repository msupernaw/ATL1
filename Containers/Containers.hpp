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


    //    template<class T2, class A>
    //    const typename atl::PromoteType<typename VectorExpression<T2, A>::RET_TYPE, typename VectorExpression<T2, A>::RET_TYPE >::return_type Norm2(const VectorExpression<T2, A> &expr) {
    //        typedef typename atl::PromoteType<typename VectorExpression<T2, A>::RET_TYPE, typename VectorExpression<T2, A>::RET_TYPE >::return_type RET;
    //       RET ret; // = TT(0.0);
    //        size_t s = expr.Size(0);
    //        for (int i = 0; i < s; i++) {
    //
    //            ret += expr(i) * expr(i);
    //        }
    //        return ret;
    //    }
    //    

   

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
        for(int i =start; i<end; i++){
            ret+= expr(i)*expr(i);
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
            ret = expr(0)*expr(0);
            for (size_t i = 1UL; i < end; i += 2UL) {
                ret += expr(i)*expr(i) + expr(i + 1)*expr(i + 1);
            }

            if (end < s) {
                ret += expr(end)*expr(end);
            }
        }
        return (ret);
    }
    
    template<class T2, class A>
    inline const typename A::BASE_TYPE Norm(const VectorExpression<T2, A> &expr, bool concurrent = false) {
        return std::sqrt(Norm2(expr,concurrent));
    }
    
    template<class T2, class A>
    inline void Sum(typename A::BASE_TYPE& ret, const VectorExpression<T2, A> &expr, bool concurrent = false) {
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
    }
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
}

#endif	/* CONTAINERS_HPP */

