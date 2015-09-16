/* 
 * File:   ConcurrentOperators.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 3:33 PM
 */

#ifndef CONCURRENTOPERATORS_HPP
#define	CONCURRENTOPERATORS_HPP

//#warning need to add friend relationships to matrix, vector, and array classes
#include "ArrayExpressionBase.hpp"
#include "Array.hpp"
#include "Matrix.hpp"
#include "MatrixMultiply.hpp"
#include "MatrixAdd.hpp"
#include <thread>
#include <vector>
#include "../Traits/Type.hpp"

namespace atl {

    /**
     * Array operations
     */
    template<class T, class T2, class RHS>
    void ConcurrentOpEquals(Array<T>& v, const ArrayExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpPlusEquals(Array<T>& v, const ArrayExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpMinusEquals(Array<T>& v, const ArrayExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpTimesEquals(Array<T>& v, const ArrayExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpDivideEquals(Array<T>& v, const ArrayExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    /**
     * Vector operations
     */

    template<class T, class T2, class RHS>
    void ConcurrentOpEquals(Vector<T>& v, const VectorExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpPlusEquals(Vector<T>& v, const VectorExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpMinusEquals(Vector<T>& v, const VectorExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpTimesEquals(Vector<T>& v, const VectorExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpDivideEquals(Vector<T>& v, const VectorExpression<T2, RHS>& expr) {
        std::cout << __FILE__ << ":" << __LINE__ << " not yet implemented!" << std::endl;
        exit(0);
    }

    /**
     * Matrix operations
     */


    template<class T, class T2, class RHS>
    void ConcurrentOpEqualsThread(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr, int start, int end) {

        for (int i = start; i < end; i++) {
            for (int j = 0; j < v.jsize; j++) {
                v.data_m[(i * v.jsize) + j] = expr(i, j);
            }
        }
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpEquals(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr) {

        std::vector<std::thread> pool;
        v.isize = expr.Size(0);
        v.jsize = expr.Size(1);
        int range = v.isize / std::thread::hardware_concurrency();
        v.data_m.resize(v.isize * v.jsize);

        for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
            if (i < (std::thread::hardware_concurrency() - 1)) {
                pool.push_back(std::thread(ConcurrentOpEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, (i + 1) * range));
            } else {
                pool.push_back(std::thread(ConcurrentOpEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, v.isize));

            }
        }

        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }


    }

    template<class T, class T2, class RHS>
    void ConcurrentOpPlusEqualsThread(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr, int start, int end) {

        for (int i = start; i < end; i++) {
            for (int j = 0; j < v.jsize; j++) {
                v.data_m[(i * v.jsize) + j] += expr(i, j);
            }
        }

    }

    template<class T, class T2, class RHS>
    void ConcurrentOpPlusEquals(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr) {
        std::vector<std::thread> pool;
        v.isize = expr.Size(0);
        v.jsize = expr.Size(1);
        int range = v.isize / std::thread::hardware_concurrency();
        v.data_m.resize(v.isize * v.jsize);

        for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
            if (i < (std::thread::hardware_concurrency() - 1)) {
                pool.push_back(std::thread(ConcurrentOpPlusEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, (i + 1) * range));
            } else {
                pool.push_back(std::thread(ConcurrentOpPlusEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, v.isize));
            }
        }

        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }

    }

    template<class T, class T2, class RHS>
    void ConcurrentOpMinusEqualsThread(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr, int start, int end) {

        for (int i = start; i < end; i++) {
            for (int j = 0; j < v.jsize; j++) {
                v.data_m[(i * v.jsize) + j] += expr(i, j);
            }
        }

    }

    template<class T, class T2, class RHS>
    void ConcurrentOpMinusEquals(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr) {
        std::vector<std::thread> pool;
        v.isize = expr.Size(0);
        v.jsize = expr.Size(1);
        int range = v.isize / std::thread::hardware_concurrency();
        v.data_m.resize(v.isize * v.jsize);

        for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
            if (i < (std::thread::hardware_concurrency() - 1)) {
                pool.push_back(std::thread(ConcurrentOpMinusEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, (i + 1) * range));
            } else {
                pool.push_back(std::thread(ConcurrentOpMinusEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, v.isize));
            }
        }

        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpTimesEqualsThread(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr, int start, int end) {
//#warning function needs review
        for (int i = start; i < end; i++) {
            for (int j = 0; j < v.jsize; j++) {
                v.data_m[(i * v.jsize) + j] = (expr * v)(i, j);
            }
        }

    }

    template<class T, class T2, class RHS>
    void ConcurrentOpTimesEquals(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr) {
        std::vector<std::thread> pool;
        v.isize = expr.Size(0);
        v.jsize = expr.Size(1);
        int range = v.isize / std::thread::hardware_concurrency();
        v.data_m.resize(v.isize * v.jsize);

        for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
            if (i < (std::thread::hardware_concurrency() - 1)) {
                pool.push_back(std::thread(ConcurrentOpEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, (i + 1) * range));
            } else {
                pool.push_back(std::thread(ConcurrentOpEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, v.isize));
            }
        }

        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }
    }

    template<class T, class T2, class RHS>
    void ConcurrentOpDivideEqualsThread(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr, int start, int end) {

        for (int i = start; i < end; i++) {
            for (int j = 0; j < v.jsize; j++) {
                v.data_m[(i * v.jsize) + j] /= expr(i, j);
            }
        }

    }

    template<class T, class T2, class RHS>
    void ConcurrentOpDivideEquals(Matrix<T>& v, const MatrixExpression<T2, RHS>& expr) {
        std::vector<std::thread> pool;
        v.isize = expr.Size(0);
        v.jsize = expr.Size(1);
        int range = v.isize / std::thread::hardware_concurrency();
        v.data_m.resize(v.isize * v.jsize);

        for (int i = 0; i < std::thread::hardware_concurrency(); i++) {
            if (i < (std::thread::hardware_concurrency() - 1)) {
                pool.push_back(std::thread(ConcurrentOpDivideEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, (i + 1) * range));
            } else {
                pool.push_back(std::thread(ConcurrentOpDivideEqualsThread<T, T2, RHS>, std::ref(v), std::ref(expr), i*range, v.isize));
            }
        }

        for (int i = 0; i < pool.size(); i++) {
            pool[i].join();
        }
    }


}



#endif	/* CONCURRENTOPERATORS_HPP */

