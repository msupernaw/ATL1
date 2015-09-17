/* 
 * File:   MatrixArrayOperators.hpp
 * Author: matthewsupernaw
 *
 * Created on October 15, 2014, 3:32 PM
 */

#ifndef MATRIXARRAYOPERATORS_HPP
#define	MATRIXARRAYOPERATORS_HPP
#include "ArrayExpressionBase.hpp"
#include "MatrixExpressionBase.hpp"
#include "../Traits/Promote.hpp"

namespace atl {
#ifdef ATL_CONTAINER_EXPERIMENTAL
    template< class T, class A>
    struct MatrixArrayExpression {
        typedef T RET_TYPE;

        const A & Cast() const {
            return static_cast<const A&> (*this);
        }

        inline const size_t Size(const int32_t & dimension) const {
            return Cast().Size(dimension);
        }

        inline const size_t Dimensions() const {
            return Cast(). Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return Cast().operator()(i, j);
        }

        inline void ExpressionLength(uint32_t& length)const {
            Cast().ExpressionLength(length);
        }

    };

    template< class LHS, class RHS>
    struct MatrixArrayAdd : MatrixArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixArrayAdd<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixArrayAdd(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            for (int i = 0; i < 7; i++) {
                assert(lhs_m.Size(i) == rhs_m.Size(i));
            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) + rhs_m(i, j);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }



    };

    template< class LHS, class RHS>
    struct MatrixArraySubtract : MatrixArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixArraySubtract<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixArraySubtract(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            for (int i = 0; i < 7; i++) {
                assert(lhs_m.Size(i) == rhs_m.Size(i));
            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) - rhs_m(i, j);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }


    };

    template< class LHS, class RHS>
    struct MatrixArrayMultiply : MatrixExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixArrayMultiply<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::intrinsic_trait<typename LHS::RET_TYPE> INTR_TRAITL;
        typedef typename atl::intrinsic_trait<typename RHS::RET_TYPE> INTR_TRAITR;

        const LHS& lhs_m;
        const RHS& rhs_m;
        size_t lrows;
        size_t lcols;
        size_t end;

        inline explicit MatrixArrayMultiply(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {
            lrows = lhs_m.Size(0);
            lcols = lhs_m.Size(1);
            end = (((lcols - 1UL) & size_t(-2)) + 1UL);
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            //                        for (int i = 0; i < 2; i++) {
            //                            assert(lhs_m.Size(i) == rhs_m.Size(i));
            //                        }
            assert(lhs_m.Size(1) == rhs_m.Size(0));
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {

            if (lhs_m.Size(1) == rhs_m.Size(0)) {

                switch (dimension) {
                    case 0:
                        return lhs_m.Size(0);
                    case 1:
                        return rhs_m.Size(1);
                    default:
                        return 0;
                }
            }
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            assert(i < lhs_m.Size(0));
            assert(j < lhs_m.Size(1));
#endif   
            RET_TYPE ret;


            if (0 != lrows) {
#ifdef ATL_HAS_SSE

                if (lcols < 10) {
                    goto compute;
                }
                if (INTR_TRAITL::is_intrinsic && INTR_TRAITR::is_intrinsic) {

                    if (INTR_TRAITL::size == 64 || INTR_TRAITR::size == 64) {

                        __declspec(align(16))double x[2];
                        __declspec(align(16)) double y[2];
                        int k;

                        __m128d X, Y; // 128-bit values
                        __m128d acc = _mm_setzero_pd(); // set to (0, 0)
                        double inner_prod, temp[4];
                        inner_prod = 0;
                        k = 0;
                        for (k = 0; k < lcols - 2; k += 2) {

                            x[0] = lhs_m(i, k);
                            x[1] = lhs_m(i, k + 1UL);


                            y[0] = rhs_m(k, j);
                            y[1] = rhs_m(k + 1UL, j);



                            X = _mm_load_pd(&x[0]); // load chunk of 2 doubles
                            Y = _mm_load_pd(&y[0]);
                            acc = _mm_add_pd(acc, _mm_mul_pd(X, Y));
                            _mm_store_pd(&temp[0], acc); // store acc into an array of floats

                        }
                        inner_prod += temp[0] + temp[1]; // + temp[2] + temp[3];
                        // add the remaining values
                        for (; k < lcols; k++) {
                            inner_prod += lhs_m(i, k) * rhs_m(k, j);
                        }

                        return inner_prod;

                    } else if (INTR_TRAITL::size == 32 || INTR_TRAITR::size == 32) {

                        __declspec(align(16)) float x[4];
                        __declspec(align(16))float y[4];
                        int k;

                        __m128 X, Y; // 128-bit values
                        __m128 acc = _mm_setzero_ps(); // set to (0, 0, 0, 0)
                        float inner_prod, temp[4];
                        inner_prod = 0;
                        k = 0;
                        for (k = 0; k < lcols - 4; k += 4) {

                            x[0] = lhs_m(i, k);
                            x[1] = lhs_m(i, k + 1UL);
                            x[2] = lhs_m(i, k + 2UL);
                            x[3] = lhs_m(i, k + 3UL);

                            y[0] = rhs_m(k, j);
                            y[1] = rhs_m(k + 1UL, j);
                            y[2] = rhs_m(k + 2UL, j);
                            y[3] = rhs_m(k + 3UL, j);


                            X = _mm_load_ps(&x[0]); // load chunk of 4 floats
                            Y = _mm_load_ps(&y[0]);
                            acc = _mm_add_ps(acc, _mm_mul_ps(X, Y));
                            _mm_store_ps(&temp[0], acc); // store acc into an array of floats

                        }
                        inner_prod += temp[0] + temp[1] + temp[2] + temp[3];
                        // add the remaining values
                        for (; k < lcols; k++) {
                            inner_prod += lhs_m(i, k) * rhs_m(k, j);
                        }

                        return inner_prod;
                    } else {
                        goto compute;
                    }
                } else {
#endif

compute:
                    ret = lhs_m(i, 0) * rhs_m(0, j);


                    for (size_t k = 1UL; k < end; k += 2UL) {
                        ret += lhs_m(i, k) * rhs_m(k, j);
                        ret += lhs_m(i, k + 1UL) * rhs_m(k + 1UL, j);
                    }
                    if (end < lcols) {
                        ret += lhs_m(i, end) * rhs_m(end, j);
                    }
#ifdef ATL_HAS_SSE
                }
#endif
            }

            return ret;
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }

    };

    template< class LHS, class RHS>
    struct MatrixArrayDivide : MatrixExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixArrayDivide<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixArrayDivide(const LHS& lhs, const RHS & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {


#ifdef ATL_ENABLE_BOUNDS_CHECKING
            for (int i = 0; i < 2; i++) {
                assert(lhs_m.Size(i) == rhs_m.Size(i));
            }
#endif

        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return lhs_m(i, j) / rhs_m(i, j);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }



    };

    template <class LHS, class RHS>
    inline const MatrixArrayAdd< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArrayAdd< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArrayAdd< LHS, RHS> operator+(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArrayAdd< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArraySubtract< LHS, RHS> operator-(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArraySubtract< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArraySubtract< LHS, RHS> operator-(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArraySubtract< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArrayMultiply< LHS, RHS> operator*(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArrayMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArrayMultiply< LHS, RHS> operator*(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArrayMultiply< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArrayDivide< LHS, RHS> operator/(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const ArrayExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArrayDivide< LHS, RHS > (a.Cast(), b.Cast());
    }

    template <class LHS, class RHS>
    inline const MatrixArrayDivide< LHS, RHS> operator/(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixArrayDivide< LHS, RHS > (a.Cast(), b.Cast());
    }

    template<class T2, class A>
    std::ostream& operator<<(std::ostream& out, const atl::MatrixArrayExpression< T2, A> &expr) {

        int dims = expr.Dimensions();
        std::stringstream ss;


        out << ss.str();
        for (int i = 0; i < expr.Size(0); i++) {
            for (int j = 0; j < expr.Size(1); j++) {

                out << expr(i, j) << " ";
            }
            out << "\n";
        }


        return out;
    }
#endif
}

#endif	/* MATRIXARRAYOPERATORS_HPP */

