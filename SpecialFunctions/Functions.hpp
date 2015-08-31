/* 
 * File:   Functions.hpp
 * Author: matthewsupernaw
 *
 * Created on May 31, 2013, 2:45 PM
 */

#ifndef FUNCTIONS_HPP
#define	FUNCTIONS_HPP

//#include "Distributions.hpp"
#include <string>
#include <iostream>
#include "../AutoDiff/AutoDiff.hpp"
#ifndef M_M_PI
// Source: http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html
#define M_M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406 
#endif

/**
 * 
 *@defgroup SpecialFunctions Special Functions
 */

namespace atl {
//
//    template<class T> T mfexp(const T & x) {
//        T b = T(60);
//        if (x <= b && x >= T(-1) * b) {
//            return std::exp(x);
//        } else if (x > b) {
//            return std::exp(b)*(T(1.) + T(2.) * (x - b)) / (T(1.) + x - b);
//        } else {
//            return std::exp(T(-1) * b)*(T(1.) - x - b) / (T(1.) + T(2.) * (T(-1) * x - b));
//        }
//    }


    /**
     *Special
     */

    /**
     *class SpecialFunctionException
     */
    class SpecialFunctionException : public std::exception {
    public:

        SpecialFunctionException(std::string ex) : error(ex) {

        }

        virtual ~SpecialFunctionException()throw () {

        }

        const char* what() {
            return this->error.c_str();
        }
    private:
        std::string error;

    };

    static size_t Factorial(size_t n) {
        if (n <= 1)
            return (1);
        else
            n = n * Factorial(n - 1);
        return (n);
    }

    /**
     * Binomial coefficient.
     * 
     * \f$
     * {n \choose k} = \frac{n!}{k!(n-k)!}
     * \f$
     * 
     * @param n
     * @param k
     * @return 
     */
    static size_t Choose(const size_t &n, const size_t &k) {
        
        if(k > n){
            return 0;
        }
        
        return Factorial(n)/(Factorial(k)*Factorial(n-k));
    }

    //
    //    template<class T>
    //    T LogFactorial(const T &n) {
    //        int y;
    //        T z;
    //        if (n == 1){
    //            return 0;
    //        }else {
    //            z = 0;
    //
    //            for (y = 2; y <= n; y++) {
    //                z = std::log(y) + z;
    //            }
    //
    //            return z;
    //        }
    //    }

    /**
     * Error function.
     * Reference: Handbook of Mathematical Functions by Abramowitz and Stegun. 
     * 
     * @param x
     * @return 
     */
    template<class T>
    static T Erf(const T &x) {
        // constants
        T a1 = T(0.254829592);
        T a2 = T(-0.284496736);
        T a3 = T(1.421413741);
        T a4 = T(-1.453152027);
        T a5 = T(1.061405429);
        T p = T(0.3275911);

        // Save the sign of x
        T sign = T(1);
        if (x < T(0))
            sign = T(-1);
        T xx = std::fabs(x);

        // A&S formula 7.1.26
        T t = T(1.0) / (T(1.0) + p * xx);
        T y = T(1.0) - T((((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * std::exp(T(-1.0) * xx * xx));

        return sign*y;
    }

    /**
     * Complimentary Error function.
     * Reference: Handbook of Mathematical Functions by Abramowitz and Stegun. 
     * 
     * @param x
     * @return 
     */
    template<class T>
    static T ErfC(T &x) {
        return T(1) - Erf<T > (x);
    }

    /**
     * Taylor series expansion.
     * probably a better method for this.(MS)
     * @param y
     * @return 
     */
    template<class T>
    static T InverseErf(T &y) {

        return T(.5) * sqrt(T(M_PI))*(y + ((T(M_PI) / T(12))*(y * y * y)) +
                ((T(7) * T(M_PI) * T(M_PI) / T(480))*(y * y * y * y * y))+
                ((T(127) * T(M_PI) * T(M_PI) * T(M_PI) / T(40320))*(y * y * y * y * y * y * y))+
                ((T(4369) * T(M_PI) * T(M_PI) * T(M_PI) * T(M_PI) / T(5806080))*(y * y * y * y * y * y * y * y * y))+
                ((T(34807) * T(M_PI) * T(M_PI) * T(M_PI) * T(M_PI) * T(M_PI) / T(182476800))*(y * y * y * y * y * y * y * y * y * y * y)));
    }

    /**
     *  Returns the natural logarithm of the gamma
     *  function using the lanczos approximation
     *  method.
     * 
     */
    template<class T>
    static T GammaLn(const T &xx) {

        T ret;


        int k;
        T Ag;
        T term1, term2;

        T x = xx - T(1.0);


        Ag = T(0.99999999999980993227684700473478) +
                (T(676.520368121885098567009190444019) / (x + T(1)))+
                (T(-1259.13921672240287047156078755283) / (x + T(2)))+
                (T(771.3234287776530788486528258894 / (x + T(3))))+
                (T(-176.61502916214059906584551354) / (x + T(4)))+
                (T(12.507343278686904814458936853) / (x + T(5)))+
                (T(-0.13857109526572011689554707) / (x + T(6)))+
                (T(9.984369578019570859563e-6) / (x + T(7)))+
                (T(1.50563273514931155834e-7) / (x + T(8)));

        term1 = (x + T(0.5)) * std::log((x + T(7.5)) / T(M_E));
        term2 = T(0.9189385332046727418)/*HALF LOG 2M_PI*/ + std::log(Ag);
        ret = term1 + (term2 - T(7.0));

        return ret;

    }

    /**
     * 
     * @param x
     * @return 
     */
    template<class T>
    static T Gamma(const T &x) {
        T ret;
        if (x <= T(0)) {
            return T(0);
        } else {
            ret = GammaLn<T > (x);
            return mfexp(ret);
        }

    }

    template<class T>
    T IncompleteGamma(T x, T alpha, T ln_gamma_alpha) {
        // (1) series expansion     if (alpha>x || x<=1)
        // (2) continued fraction   otherwise
        // RATNEST FORTRAN by
        // Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
        // 19: 285-287 (AS32)

        T accurate = 1e-8, overflow = 1e30;
        T factor, gin, rn, a, b, an, dif, term;
        T pn0, pn1, pn2, pn3, pn4, pn5;

        if (x == 0.0) {
            return 0.0;
        }
        if (x < 0.0 || alpha <= 0.0) {
            std::cout << ("Arguments out of bounds");
        }

        factor = std::exp(alpha * std::log(x) - x - ln_gamma_alpha);

        if (x > 1 && x >= alpha) {
            // continued fraction
            a = 1 - alpha;
            b = a + x + 1;
            term = 0;
            pn0 = 1;
            pn1 = x;
            pn2 = x + 1;
            pn3 = x * b;
            gin = pn2 / pn3;

            do {
                a++;
                b += 2;
                term++;
                an = a * term;
                pn4 = b * pn2 - an * pn0;
                pn5 = b * pn3 - an * pn1;

                if (pn5 != 0) {
                    rn = pn4 / pn5;
                    dif = std::fabs(gin - rn);
                    if (dif <= accurate) {
                        if (dif <= accurate * rn) {
                            break;
                        }
                    }

                    gin = rn;
                }
                pn0 = pn2;
                pn1 = pn3;
                pn2 = pn4;
                pn3 = pn5;
                if (std::fabs(pn4) >= overflow) {
                    pn0 /= overflow;
                    pn1 /= overflow;
                    pn2 /= overflow;
                    pn3 /= overflow;
                }
            } while (true);
            gin = 1 - factor * gin;
        } else {
            // series expansion
            gin = 1;
            term = 1;
            rn = alpha;
            do {
                rn++;
                term *= x / rn;
                gin += term;
            } while (term > accurate);
            gin *= factor / alpha;
        }
        return gin;
    }

    /**
     * 
     * @param x
     * @return 
     */
    template<class T>
    static T LogFactorial(const T &x) {
        return GammaLn(x + T(1.0));
    }

    template <class T>
    static T UpperIncompleteGamma(const T& s, const T& x) {
        static const T T_0 = T(0), T_1 = T(1), T_2 = T(2), T_3 = T(3);

        T A_prev = T_0;
        T B_prev = T_1;
        T A_cur = pow(x, s) / exp(x);
        T B_cur = x - s + T_1;

        T a = s - T_1;
        T b = B_cur + T_2;
        T n = s - T_3;

        for (;;) {
            const T A_next = b * A_cur + a * A_prev;
            const T B_next = b * B_cur + a * B_prev;

            if (A_next * B_cur == A_cur * B_next) {
                return A_cur / B_cur;
            }

            A_prev = A_cur;
            A_cur = A_next;

            B_prev = B_cur;
            B_cur = B_next;

            a = a + n;
            b = b + T_2;
            n = n - T_2;
        }
    }

    /**
     * 
     * @param a
     * @param x
     * @return 
     */
    template<class T>
    static T LowerIncompleteGamma(const T &x, const T &alpha) {


        T sum = T(0);
        T term = 1.0 / alpha;
        T n = 1;
        while (term != 0) {
            sum += term;
            term *= (x / (alpha + n));
            n++;
        }
        return std::pow(x, alpha) * std::exp(-1 * x) * sum;
    }


    //return IncompleteGamma(x,alpha ,atl::GammaLn(alpha));
    //    // (1) series expansion     if (alpha>x || x<=1)
    //    // (2) continued fraction   otherwise
    //    // RATNEST FORTRAN by
    //    // Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
    //    // 19: 285-287 (AS32)
    //
    //    T ln_gamma_alpha = std::log(Gamma(alpha));
    //    T accurate = 1e-8, overflow = 1e30;
    //    T factor, gin, rn, a, b, an, dif, term;
    //    T pn0, pn1, pn2, pn3, pn4, pn5;
    //
    //    if (x == 0.0) {
    //        return 0.0;
    //    }
    //    if (x < 0.0 || alpha <= 0.0) {
    //        throw SpecialFunctionException("SpecialFunction exception. call to IncompleteGamma(x,a,b), arguments out of bounds!");
    //    }
    //
    //    factor = std::exp(alpha * std::log(x) - x - ln_gamma_alpha);
    //
    //    if (x > 1 && x >= alpha) {
    //        // continued fraction
    //        a = 1 - alpha;
    //        b = a + x + 1;
    //        term = 0;
    //        pn0 = 1;
    //        pn1 = x;
    //        pn2 = x + 1;
    //        pn3 = x * b;
    //        gin = pn2 / pn3;
    //
    //        do {
    //            a++;
    //            b += 2;
    //            term++;
    //            an = a * term;
    //            pn4 = b * pn2 - an * pn0;
    //            pn5 = b * pn3 - an * pn1;
    //
    //            if (pn5 != 0) {
    //                rn = pn4 / pn5;
    //                dif = std::fabs(gin - rn);
    //                if (dif <= accurate) {
    //                    if (dif <= accurate * rn) {
    //                        break;
    //                    }
    //                }
    //
    //                gin = rn;
    //            }
    //            pn0 = pn2;
    //            pn1 = pn3;
    //            pn2 = pn4;
    //            pn3 = pn5;
    //            if (std::fabs(pn4) >= overflow) {
    //                pn0 /= overflow;
    //                pn1 /= overflow;
    //                pn2 /= overflow;
    //                pn3 /= overflow;
    //            }
    //        } while (true);
    //        gin = 1 - factor * gin;
    //    } else {
    //        // series expansion
    //        gin = 1;
    //        term = 1;
    //        rn = alpha;
    //        do {
    //            rn++;
    //            term *= x / rn;
    //            gin += term;
    //        } while (term > accurate);
    //        gin *= factor / alpha;
    //    }
    //    return gin;
    //}

    /**
     * beta function, a.k.a. Euler integral. 
     * @param a
     * @param b
     * @return 
     */
    template<class T>
    static T Beta(const T &a, const T &b) {
        return Gamma<T > (a) * Gamma<T > (b) / Gamma<T > (a + b);
        /*
                    T ret;
                    T aa = a;
                    T bb = b;
                    T ga;
                    T gb;
                    T gab;

                    ga = Gamma(aa);
                    gb = Gamma(bb);
                    gab = Gamma(aa + bb);
                    ret = ga * gb / gab;

                    return ret;*/
    }

    /**
     * 
     * @param a
     * @param b
     * @return 
     */
    template<class T>
    static T BetaLn(const T &a, const T &b) {

        T ret;
        T aa = a;
        T bb = b;
        T ga;
        T gb;
        T gab;

        ga = GammaLn(aa);
        gb = GammaLn(bb);
        gab = GammaLn(aa + bb);
        ret = ga + gb - gab;

        return ret;
    }

    /**
     * A continued fraction representation of the beta function.
     * 
     * @param a
     * @param b
     * @param x
     * @param max_iterations
     * @return 
     */
    template<class T>
    static T BetaContinuedFraction(const T &a, const T &b, const T &x, int max_iterations = 50) {


        T m = 1;
        T eps = T(3e-5);
        T am = T(1);
        T bm = T(1);
        T az = T(1);
        T qab = a + b;
        T qap = a + T(1);
        T qam = a - T(1);
        T bz = T(1) - qab * x / qap;
        T aold = T(0);
        T em, tem, d, ap, bp, app, bpp;

        while ((m < max_iterations) && (std::fabs(az - aold) >= eps * std::fabs(az))) {
            em = m;
            tem = em + em;
            d = em * (b - m) * x / ((qam + tem)*(a + tem));
            ap = az + d*am;
            bp = bz + d*bm;
            d = T(-1)*(a + em)*(qab + em) * x / ((a + tem)*(qap + tem));
            app = ap + d*az;
            bpp = bp + d*bz;
            aold = az;
            am = ap / bpp;
            bm = bp / bpp;
            az = app / bpp;
            bz = 1;
            m++;
        }
        return az;
    }

    /**
     * 
     * @param x
     * @param a
     * @param b
     * @return 
     */
    template<class T>
    static T IncompleteBeta(const T &x, const T &a, const T &b) {

        // the incomplete beta function from 0 to x with parameters a, b
        // x must be in (0,1) (else returns error)

        if (x < T(0) || x > T(1.0)) {
            throw SpecialFunctionException("SpecialFunction exception. call to IncompleteBeta(x,a,b), x must be in range of 0-1!");
        }

        T er = T(0);
        T bt = T(0); //, 
        T beta = T(0);
        if (x == T(0) || x == T(1)) {
            bt = 0;
        } else if ((x > T(0)) && (x < T(1))) {
            bt = Gamma<T > (a + b) * std::pow(x, a) * std::pow(T(1) - x, b) / (Gamma<T > (a) * Gamma<T > (b));
        }
        if (x < (a + T(1)) / (a + b + T(2))) {
            beta = bt * BetaContinuedFraction<T > (a, b, x) / a;
        } else {
            beta = T(1) - bt * BetaContinuedFraction<T > (b, a, T(1) - x) / b;
        }
        return beta;
    }

    template<class T>
    static T RegularizedBeta(const T &x, const T &a, const T &b) {

        T result;
        T rr;
        T bf;
        T lnb;
        // BetaContinuedFractionDoubleRef bcf = new BetaContinuedFraction<double>(a, b);
        T ret = T(0);
      

        if (/*isNaN<T > (x) || isNaN<T > (a) || isNaN<T > (b) || */(x < T(0)) || (x > T(1)) || (a <= T(0.0)) || (b <= T(0.0))) {
            //  ret = T(std::numeric_limits<double>::quiet_NaN());
            //            result = T(0);
            return ret;

        } else if (x > (a + T(1.0)) / (a + b + T(2.0))) {
            rr = RegularizedBeta<T > (T(1.0) - x, b, a);
            result = (T(1.0) - rr);
            return result;
        } else {

            bf = BetaContinuedFraction<T > (a, b, x);
            lnb = BetaLn<T > (a, b);
            // ret = exp((a * log(x)) + (b * log(1.0 - x)) - log(a) - lnb->getValue()) * 1.0 / bf->getValue();
            ret = std::exp((a * std::log(x)) + (b * std::log(T(1.0) - x)) -
                    std::log(a) - BetaLn<T > (a, b)) *
                    T(1.0) / BetaContinuedFraction<T > (a, b, x);

            return ret;
        }

    }

    template<class T>
    static T IncompleteReqularizedBeta(const T &x, const T &a, const T &b) {

        return atl::IncompleteBeta<T > (x, a, b) / atl::RegularizedBeta<T > (T(1), a, b);
    }



    




}//atl

#endif	/* FUNCTIONS_HPP */

