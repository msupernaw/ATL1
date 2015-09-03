/* 
 * File:   Distribution.hpp
 * Author: matthewsupernaw
 *
 * Created on June 3, 2013, 8:46 PM
 */

#ifndef Distribution_HPP
#define	Distribution_HPP

#include <cmath>
#include "../SpecialFunctions/Functions.hpp"


#define T_NAN (T(0)/T(0))

#ifdef  DEFAULT_EPSILON 

#undef  DEFAULT_EPSILON

#define DEFAULT_EPSILON 10e-7

#else

#define DEFAULT_EPSILON 10e-7

#endif

#ifndef MAX_ITERATIONS
#define MAX_ITERATIONS 1000
#endif
/**
 * @defgroup Distributions
 * 
 * @Distributions Every distribution in ATL has four functions and a class derived from
 * class atl::DistributionBase.<br> Distribution function root names are preceeded 
 * by a letter indicating the function type:<br><br>
 * <ul>
 * <li><b>p</b> for "probability", the cumulative distribution function(cdf)</li>
 * <li><b>q</b> for "quantile", the inverse of the the cumulative distribution function(cdf)</li>
 * <li><b>d</b> for "density", the probability density function</li>
 * <li><b>r</b> for "random", a random variable from the specified distribution</li>
 * </ul>
 * 
 *<br>Available Functions:<br><br>
 *<br> Distribution	            Functions                      Class name
 *<br>Beta              pbeta    qbeta	  dbeta	  rbeta         BetaDistribution
 *<br>Binomial          pbinom	 qbinom	  dbinom  rbinom        BinomialDistribution
 *<br>Cauchy            pcauchy	 qcauchy  dcauchy rcauchy       CauchyDistribution
 *<br>Chi-Square	pchisq	 qchisq	  dchisq  rchisq        ChiSquaredDistribution
 *<br>Exponential	pexp	 qexp	  dexp	  rexp          ExponentialDistribution
 *<br>F                 pf	 qf	  df	  rf            FDistribution
 *<br>Gamma             pgamma	 qgamma	  dgamma  rgamma        GammaDistribution
 *<br>Geometric         pgeom	 qgeom	  dgeom	  rgeom         GeometricDistribution
 *<br>Hypergeometric	phyper	 qhyper	  dhyper  rhyper        HypergeometricDistribution
 *<br>Logistic          plogis	 qlogis	  dlogis  rlogis        LogisticDistribution             
 *<br>Log Normal	plnorm	 qlnorm	  dlnorm  rlnorm        LogNormalDistribution
 *<br>Negative Binomial	pnbinom	 qnbinom  dnbinom rnbinom       NegativeBinomialDistribution
 *<br>Normal            pnorm	 qnorm	  dnorm	  rnorm         NormalDistribution
 *<br>Poisson           ppois	 qpois	  dpois	  rpois         PoissonDistribution     
 *<br>Student's t         pt	 qt	  dt	   rt           StudentsTDistribution
 *<br>Uniform           punif	 qunif	  dunif	   runif        UniformDistribution
 *<br>Weibull           pweibull qweibull dweibull rweibull     WeibullDistribution
 *<br>
 *<br>
 * 
 * @ingroup Statistics
 */



namespace atl {

    /**
     * @ingroup Distributions
     *
     *  @brief An abstract class where polymorphism is desired. Classes that derive from
     * this class overload the pure virtual functions:<br><br>
     *  const T Probability(const T &x) const<br>
     *  const T Cumulative(const T &x) const<br>
     *  const T Inverse(const T &probability) const<br>
     *  const T Random() const<br><br>
     * 
     * @param x
     * @return 
     */
    template<class T>
    class DistributionBase {
    public:

        /**
         * The probability density function for this distribution.
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        virtual const T Probability(const T &x) const = 0;

        /**
         * The cumulative distribution function for this distribution.
         * 
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        virtual const T Cumulative(const T &x) const = 0;

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        virtual const T Inverse(const T &probability) const = 0;

        /**
         * Random number.
         * 
         * @return a random deviate from this distribution.
         */
        virtual const T Random() const = 0;




    };


    /**
     *@defgroup NormalDistribution Normal Distribution
     * 
     * 
     * 
     *@ingroup Distributions  
     *
     */

    /**
     * 
     * @ingroup NormalDistribution
     * 
     * @brief Cumulative distribution function for a standard normal distribution.
     *
     * \f$
     * f(x) = \frac{1}{2} [1+ erf(\frac{x - \mu}{\sqrt{2\theta^2}})]
     * \f$
     * 
     * <b>Where,</b><br>
     *  \f$\ \mu \f$  is the mean.<br>
     *  \f$\ \theta\f$ is the standard deviation.<br>
     * 
     *  @image html http://upload.wikimedia.org/wikipedia/commons/c/ca/Normal_Distribution_CDF.svg
     * 
     *  <br> <b>Source:</b> http://en.wikipedia.org/wiki/Normal_distribution<br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013
     *   
     * @param x
     * @param mean
     * @param sd 
     * @return the probability that a stochastic variable x is less then X.
     */
    template<class T>
    const T pnorm(const T &x, const T &mean, const T &sd) {
        return T(0.5)*(T(1.0) + atl::Erf<T > ((x - mean) / std::sqrt(T(2.0) * sd * sd)));
    }

    /**
     * \ingroup NormalDistribution
     * 
     * @brief The probability density function for a standard normal distribution.
     *
     * \f$
     * f(x) = \frac{1}{\sqrt{2\pi\theta^2}}exp(-\frac{(x-\mu)^2}{2\theta^2})
     * \f$<br>
     * <b>Where,</b><br>
     *  \f$\ \mu \f$  is the mean.<br>
     *  \f$\ \theta\f$ is the standard deviation.<br>
     * 
     *  @image html http://upload.wikimedia.org/wikipedia/commons/7/74/Normal_Distribution_PDF.svg
     * 
     *  <br> <b>Source:</b> http://en.wikipedia.org/wiki/Normal_distribution<br>
     * 
     * 
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013
     *  
     * @param x
     * @param mean
     * @param sd
     * @return relative likelihood for this random variable to have value x.
     */
    template<class T>
    const T dnorm(const T &x, const T &mean, const T &sd) {
        return (T(1.0) / (sd * std::sqrt(T(2.0) * T(M_PI))))*
                std::exp((T(-1.0)*(x - mean)*(x - mean))
                / (T(2.0) * sd * sd));
    }

    template<class T>
    size_t Size(const T &x) {
        return sizeof (T);
    }

    /**
     * \ingroup NormalDistribution
     * 
     * 
     * @brief Inverse Cumulative distribution function for a standard normal distribution.
     *
     * <br><b>Source</b> http://en.wikipedia.org/wiki/Normal_distribution<br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013
     * 
     *  
     * 
     * @param probability
     * @param mean
     * @param sd
     * @return the inverse of the probability
     */
    template<class T>
    const T qnorm(const T &probability, const T &mean, const T &sd) {
//          return (T(22619537.0) * ::fabs(sd) * atl::InverseErf(T(2.0) * probability-T(1.0)) + T(15994428.0) * mean) / T(15994428.0);
        double guess = mean;
        double high = (100);
        double low = (-100);
        double x = guess;
        double xNew = guess;
        double error, pdf, dx = (1.0);
        int i = 0;
     
        while (std::fabs(dx) > T(DEFAULT_EPSILON) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pnorm<double>(x, mean, sd) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dnorm<T>(x, mean, sd);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }
       
        T e2 =pnorm<T>(x, mean, sd) - probability;
        T pdf2 = dnorm<T>(x, mean, sd);
        T xx = x;
        return (xx-(e2/pdf2));

    }

    /**
     * \ingroup NormalDistribution
     * 
     * @brief Random deviate
     * 
     * Returns a random number within a normal distribution defined by the
     * mean and standard deviation.
     * 
     * <br><b>Source</b> http://en.wikipedia.org/wiki/Normal_distribution<br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013
     * 
     * @param mean
     * @param standard_deviation
     * @return random number x.
     */
    template<class T>
    const T rnorm(const T &mean, const T &standard_deviation) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qnorm<T> (x, mean, standard_deviation);
        return ret;
    }

    /**
     * \ingroup NormalDistribution
     * @brief Derived from the abstract DistributionBase class.
     * Class NormalDistribution.
     */
    template<class T>
    class NormalDistribution : public DistributionBase<T> {
    public:

        NormalDistribution(T mean = T(0.0), T standard_deviation = T(1.0)) :
        mean(mean),
        standard_deviation(standard_deviation) {

        }

        /**
         * Return the mean.
         * 
         * @return 
         */
        T GetMean() const {
            return mean;
        }

        /**
         * Set the mean.
         * 
         * @param mean
         */
        void SetMean(T mean) {
            this->mean = mean;
        }

        /**
         * Return the standard deviation.
         * @return 
         */
        T GetStandardDeviation() const {
            return standard_deviation;
        }

        /**
         * Set the standard deviation.
         * 
         * @param standard_deviation
         */
        void SetStandardDeviation(T standard_deviation) {
            this->standard_deviation = standard_deviation;
        }

        /**
         * The probability density function for this normal distribution.
         * 
         * see T dnorm(const T &x, const T &mean, const T &sd)
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return dnorm<T > (x, this->GetMean(), this->GetStandardDeviation());
        }

        /**
         * The cumulative distribution function for this normal distribution.
         * 
         * see T pnorm(const T &x, const T &mean, const T &sd)
         * 
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return pnorm<T > (x, this->GetMean(), this->GetStandardDeviation());
        }

        /**
         * The inverse of the cumulative normal distribution function for 
         * this normal distribution.
         * 
         * see T qnorm(const T &probability, const T &mean, const T &sd)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return qnorm<T > (probability, this->GetMean(), this->GetStandardDeviation());
        }

        /**
         * Random number.
         * 
         * see T rnorm(const T &mean, const T &sd)
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return rnorm<T > (this->GetMean(), this->GetStandardDeviation());
        }


    private:
        T mean;
        T standard_deviation;

    };


    /**
     * @defgroup BetaDistribution Beta Distribution
     * 
     *@ingroup Distributions  
     *
     */

    /**
     * \ingroup BetaDistribution
     * 
     * @brief The probability density function for a beta distribution.
     * 
     * \f$
     * f(x) = exp(-BetaLn(\alpha,\beta)+(\alpha-1)ln(x)+(\beta-1)ln(1-x))
     * \f$
     * <br><br><b>Where,</b><br>
     * \f$\alpha\f$ is a shape parameter.<br>
     * \f$\beta\f$ is a shape parameter.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/f/f3/Beta_distribution_pdf.svg
     * 
     * <br><b>Source:</b><br> http://en.wikipedia.org/wiki/Beta_distribution <br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013
     * 
     * @param x
     * @param a
     * @param b
     * @return relative likelihood for this random variable to have value x.
     */
    template<class T>
    const T dbeta(const T &x, const T &a, const T &b) {
        if (x <= T(0.0) || x >= T(1.0)) {
            return 0.0;
        } else {
            T ret = T(-1) * atl::BetaLn<T > (a, b) + (a - T(1.0)) * std::log(x)+(b - T(1.0)) * std::log(T(1.0) - x);
            return (T) std::exp(ret);
        }
    }

    /**
     * \ingroup BetaDistribution
     * 
     * @brief Beta cumulative distribution function.
     * 
     * \f$
     * f(x) = I_x(\alpha,\beta)
     * \f$
     * <br><br><b>Where,</b><br>
     * \f$\alpha\f$ is a shape parameter.<br>
     * \f$\beta\f$ is a shape parameter.<br>
     * \f$I_x\f$ is the incomplete beta function.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/1/11/Beta_distribution_cdf.svg
     * 
     * <br><b>Source</b> http://en.wikipedia.org/wiki/Beta_distribution<br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013
     * 
     * @param x
     * @param a
     * @param b
     * @return the probability that a stochastic variable x is less then X.
     */
    template<class T>
    const T pbeta(const T &x, const T &a, const T &b) {
        if (x <= T(0)) {

            return T(0);
        } else if (x >= T(1)) {

            return T(1);
        } else {
            T ret = atl::IncompleteBeta<T > (x, a, b);
            return ret;
        }

    }

    /**
     * \ingroup BetaDistribution
     * 
     * @brief Inverse of the beta cumulative distribution function.
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013 
     * 
     * @param probability
     * @param a
     * @param b
     * @return 
     */
    template<class T>
    const T qbeta(const T &probability, const T &a, const T &b) {
        T result = T(0);
        T guess = T(0.1);
        T high = T(100);
        T low = T(-100);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;
        while (std::fabs(dx) > T(DEFAULT_EPSILON) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pbeta<T> (x, a, b) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dbeta<T> (x, a, b);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup BetaDistribution
     * 
     * Random deviate from a beta distribution.
     * @param a
     * @param b
     * @return random number x.
     */
    template<class T>
    const T rbeta(const T &a, const T &b) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qbeta<T> (x, a, b);
        return ret;

    }

    /**
     * \ingroup BetaDistribution
     * @brief Derived from the abstract DistributionBase class.
     * Class BetaDistribution.
     */
    template<class T>
    class BetaDistribution : public DistributionBase<T> {
    public:

        BetaDistribution(T alpha = T(0), T beta = T(0)) :
        alpha_m(alpha),
        beta_m(beta) {

        }

        /**
         * Return the alpha shape parameter.
         * @return 
         */
        T GetAlpha() const {
            return alpha_m;
        }

        /**
         * Set the alpha shape parameter.
         *  
         */
        void SetAlpha(T alpha) {
            this->alpha_m = alpha;
        }

        /**
         * Return the beta shape parameter.
         * @return 
         */
        T GetBeta() const {
            return beta_m;
        }

        /**
         * Set the alpha shape parameter.
         * 
         */
        void SetBeta(T beta) {
            this->beta_m = beta;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dbeta(const T &x, const T &a, const T &b)
         * 
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dbeta(x, this->GetAlpha(), this->GetBeta());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T qbeta(const T &x, const T &a, const T &b)
         * 
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pbeta(x, this->GetAlpha(), this->GetBeta());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * see T qbeta(const T &robability, const T &a, const T &b)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qbeta(probability, this->GetAlpha(), this->GetBeta());
        }

        /**
         * Random number.
         * 
         * see T rbeta(const T &a, const T &b)
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rbeta(this->GetAlpha(), this->GetBeta());
        }


    private:
        T alpha_m;
        T beta_m;
    };



    /**
     *@defgroup BinomialDistribution Binomial Distribution
     * 
     *@ingroup Distributions  
     *
     */

    /**
     * \ingroup BinomialDistribution
     * 
     * @brief The probability mass function for a binomial distribution.
     * 
     * An approximation to the binomial probability density function is used to accommodate
     * large n values using:<br><br>
     * 
     * \f$
     * 
     * p(x) = exp(ln(n)! -ln(x)! -  ln(n-x)! + xln(p)+(n-x)ln(1-p))
     * 
     * \f$
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/7/75/Binomial_distribution_pmf.svg
     * 
     * <br><b>Source:</b> http://mathworld.wolfram.com/BinomialDistribution.html<br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013 
     * 
     * @param x
     * @param n
     * @param p
     * @return  relative likelihood for this random variable to have value x.
     */
    template<class T>
    const T dbinom(const T &x, const T &n, const T &p) {

        if (x < T(0.0) || x > n) {
            return T(0);
        }

        //use log space.
        T lbicof = LogFactorial(n) - LogFactorial(x) - LogFactorial(n - x);
        T ret = lbicof + (x) * std::log(p) + (n - x) * std::log(T(1.0) - p);

        // }
        return std::exp(ret);
    }

    /**
     * @ingroup BinomialDistribution
     * 
     * @brief The binomial cumulative distribution function.
     * 
     * An approximation to the binomial cumulative distribution function is used to accommodate
     * large n values using:<br><br>
     * \f$
     * f(x) = \sum\limits_{i=0}^x exp(ln(n)! -ln(i)! -  ln(n-i)! + iln(p)+(n-i)ln(1-p))
     * \f$
     * <br><b>Where,</b><br>
     * \f$x\f$ is the number of successes.<br>
     * \f$n\f$ is the number of trials.<br>
     * \f$p\f$ is the probability of success in each trial.<br>
     *
     * @image html http://upload.wikimedia.org/wikipedia/commons/5/55/Binomial_distribution_cdf.svg
     * 
     * <br><b>Source:</b> http://en.wikipedia.org/wiki/Binomial_distribution<br>
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013 
     * 
     * @param x
     * @param n
     * @param p
     * @return the probability that a stochastic variable x is less then X.
     */
    template<class T>
    const T pbinom(const T &x, const T &n, const T &p) {
        if (n <= x) {
            return T(1);
        }

        T sum = T(0);
        for (T i = 0; i < x + 1; i++) {
            sum += dbinom(i, n, p);
        }
        return sum;

    }

    /**
     * \ingroup BinomialDistribution
     * 
     * @brief Inverse of the binomial cumulative distribution function.
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013 
     * 
     * 
     * @param probability
     * @param n
     * @param p
     * @return 
     */
    template<class T>
    const size_t qbinom(const T &probability, const T &n, const T &p) {

        T result = T(0);
        T guess = probability*n;
        T high = T(100);
        T low = T(-100);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;
        T std = sqrt((n * p * (1 - p)));
        T np = n*p;

        do {//need to work this out beter.



            error = pnorm<T > (x, np, std) - probability;

            if (std::fabs(error) <= T(DEFAULT_EPSILON)) {
                break;
            }


            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dnorm<T > (x, np, std);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        } while (std::fabs(dx) > T(DEFAULT_EPSILON) && i++ < MAX_ITERATIONS);

        result = x;

        return result;


    }

    /**
     *  \ingroup BinomialDistribution
     * 
     * @brief Returns a random number within a binomial distribution.
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     * @param n
     * @param p
     * @return random number x.
     */
    template<class T>
    const T rbinom(const T &n, const T &p) {
        T x = ((T) rand() / (T) RAND_MAX);

        T ret = qbinom<T> (x, n, p);
        return ret;
    }

    /**
     * @ingroup BinomialDistribution
     * @brief Derived from the abstract DistributionBase class.
     * 
     * Class BinomialDistribution
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     *
     */
    template<class T>
    class BinomialDistribution : public DistributionBase<T> {
    public:

        BinomialDistribution(T n = T(10), T p = T(.5)) :
        n_m(n),
        p_m(p) {

        }

        T GetNumberOfTrials() const {
            return n_m;
        }

        void SetNumberOfTrials(T n) {
            this->n_m = n;
        }

        T GetProbabilityOfSuccess() const {
            return p_m;
        }

        void SetProbabilityOfSuccess(T p) {
            this->p_m = p;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dbinom(const T &x, const T &n, const T &p)
         
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dbinom<T > (x, this->GetNumberOfTrials(), this->GetProbabilityOfSuccess());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pbinom(const T &x, const T &n, const T &p) 
         *
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pbinom<T > (x, this->GetNumberOfTrials(), this->GetProbabilityOfSuccess());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * see T qbinom(const T &probability, const T &n, const T &p) 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qbinom<T > (probability, this->GetNumberOfTrials(), this->GetProbabilityOfSuccess());
        }

        /**
         * Random number.
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * 
         * see T rbeta(const T &a, const T &b)
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rbinom<T > (this->GetNumberOfTrials(), this->GetProbabilityOfSuccess());
        }

    private:
        T n_m;
        T p_m;

    };


    /**
     * @defgroup Cauchy Cauchy Distribution
     * 
     * @ingroup Distributions
     */

    /**
     * @ingroup Cauchy
     * 
     * @brief The cumulative density function for a Cauchy distribution.
     * 
     * \f$
     * 
     * f(x) = \frac{1}{\pi} arctan(\frac{x-x_0}{\gamma}) + \frac{1}{2}
     * \f$
     * 
     * <br><b>Where,</b><br><br>
     * \f$x_0\f$ is location.<br>
     * \f$\gamma\f$ is scale.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/5/5b/Cauchy_cdf.svg
     * 
     * <br><b>Source:</b> http://en.wikipedia.org/wiki/Cauchy_distribution<br>
     *
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     *  
     * @param x
     * @param y
     * @param z
     * @return 
     */
    template<class T>
    const T pcauchy(const T &x, const T &location, const T &scale) {
        return (T(1.0) / T(M_PI))*std::atan((x - location) / scale) + T(.5);
    }

    /**
     * @ingroup Cauchy
     * 
     * @brief The probability density function.
     * 
     * \f$
     * f(x) = \frac{1}{\pi\gamma[1 + (\frac{x-x_0}{\gamma})^2]}
     * \f$
     * 
     * <br><b>Where,</b><br><br>
     * \f$x_0\f$ is location.<br>
     * \f$\gamma\f$ is scale.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/5/5b/Cauchy_cdf.svg
     * 
     * <br><b>Source:</b> http://en.wikipedia.org/wiki/Cauchy_distribution<br>
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param x
     * @param y
     * @param z
     * @return 
     */
    template<class T>
    const T dcauchy(const T &x, const T &location, const T &scale) {
        return T(1.0) / (T(M_PI) * scale * (T(1) +(x - location)*(x - location)));
    }

    /**
     * @ingroup Cauchy
     * 
     * @brief The inverse cumulative distribution function for a Cauchy distribution.
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param probability
     * @param location
     * @param scale
     * @return 
     */
    template<class T>
    const T qcauchy(const T &probability, const T &location, const T &scale) {

        if (probability == T(1.0)) {
            return T(std::numeric_limits<T>::infinity());
        }

        return location + scale / std::tan(T(M_PI) * probability);

    }

    /**
     * @ingroup Cauchy 
     * 
     * @brief Returns a random number within a Cauchy distribution.
     * 
     * @param location
     * @param scale
     * @return a random number x.
     */
    template<class T>
    const T rcauchy(const T &location, const T &scale) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = T(qcauchy<T> (x, location, scale));
        return ret;
    }

    /**
     *@ingroup Cauchy
     * 
     * @brief Derived from the abstract DistributionBase class.
     * 
     * Class CauchyDistribution
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     */
    template<class T>
    class CauchyDistribution : public DistributionBase<T> {
    public:

        CauchyDistribution(T location = T(0), T scale = T(.5)) :
        location_m(location),
        scale_m(scale) {
        }

        T GetLocation() const {
            return this->location_m;
        }

        void SetLocation(T location) {
            this->location_m = location;
        }

        T GetScale()const {
            return this->scale_m;
        }

        void SetScale(T scale) {
            this->scale_m = scale;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dcauchy(const T &x, const T &location, const T &scale)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dcauchy<T > (x, this->GetLocation(), this->GetScale());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pcauchy(const T &x, const T &location, const T &scale)
         *
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pcauchy<T > (x, this->GetLocation(), this->GetScale());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * see T qbinom(const T &probability, const T &n, const T &p) 
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qcauchy<T > (probability, this->GetLocation(), this->GetScale());
        }

        /**
         * Random number.
         * 
         * see T rbinom(const T &n, const T &p) 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rcauchy<T > (this->GetLocation(), this->GetScale());
        }


    private:
        T location_m;
        T scale_m;

    };


    /**
     * @defgroup Poisson Poisson Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Poisson
     * 
     * @brief The probability mass function.
     * 
     * \f$
     * f(x) = \frac{\lambda^x}{x!}exp(-\lambda)
     * \f$
     * 
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/1/16/Poisson_pmf.svg
     * 
     * <br><b>Source: </b> http://en.wikipedia.org/wiki/Poisson_distribution<br>
     * 
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     * 
     * @param x
     * @param lambda
     * @return 
     */
    template<class T>
    const T dpois(const T &x, const T &lambda) {
        return (exp (x * std::log(lambda) - lambda - atl::GammaLn(x + T(1))));
    }

    /**
     * @ingroup Poisson
     * 
     * @brief The cumulative density function for a Cauchy distribution.
     * 
     * \f$
     * f(x) = exp(-\lambda)\sum\limits_{i=0}^x\frac{\lambda^i}{i!}
     * \f$
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/7/7c/Poisson_cdf.svg
     * 
     * <br><b>Source: </b> http://en.wikipedia.org/wiki/Poisson_distribution<br>
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param x
     * @param lambda
     * @return 
     */
    template<class T>
    const T ppois(const T &x, const T &lambda) {
        T sum = T(0.0);
        for (T i = T(0.0); i < std::floor(x) + T(1); i++) {
            sum += (std::pow(lambda, i) / T(atl::Factorial(i)));
        }

        return exp (T(-1) * lambda) * sum;
    }

    /**
     * @ingroup Poisson
     * 
     * @brief The inverse cumulative distribution function for a Poisson distribution.
     *  
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param probability
     * @param lambda
     * @return 
     */
    template<class T>
    const int qpois(const T &probability, const T &lambda) {
        T result = T(0);
        T guess = 4;
        T high = T(100 * lambda);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = ppois<T > (x, lambda) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dpois<T > (x, lambda);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup Poisson
     * 
     * @brief Returns a random number within a Poisson distribution.
     *  
     *
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     *  
     * @param lambda
     * @return a random value x.
     */
    template<class T>
    const T rpois(const T &lambda) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qpois<T> (x, lambda);
        return ret;
    }

    /**
     * @ingroup Poisson
     * 
     * @brief Derived from the abstract DistributionBase class.
     * 
     * Class PoissonDistribution
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     */
    template<class T>
    class PoissonDistribution : public DistributionBase<T> {
    public:

        PoissonDistribution(T lambda = T(1)) :
        lambda_m(lambda) {

        }

        T GetLambda()const {
            return this->lambda_m;
        }

        void SetLambda(T lambda) {
            this->lambda_m = lambda;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dpois(const T &x, const T &lambda)
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dpois<T > (x, this->GetLambda());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T ppois(const T &x, const T &lambda)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::ppois<T > (x, this->GetLambda());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * see T qpois(const T &probability, const T &lambda)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qpois<T > (probability, this->GetLambda());
        }

        /**
         * Random number.
         * 
         * see T rpois(const T &lambda)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rpois<T > (this->GetLambda());
        }


    private:
        T lambda_m;

    };

    /**
     * @defgroup GammaDist Gamma Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup GammaDist
     * 
     * @brief The cumulative density function for a gamma distribution.
     * 
     * \f$
     * f(x) = \frac{1}{\Gamma(k)}\gamma(k,\frac{x}{\theta})
     * \f$
     * <br><b>Where,<br></b>
     * \f$k\f$ is the shape.<br>
     * \f$\theta\f$ is the scale.<br>
     * \f$\Gamma\f$ is the gamma function.<br>
     * \f$\gamma\f$ is the incomplete gamma function.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/8/8d/Gamma_distribution_cdf.svg/500px-Gamma_distribution_cdf.svg.png
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param x
     * @param shape
     * @param scale
     * @return 
     */
    template<class T>
    const T pgamma(const T &x, const T &shape, const T &scale) {

        if (shape > T(0) && scale > T(0)) {
            // return   (T(1.0)/atl::Gamma(shape)) *atl::IncompleteGamma(shape,x/scale);

            T sum = T(0);

            for (int i = 0; i < shape; i++) {
                sum += T(1.0) / atl::Factorial(i) * std::pow((x / scale), (T(i))) * std::exp(T(-1)*(x / scale));

            }

            return T(1) - sum;

        }

        return T(0);
    }

    /**
     *    
     * @ingroup GammaDist
     * 
     * @brief The probability density function.
     * 
     * \f$
     * f(x) = \frac{1}{\Gamma(k)\theta^k}x^{k-1}e^{-\frac{x}{\theta}}
     * \f$
     * 
     * <br><b>Where,<br></b>
     * \f$k\f$ is the shape.<br>
     * \f$\theta\f$ is the scale.<br>
     * \f$\Gamma\f$ is the gamma function.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/e/e6/Gamma_distribution_pdf.svg/500px-Gamma_distribution_pdf.svg.png
     * 
     * <br><b>Source:</b> http://en.wikipedia.org/wiki/Gamma_distribution<br>
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param x
     * @param shape
     * @param scale
     * @return the probability that a stochastic variable x has the value X.
     */
    template<class T>
    const T dgamma(const T &x, const T &shape, const T &scale) {
        return (T(1.0) / std::pow(scale, shape))*(T(1.0) / atl::Gamma<T > (shape)) * std::pow(x, shape - T(1)) * std::exp(T(-1) * x / scale);
    }

    /**
     * @ingroup GammaDist
     * 
     * @brief The inverse cumulative distribution function for a gamma distribution.
     *  
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param probability
     * @param shape
     * @param scale
     * @return 
     */
    template<class T>
    const T qgamma(const T &probability, const T &shape, const T &scale) {


        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pgamma<T > (x, shape, scale) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dgamma<T > (x, shape, scale);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }



        return x;
    }

    /**
     * @ingroup GammaDist
     * 
     * @brief Returns a random number within a gamma distribution.
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     * 
     * @param shape
     * @param scale
     * @return a random value x. 
     */
    template<class T>
    const T rgamma(const T &shape, const T &scale) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qgamma<T> (x, shape, scale);
        return ret;
    }

    /**
     * @ingroup GammaDist
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class GammaDistribution
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     * 
     */
    template<class T>
    class GammaDistribution : public DistributionBase<T> {
    public:

        GammaDistribution(T shape = T(1.0), T scale = T(2.0)) :
        shape_m(shape),
        scale_m(scale) {

        }

        T GetShape() const {
            return this->shape_m;
        }

        void SetShape(T shape) {
            this->shape_m = shape;
        }

        T GetScale() const {
            return this->scale_m;
        }

        void SetScale(T scale) {
            this->scale_m = scale;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dgamma(const T &x, const T &shape, const T &scale)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dgamma<T > (x, this->GetShape(), this->GetScale());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pgamma(const T &x, const T &shape, const T &scale)
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pgamma<T > (x, this->GetShape(), this->GetScale());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qgamma(const T &probability, const T &shape, const T &scale)
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qgamma<T > (probability, this->GetShape(), this->GetScale());
        }

        /**
         * Random number.
         * 
         * see T rgamma(const T &shape, const T &scale)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rgamma<T > (this->GetShape(), this->GetScale());
        }



    private:
        T shape_m;
        T scale_m;
    };

    /**
     * @defgroup Chi-Square Chi-Square Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Chi-Square
     * 
     * @brief The cumulative distribution function for a chi-squared distribution.
     * 
     * \f$
     * f(x) = \frac{1}{\Gamma(\frac{k}{2})}\gamma(\frac{k}{2},\frac{x}{2})
     * \f$
     * 
     * <br><b>Where,<br></b>
     * \f$k\f$ is "degrees of freedom"
     * \f$\Gamma\f$ is the Gamma function. 
     * \f$\gamma\f$ is the lower incomplete Gamma function. 
     *  
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/3/35/Chi-square_pdf.svg/500px-Chi-square_pdf.svg.png
     * 
     * <br><b>Source: </b> http://en.wikipedia.org/wiki/Chi-squared_distribution<br>
     * 
     *
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     *   
     * @param x
     * @param df
     * @return 
     */
    template<class T>
    const T pchisq(const T &x, const T &df) {

        return (T(1.0) / (atl::Gamma(df / T(2))))*atl::LowerIncompleteGamma(df / T(2.0), x / T(2.0));
        //return T(1.0) - (T(1.0) / atl::Gamma(df / T(2.0))) * atl::LowerIncompleteGamma(df / T(2.0), x / T(2.0)); alternative version
    }

    /**
     * @ingroup Chi-Square
     * 
     * @brief The probability density function for a chi-squared distribution.
     * 
     * \f$
     * f(x) = \frac{1}{2\frac{k}{2}\Gamma(\frac{k}{2})}x^{\frac{k}{2}-1}e^{-\frac{x}{2}}
     * \f$
     *<br><b>Where,<br></b>
     * 
     * \f$k\f$ is "degrees of freedom"
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/0/01/Chi-square_cdf.svg/500px-Chi-square_cdf.svg.png
     * 
     * <br><b>Source: </b> http://en.wikipedia.org/wiki/Chi-squared_distribution<br>
     * 
     * @param x
     * @param df
     * @return 
     */
    template<class T>
    const T dchisq(const T &x, const T &df) {
        if (x >= T(0)) {
            //            return (T(1.0) / (std::pow(T(2.0), df / T(2.0)) * atl::Gamma<T > (df / T(2.0))))*
            //                    std::pow(x, (df / T(2.0)) - T(1)) * std::exp(T(-1) * x / T(2.0)); //alternative version
            return ((T(1) / std::pow(T(2), (df / T(2.0))) * atl::Gamma(df / T(2.0))) * std::pow(x, ((df / T(2.0)) - T(1.0))) * std::exp(T(-1.0) * x / T(2)));
        }

        return T(0);
    }

    /**
     * @ingroup Chi-Square
     * 
     * @brief The inverse to the cumulative distribution function for a chi-squared distribution.
     * 
     * @param probability
     * @param df
     * @return 
     */
    template<class T>
    const T qchisq(const T &probability, const T &df) {

        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pchisq<T > (x, df) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dchisq<T > (x, df);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }


        return x;
    }

    /**
     * @ingroup Chi-Square
     
     * @brief Returns a random number within a chi-squared distribution.
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     * 
     * @param df
     * @return 
     */
    template<class T>
    const T rchisq(const T &df) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qchisq<T> (x, df);
        return ret;
    }

    /**
     *@ingroup Chi-Square
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class ChiSquaredDistribution
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     */
    template<class T>
    class ChiSquaredDistribution : public DistributionBase<T> {
    public:

        ChiSquaredDistribution(T degrees_of_freedom = T(1)) :
        df_m(degrees_of_freedom) {

        }

        T GetDegreesOfFreedom()const {
            return this->df_m;
        }

        void SetDegreesOfFreedom(T degrees_of_freedom) {
            this->df_m = degrees_of_freedom;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dchisq(const T &x, const T &df)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dchisq(x, this->GetDegreesOfFreedom());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pchisq(const T &x, const T &df)
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pchisq(x, this->GetDegreesOfFreedom());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qchisq(const T &probability, const T &df)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qchisq(probability, this->GetDegreesOfFreedom());
        }

        /**
         * Random number.
         * 
         * see T rchisq(const T &df)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rchisq(this->GetDegreesOfFreedom());
        }



    private:
        T df_m;

    };

    /**
     * @defgroup Exponential Exponential Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Exponential
     * 
     * @brief The cumulative distribution function for a exponential distribution.
     * 
     * \f$
     * f(x) = 1- e^{-\gamma x}
     * \f$
     * <br><br><b>Where,</b><br>
     * \f$\gamma\f$ is the rate (\f$\gamma\f$ > 0).<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Exponential_cdf.svg/500px-Exponential_cdf.svg.png
     * 
     * <br><b>Source: </b> http://en.wikipedia.org/wiki/Exponential_distribution<br>
     * 
     *
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013   
     * @param x
     * @param rate
     * @return 
     */
    template<class T>
    const T pexp(const T &x, const T &rate) {
        if (rate >= T(0)) {
            return T(1) - std::exp(T(-1) * rate * x);
        }

        return T(0);

    }

    /**
     * @ingroup Exponential
     * 
     * @brief The probability density function for a exponential distribution.
     * 
     * \f$
     * f(x) = \gamma e^{-\gamma x}
     * \f$
     *<br><br><b>Where,</b><br>
     * \f$\gamma\f$ is the rate (\f$\gamma\f$ > 0).<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/e/ec/Exponential_pdf.svg/500px-Exponential_pdf.svg.png
     
     * 
     * <br><b>Source: </b> http://en.wikipedia.org/wiki/Exponential_distribution<br>
     * 
     *
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *  
     * @param x
     * @param rate
     * @return 
     */
    template<class T>
    const T dexp(const T &x, const T &rate) {
        if (rate >= T(0)) {
            return rate * std::exp(T(-1) * rate * x);
        }

        return T(0);
    }

    /**
     * @ingroup Exponential
     * 
     * @The Inverse cumulative distribution function for the exponential distribution. 
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     * @param probability
     * @param rate
     * 
     * @return 
     */
    template<class T>
    const T qexp(const T &probability, const T &rate) {
        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pexp<T > (x, rate) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dexp<T> (x, rate);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup Exponential
     * 
     * @brief Returns a random number within a exponential distribution.
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * @param rate
     * @return 
     */
    template<class T>
    const T rexp(const T &rate) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qexp<T> (x, rate);
        return ret;
    }

    /**
     *@ingroup Exponential
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class ExponentialDistribution
     * 
     *    
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *
     */
    template<class T>
    class ExponentialDistribution : public DistributionBase<T> {
    public:

        ExponentialDistribution(T rate = T(1)) :
        rate_m(rate) {

        }

        T GetRate()const {
            return this->rate_m;
        }

        void SetRate(T rate) {
            this->rate_m = rate;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dexp(const T &x, const T &rate)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dexp(x, this->GetRate());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pexp(const T &x, const T &rate)
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pexp(x, this->GetRate());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qexp(const T &probability, const T &rate)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qexp(probability, this->GetRate());
        }

        /**
         * Random number.
         * 
         * see T rexp(const T &rate)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rexp(this->GetRate());
        }



    private:
        T rate_m;

    };

    /**
     * @defgroup FDist F Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup FDist
     * 
     * @brief The cumulative distribution function for a F distribution.
     * 
     * \f$
     * f(x) = I_{\frac{d_1 x}{d_1 x +d_2}}(\frac{d_1}{2},\frac{d_2}{2})
     * \f$
     * 
     * <b><br>Where,</b><br><br>
     * \f$d_1, d_2\f$ > 0 degrees of freedom.<br>
     * \f$I\f$ is the regularized incomplete beta function.<br><br>
     *  
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/d/df/F_distributionCDF.png/320px-F_distributionCDF.png
     * 
     * 
     * <b>Source:</b> http://en.wikipedia.org/wiki/F-distribution<br>
     * 
     * @author Matthew Supernaw
     * 
     * @param x
     * @param d1
     * @param d2
     * @return the probability that a stochastic variable x is less then X.
     */
    template<class T>
    const T pf(const T &x, const T &d1, const T &d2) {
        return atl::IncompleteReqularizedBeta((d1 * x) / (d1 * x + d2), d1 / T(2), d2 / T(2));
    }

    /**
     * @ingroup FDist
     * 
     * @brief The probability density function for a F distribution.
     * 
     * \f$
     * f(x) = \frac{\sqrt{\frac{(d_1 x)^{d_1} d_2^{d_2}}{(d_1 x +d_2)^{d_1+d_2}}}}{x B(\frac{d_1}{2},\frac{d_2}{2})}
     * \f$
     * 
     * <b><br>Where,</b><br>
     * \f$d_1, d_2\f$ > 0 degrees of freedom.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/f/f7/F_distributionPDF.png/320px-F_distributionPDF.png
     * 
     * <b>Source:</b> http://en.wikipedia.org/wiki/F-distribution<br>
     *
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013  
     *  
     * @param x
     * @param d1
     * @param d2
     * @return 
     */
    template<class T>
    const T df(const T &x, const T &d1, const T &d2) {
        if (x >= T(0) && d1 > T(0) && d2 > T(0)) {
            return std::sqrt((std::pow(x*d1, d1)*(d2 * d2)) / std::pow(x * d1 + d2, d1 + d2)) / (x * atl::Beta<T > (d1 / T(2.0), d2 / T(2.0)));
        }
        return T(0);
    }

    /**
     * @ingroup FDist
     * 
     * @brief The inverse to the cumulative distribution function.
     * 
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     *  
     * @param probability
     * @param d1
     * @param d2
     * @return 
     */
    template<class T>
    const T qf(const T &probability, const T &d1, const T &d2) {

        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pf<T > (x, d1, d2) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = df<T > (x, d1, d2);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;

    }

    /**
     * @ingroup FDist
     * 
     * @brief Returns a random number within this distribution.
     * 
     * @author Matthew Supernaw
     * @date June 6, 2013  
     * 
     * @param d1
     * @param d2
     * @return 
     */
    template<class T>
    const T rf(const T &d1, const T &d2) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qf<T>(x, d1, d2);
        return ret;
    }

    /**
     * @ingroup FDist
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class FDistribution
     */
    template<class T>
    class FDistribution : public DistributionBase<T> {
    public:

        FDistribution(T degree_of_freedom1 = T(1), T degree_of_freedom2 = T(1)) :
        d1_m(degree_of_freedom1),
        d2_m(degree_of_freedom2) {

        }

        T GetDegreesOfFreedom1()const {
            return this->d1_m;
        }

        void SetDegreesOfFreedom1(T degrees_of_freedom) {
            this->d1_m = degrees_of_freedom;
        }

        T GetDegreesOfFreedom2()const {
            return this->d2_m;
        }

        void SetDegreesOfFreedom2(T degrees_of_freedom) {
            this->d2_m = degrees_of_freedom;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T df(const T &x, const T &d1,const T &d2)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::df(x, this->GetDegreesOfFreedom1(), this->GetDegreesOfFreedom2());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pf(const T &x, const T &d1,const T &d2)
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pf(x, this->GetDegreesOfFreedom1(), this->GetDegreesOfFreedom2());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qf(const T &x, const T &d1,const T &d2)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qf(probability, this->GetDegreesOfFreedom1(), this->GetDegreesOfFreedom2());
        }

        /**
         * Random number.
         * 
         * see T rf(const T &x, const T &d1,const T &d2)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rf(this->GetDegreesOfFreedom1(), this->GetDegreesOfFreedom2());
        }


    private:

        T d1_m;
        T d2_m;
    };


    /**
     * @defgroup Geometric Geometric Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Geometric
     * 
     * @brief The cumulative distribution function for the geometric distribution.
     * 
     * \f$
     * f(x) = 1-(1-p)^x p
     * \f$
     * <br><br><b>Where,<br></b>
     * \f$0 < p \le 1\f$  success probability
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/6/6f/Geometric_cdf.svg/500px-Geometric_cdf.svg.png
     * 
     * <br><b>Source:</b> http://en.wikipedia.org/wiki/Geometric_distribution<br>
     * 
     * 
     * @param x
     * @param p 
     * @return 
     */
    template<class T>
    const T pgeom(const T &x, const T &p) {
        if (x >= T(0) && p > T(0) && p < T(1)) {
            return T(1) - std::pow((T(1) - p), x + T(1));
        }
        return T_NAN;
    }

    /**
     * @ingroup Geometric
     * 
     * @brief The probability mass function for the geometric distribution.
     * 
     * \f$
     * f(x) = (1-p)^{x-1}p
     * \f$
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/4/4b/Geometric_pmf.svg/500px-Geometric_pmf.svg.png
     * 
     * <br><b>Source:</b> http://en.wikipedia.org/wiki/Geometric_distribution<br>
     *  
     * @param x
     * @param p
     * @return 
     */
    template<class T>
    const T dgeom(const T &x, const T &p) {
        if (x >= T(0) && p > T(0) && p < T(1)) {
            return p * std::pow((T(1) - p), x);
        }
        return T_NAN;
    }

    /**
     * @ingroup Geometric
     * 
     * @brief The inverse cumulative distribution function.
     * 
     * 
     * @param probability
     * @param p
     * @return 
     */
    template<class T>
    const T qgeom(const T &probability, const T &p) {

        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pgeom<T > (x, p) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dgeom<T > (x, p);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup Geometric
     * 
     * @brief Returns a random value in this distribution.
     * 
     * @param p
     * @return 
     */
    template<class T>
    const T rgeom(const T &p) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qgeom<T>(x, p);
        return ret;
    }

    /**
     * @ingroup Geometric
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class GeometricDistribution
     * 
     */
    template<class T>
    class GeometricDistribution : public DistributionBase<T> {
    public:

        GeometricDistribution(T probability = T(.2)) :
        probability_m(probability) {

        }

        T GetProbability()const {
            return this->probability_m;
        }

        void SetProbability(T probability) {
            this->probability_m = probability;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dgeom(const T &x, const T &p)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dgeom(x, this->GetProbability());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T pgeom(const T &x, const T &p)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::pgeom(x, this->GetProbability());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qgeom(const T &x, const T &p)
         * 
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qgeom(probability, this->GetProbability());
        }

        /**
         * Random number.
         * 
         * see T rgeom(const T &p)
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rgeom(this->GetProbability());
        }


    private:
        T probability_m;
    };

    /**
     * @defgroup Hypergeometric Hypergeometric Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Hypergeometric
     * 
     * @brief The probability mass function for the hypergeometric distribution.
     * 
     * \f$
     * f(x) = \frac{{m \choose x} {n \choose k-x}}{{m+n \choose x}}
     * \f$
     * 
     * <b><br>Where,</b><br><br>
     * \f$n\f$ is the number of black balls in the urn<br>
     * \f$m\f$ is the number of white balls in the urn.<br>
     * \f$k\f$ is the number of draws.<br>
     * \f$x\f$ the number of white balls drawn without replacement from an urn which contains both black and white balls.<br>
     * 
     * 
     * @param x
     * @param m
     * @param n
     * @param k
     * @return 
     */
    template<class T>
    const T dhyper(const T &x, const T &m, const T &n, const T &k) {
        return (T(atl::Choose(m, x)) * T(atl::Choose(n, (k - x)))) / T(atl::Choose((m + n), k));
    }

    /**
     * @ingroup Hypergeometric
     * 
     * @brief The cumulative distribution function for a hypergeometric distribution.
     * 
     * 
     * \f$
     * f(x) = \sum\limits_{i=0}^x \frac{{m \choose i} {n \choose k-i}}{{m+n \choose i}}
     * \f$
     * 
     * <b><br>Where,</b><br><br>
     * \f$n\f$ is the number of black balls in the urn<br>
     * \f$m\f$ is the number of white balls in the urn.<br>
     * \f$k\f$ is the number of draws.<br>
     * \f$x\f$ the number of white balls drawn without replacement from an urn which contains both black and white balls.<br>
     * 
     * 
     * @param x
     * @param m
     * @param n
     * @param k
     * @return 
     */
    template<class T>
    const T phyper(const size_t &x, const T &m, const T &n, const T &k) {

        T sum = T(0);

        for (T i = T(0); i < x + 1; i++) {
            sum += dhyper<T > (i, m, n, k);
        }
        return sum;


    }

    /**
     * @ingroup Hypergeometric
     * 
     * @brief The inverse cumulative distribution function for a hypergeometric distribution.
     * 
     * @param probability
     * @param m
     * @param n
     * @param k
     * @return 
     */
    template<class T>
    const T qhyper(const T &probability, const T &m, const T &n, const T &k) {
        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = phyper<T > (x, m, n, k) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dhyper<T > (x, m, n, k);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = round(x);

        return result;
    }

    /**
     * @ingroup Hypergeometric
     * 
     * @brief Returns a random value within a hypergeometric distribution.
     * 
     * @param m
     * @param n
     * @param k
     * @return 
     */
    template<class T>
    const T rhyper(const T &m, const T &n, const T &k) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qhyper<T> (x, m, n, k);
        return ret;
    }

    /**
     *@ingroup Hypergeometric
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class HypergeometricDistribution
     * 
     */
    template<class T>
    class HypergeometricDistribution : public DistributionBase<T> {
    public:

        HypergeometricDistribution(T m = T(1), T n = T(1), T k = T(1)) :
        m_m(m),
        n_m(n),
        k_m(k) {

        }

        void SetM(T m) {
            this->m_m = m;
        }

        T GetM()const {
            return this->m_m;
        }

        void SetN(T n) {
            this->n_m = n;
        }

        T GetN()const {
            return this->n_m;
        }

        void SetK(T k) {
            this->k_m = k;
        }

        T GetK()const {
            return this->k_m;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dhyper(const T &x, const T &m, const T &n, const T &k) 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dhyper(x, this->GetM(), this->GetN(), this->GetK());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T phyper(const T &x, const T &m, const T &n, const T &k) 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::phyper(x, this->GetM(), this->GetN(), this->GetK());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qhyper(const T &probability, const T &m, const T &n, const T &k) 
         * 
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qhyper(probability, this->GetM(), this->GetN(), this->GetK());
        }

        /**
         * Random number.
         * 
         * see T rhyper( const T &m, const T &n, const T &k) 
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rhyper(this->GetM(), this->GetN(), this->GetK());
        }



    private:
        T m_m;
        T n_m;
        T k_m;
    };


    /**
     * @defgroup Logistic Logistic Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Logistic
     * 
     * 
     * @brief The cumulative distribution function for a logistic distribution.
     * 
     * \f$
     * f(x) = \frac{1}{1+e^{-\frac{x -\mu}{s}}}
     * \f$
     * <b><br>Where,</b><br><br>
     * \f$\mu\f$ is location.<br>
     * \f$s\f$ is scale.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/1/1a/Logistic_cdf.svg/500px-Logistic_cdf.svg.png
     * 
     * <b>Source:</b> http://en.wikipedia.org/wiki/Logistic_distribution<br><br>
     * @param x
     * @param location
     * @param scale
     * @return 
     */
    template<class T>
    const T plogis(const T &x, const T &location, const T &scale) {
        return T(1.0) / (T(1) + std::exp(T(-1)*((x - location) / scale)));
    }

    /**
     * 
     * @ingroup Logistic
     *
     * @brief The probability density function for a logistic distribution.
     * 
     * \f$
     * f(x) = \frac{e^{-\frac{x -\mu}{s}}}{s(1+ e^{-\frac{x -\mu}{s}})^2}
     * \f$
     * <b><br>Where,</b><br><br>
     * \f$\mu\f$ is location.<br>
     * \f$s\f$ is scale.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/6/66/Logisticpdfunction.svg/500px-Logisticpdfunction.svg.png
     *
     * <b>Source:</b> http://en.wikipedia.org/wiki/Logistic_distribution<br><br>
     *  
     * @param x
     * @param location
     * @param scale
     * @return 
     */
    template<class T>
    const T dlogis(const T &x, const T &location, const T &scale) {
        return std::exp(T(-1)*((x - location) / scale)) / (scale * (T(1) + std::exp(T(-1)*
                ((x - location) / scale)))*
                (T(1) + std::exp(T(-1)*((x - location) / scale))));
    }

    /**
     * @ingroup Logistic
     * 
     * @brief The inverse of the cumulative distribution function for a logistic distribution.
     * 
     * @param probability
     * @param location
     * @param scale
     * @return 
     */
    template<class T>
    const T qlogis(const T &probability, const T &location, const T &scale) {

        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = plogis<T > (x, location, scale) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dlogis<T > (x, location, scale);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup Logistic
     * 
     * @Returns a random value within a logistic distribution.
     * 
     * @param location
     * @param scale
     * @return 
     */
    template<class T>
    const T rlogis(const T &location, const T &scale) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qlogis<T>(x, location, scale);
        return ret;
    }

    /**
     *@ingroup Logistic
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * Class LogisticDistribution
     * 
     */
    template<class T>
    class LogisticDistribution : public DistributionBase<T> {
    public:

        LogisticDistribution(T scale = T(2), T location = T(1)) :
        scale_m(scale),
        location_m(location) {

        }

        T GetScale()const {
            return this->scale_m;
        }

        void SetScale(T scale) {
            this->scale_m = scale;
        }

        T GetLocation()const {
            return this->location_m;
        }

        void SetLocation(T location) {
            this->location_m = location;
        }

        /**
         * The probability density function for this distribution.
         * 
         * see T dlogis(const T &x, const T &location, const T &scale)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return atl::dlogis(x, this->GetLocation(), this->GetScale());
        }

        /**
         * The cumulative distribution function for this distribution.
         * 
         * see T plogis(const T &x, const T &location, const T &scale)
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         *  
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return atl::plogis(x, this->GetLocation(), this->GetScale());
        }

        /**
         * The inverse of the cumulative distribution function for this distribution.
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         * 
         * see T qlogis(const T &x, const T &location, const T &scale)
         * 
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return atl::qlogis(probability, this->GetLocation(), this->GetScale());
        }

        /**
         * Random number.
         * 
         * see T rlogis( const T &location, const T &scale)
         * 
         * 
         *    
         * @author Matthew Supernaw
         * @date June 6, 2013  
         *
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return atl::rlogis(this->GetLocation(), this->GetScale());
        }



    private:
        T scale_m;
        T location_m;

    };

    /**
     * @defgroup LogNormal Log Normal Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup LogNormal
     * 
     * @brief The cumulative distribution function for a log-normal distribution.
     * 
     * \f$
     * f(x) = \frac{1}{2}+\frac{1}{2}erf(\frac{ln x - \mu}{\sqrt{2\theta}})
     * \f$
     * 
     * <br><b>Where,</b><br>
     * \f$\mu\f$ is the log of the mean.<br>
     * \f$\theta\f$ is the log of the standard deviation.<br>
     * \f$erf\f$ is the error function.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/b/b8/Lognormal_distribution_CDF.svg/500px-Lognormal_distribution_CDF.svg.png
     * 
     * <b><br>Source:</b>http://en.wikipedia.org/wiki/Log-normal_distribution<br><br>
     * 
     * 
     * 
     * @param x
     * @param meanLog
     * @param sdLog
     * @return 
     */
    template<class T>
    const T plnorm(const T &x, const T &meanLog, const T &sdLog) {
        if (x >= T(0)) {
            return T(.5) + T(.5) * atl::Erf<T > ((std::log(x) - meanLog) /
                    std::sqrt(T(2.0) * sdLog));
        }

        return T(0);
    }

    /**
     * @ingroup LogNormal
     * 
     * @brief The probability distribution function for a log-normal distribution.
     * 
     * \f$
     * f(x) = \frac{1}{x\sqrt{2\pi}\theta}e^{-\frac{(ln x - \mu)^{2}}{2\theta^{2}}}
     * \f$
     * 
     * <br><b>Where,</b><br>
     * \f$\mu\f$ is the log of the mean.<br>
     * \f$\theta\f$ is the log of the standard deviation.<br>
     * 
     * @image html http://upload.wikimedia.org/wikipedia/commons/thumb/8/80/Some_log-normal_distributions.svg/500px-Some_log-normal_distributions.svg.png
     * 
     * <b><br>Source:</b>http://en.wikipedia.org/wiki/Log-normal_distribution<br><br>
     * 
     * 
     * @param x
     * @param meanLog
     * @param sdLog
     * @return 
     */
    template<class T>
    const T dlnorm(const T &x, const T &meanLog, const T &sdLog) {
        if (x > T(0)) {
            return (T(1) / (x * sdLog * std::sqrt(T(2) * T(M_PI))))*
                    std::exp(T(-1)*((std::log(x) - meanLog)*
                    (std::log(x) - meanLog)) / (T(2) * sdLog));
        }
        return T(0);
    }

    /**
     * @ingroup LogNormal
     * 
     * @brief The inverse cumulative distribution function for the log-normal distribution.
     * 
     * @param probability
     * @param meanLog
     * @param sdLog
     * @return 
     */
    template<class T>
    const T qlnorm(const T &probability, const T &meanLog, const T &sdLog) {

        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = plnorm<T > (x, meanLog, sdLog) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dlnorm<T > (x, meanLog, sdLog);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup LogNormal
     * 
     * @brief Returns a random value from a log-normal distribution.
     * 
     * @param meanLog
     * @param sdLog
     * @return 
     */
    template<class T>
    const T rlnorm(const T &meanLog, const T &sdLog) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qlnorm<T > (x, meanLog, sdLog);
        return ret;
    }

    /**
     *ingroup LogNormal
     * 
     * @brief Derived from the class DistributionBase.
     * 
     * 
     * class LogNormalDistribution.
     * 
     */
    template<class T>
    class LogNormalDistribution : public DistributionBase<T> {
    public:

        /**
         * 
         * @param mean on log scale
         * @param standard_deviation on log scale
         */
        LogNormalDistribution(T mean = T(0), T standard_deviation = T(1)) :
        mean(mean),
        standard_deviation(standard_deviation) {

        }

        /**
         * Return the mean on log scale.
         * 
         * @return 
         */
        T GetMean() const {
            return mean;
        }

        /**
         * Set the mean on log scale .
         * 
         * @param mean
         */
        void SetMean(T mean) {
            this->mean = mean;
        }

        /**
         * Return the standard deviation on log scale.
         * @return 
         */
        T GetStandardDeviation() const {
            return standard_deviation;
        }

        /**
         * Set the standard deviation on log scale.
         * 
         * @param standard_deviation
         */
        void SetStandardDeviation(T standard_deviation) {
            this->standard_deviation = standard_deviation;
        }

        /**
         * The probability density function for this normal distribution.
         * 
         * see T dnorm(const T &x, const T &mean, const T &sd)
         * 
         * @param x
         * @return the probability that a stochastic variable x has the value X.
         */
        const T Probability(const T &x) const {
            return dlnorm<T > (x, this->GetMean(), this->GetStandardDeviation());
        }

        /**
         * The cumulative distribution function for this normal distribution.
         * 
         * see T pnorm(const T &x, const T &mean, const T &sd)
         * 
         * @param x
         * @return the probability that a stochastic variable x is less then X.
         */
        const T Cumulative(const T &x) const {
            return plnorm<T > (x, this->GetMean(), this->GetStandardDeviation());
        }

        /**
         * The inverse of the cumulative normal distribution function for 
         * this normal distribution.
         * 
         * see T qnorm(const T &probability, const T &mean, const T &sd)
         * 
         * @param probability
         * @return the value X for which P(x<X).
         */
        const T Inverse(const T &probability) const {
            return qlnorm<T > (probability, this->GetMean(), this->GetStandardDeviation());
        }

        /**
         * Random number.
         * 
         * see T rnorm(const T &mean, const T &sd)
         * 
         * @return a random deviate from this distribution.
         */
        const T Random() const {
            return rlnorm<T > (this->GetMean(), this->GetStandardDeviation());
        }


    private:
        T mean;
        T standard_deviation;


    };

    /**
     * @defgroup NegativeBinomial Negative Binomial Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup NegativeBinomial
     * 
     * @param x
     * @param n
     * @param p
     * @return 
     */
    template<class T>
    const T dnbinom(const T &x, const size_t &n, const T &p) {

        size_t top = x + n - 1;
        size_t bottom = x;
        size_t combinations = atl::Factorial(top) / (atl::Factorial(bottom) * atl::Factorial(top - bottom));

        return combinations * std::pow((T(1) - p), x) * std::pow(p, T(n));

    }

    /**
     * @ingroup NegativeBinomial
     * 
     * @param x
     * @param n
     * @param p
     * @return 
     */
    template<class T>
    const T pnbinom(const T &x, const size_t &n, const T &p) {




        if (p == T(0.0)) return T(0.0);
        if (p == T(1.0)) return (x == T(0)) ? T(1.0) : T(0.0);
        T sum = T(0);

        for (int i = 0; i < x + 1; i++) {
            sum += dnbinom(T(i), n, p);

        }
        return sum;
        //        return atl::IncompleteBeta<T > (p, x, n);


    }

    /**
     * @ingroup NegativeBinomial
     * 
     * @param probability
     * @param n
     * @param p
     * @return 
     */
    template<class T>
    const T qnbinom(const T &probability, const size_t &n, const T &p) {

        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        do {//need to work this out beter.

            error = pnbinom<T > (x, n, p) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dnbinom<T > (x, n, p);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        } while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS);

        result = std::ceil(x);

        return result;
    }

    /**
     * @ingroup NegativeBinomial
     * 
     * @param n
     * @param p
     * @return 
     */
    template<class T>
    const T rnbinom(const T &n, const T &p) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qnbinom<T>(x, n, p);
        return ret;
    }

    template<class T>
    class NegativeBinomialDistribution : public DistributionBase<T> {
    };

    /**
     * @defgroup StudentT Student's T Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup StudentT
     * 
     * @param x
     * @param df
     * @return 
     */
    template<class T>
    const T pt(const T &x, const T &df) {
        T xt = df / (x * x + df);
        return T(1) - T(.5) * atl::IncompleteBeta<T > (xt, df / T(2), T(.5));

    }

    /**
     * @ingroup StudentT
     * 
     * @param x
     * @param df
     * @return 
     */
    template<class T>
    const T dt(const T &x, const T &df) {
        return (atl::Gamma<T > ((df + T(1)) / T(2)) /
                (std::sqrt(df * (T(M_PI))) * atl::Gamma<T > (df / T(2))))*
                std::pow((T(1)+(x * x) / df), (T(-1)*(df + T(1)) / T(2)));

    }

    /**
     * @ingroup StudentT
     * 
     * @param probability
     * @param df
     * @return 
     */
    template<class T>
    const T qt(const T &probability, const T &df) {

        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pt<T > (x, df) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dt<T > (x, df);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;

    }

    /**
     * @ingroup StudentT
     * 
     * @param df
     * @return 
     */
    template<class T>
    const T rt(const T &df) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qt<T>(x, df);
        return ret;
    }

    template<class T>
    class StudentsTDistribution : public DistributionBase<T> {
    };

    //Studentized Range

    //    template<class T>
    //    T ptukey(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T qtukey(const T &x, const T &y, const T &z) {
    //    }

    //    template<class T>
    //    T dtukey(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T rtukey(const T &x, const T &y, const T &z) {
    //    }

    /**
     * @defgroup UniformDist Uniform Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup UniformDist
     * 
     * @param x
     * @param min
     * @param max
     * @return 
     */
    template<class T>
    const T punif(const T &x, const T &min, const T &max) {
        if (x < min) {
            return T(0);
        } else if (x >= max) {
            return T(1);
        } else {
            return (x - min) / (max - min);
        }

    }

    /**
     * @ingroup UniformDist
     * 
     * @param x
     * @param min
     * @param max
     * @return 
     */
    template<class T>
    const T dunif(const T &x, const T &min, const T &max) {
        if (min <= x && x <= max) {
            return T(1.0) / (max - min);
        } else {
            return T(0);
        }
    }

    /**
     * 
     * @ingroup UniformDist
     * 
     * @param probability
     * @param min
     * @param max
     * @return 
     */
    template<class T>
    const T qunif(const T &probability, const T &min, const T &max) {
        return probability * (max - min);
    }

    /**
     * @ingroup UniformDist
     * 
     * @param min
     * @param max
     * @return 
     */
    template<class T>
    const T runif(const T &min, const T &max) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qunif<T> (x, min, max);
        return ret;
    }

    template<class T>
    class UniformDistribution : public DistributionBase<T> {
        T min_m;
        T max_m;

    public:

        UniformDistribution(T min = T(0), T max = T(1)) :
        min_m(min), max_m(max) {
        }

        const T Probability(const T &x) const {
            return dunif<T>(x, min_m, max_m);
        }

        const T Cumulative(const T &x) const {
            return punif<T > (x, min_m, max_m);
        }

        const T Inverse(const T &probability) const {
            return qunif<T > (probability, min_m, max_m);
        }

        const T Random() const {
            return runif<T > (min_m, max_m);
        }

        T GetMax() const {
            return max_m;
        }

        void SetMax(T max) {
            this->max_m = max;
        }

        T GetMin() const {
            return min_m;
        }

        void SetMin(T min) {
            this->min_m = min;
        }


    };

    /**
     * @defgroup Weibull Weibull Distribution
     * 
     * @ingroup Distributions
     * 
     */

    /**
     * @ingroup Weibull 
     * 
     * @param x
     * @param shape
     * @param scale
     * @return 
     */
    template<class T>
    const T pweibull(const T &x, const T &shape, const T &scale) {


        if (x >= T(0)) {
            return T(1.0) - std::exp(T(-1) * std::pow((x / scale), shape));
        }

        return T(0);
    }

    /**
     * 
     * @ingroup Weibull 
     * 
     * @param x
     * @param shape
     * @param scale
     * @return 
     */
    template<class T>
    const T dweibull(const T &x, const T &shape, const T &scale) {

        if (x >= T(0)) {
            return (shape / scale)*std::pow(x / scale, shape - T(1)) * std::exp(T(-1) * std::pow((x / scale), shape));
        }

        return T(0);
    }

    /**
     * @ingroup Weibull 
     * 
     * @param probability
     * @param shape
     * @param scale
     * @return 
     */
    template<class T>
    const T qweibull(const T &probability, const T &shape, const T &scale) {
        T result = T(0);
        T guess = 4;
        T high = T(100);
        T low = T(0);
        T x = guess;
        T xNew = guess;
        T error, pdf, dx = T(1.0);
        int i = 0;


        while (std::fabs(dx) > T(1e-10) && i++ < MAX_ITERATIONS) {//need to work this out beter.

            error = pweibull<T > (x, shape, scale) - probability;

            if (error < T(0.0)) {
                low = x;
            } else {
                high = x;
            }

            pdf = dweibull<T > (x, shape, scale);

            if (pdf != T(0.0)) {
                dx = error / pdf;
                xNew = x - dx;
            }

            if (xNew < low || xNew > high || pdf == T(0.0)) {
                xNew = (low + high) / T(2.0);
                dx = xNew - x;
            }
            x = xNew;

        }

        result = x;

        return result;
    }

    /**
     * @ingroup Weibull 
     * 
     * @param shape
     * @param scale
     * @return 
     */
    template<class T>
    const T rweibull(const T &shape, const T &scale) {
        T x = ((T) rand() / (T) RAND_MAX);
        T ret = qweibull<T> (x, shape, scale);
        return ret;
    }

    template<class T>
    class WiebullDistribution : public DistributionBase<T> {
    };

    //Wilcoxon Rank Sum Statistic

    //    template<class T>
    //    T pwilcox(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T qwilcox(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T dwilcox(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T rwilcox(const T &x, const T &y, const T &z) {
    //    }
    //
    //
    //    //Wilcoxon Signed Rank Statistic
    //
    //    template<class T>
    //    T psignrank(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T qsignrank(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T dsignrank(const T &x, const T &y, const T &z) {
    //    }
    //
    //    template<class T>
    //    T rsignrank(const T &x, const T &y, const T &z) {
    //    }




}//atl


#undef  DEFAULT_EPSILON



#endif	/* Distribution_HPP */

