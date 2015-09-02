/* 
 * File:   ATL.hpp
 * Author: matthewsupernaw
 *
 * Created on November 14, 2014, 12:09 PM
 */

#ifndef ATL_HPP
#define	ATL_HPP

#include "Containers/Containers.hpp"
#include "Containers/Array.hpp"
#include "Containers/Vector.hpp"
#include "Containers/VectorAdd.hpp"
#include "Containers/VectorSubtract.hpp"
#include "Containers/VectorMultiply.hpp"
#include "Containers/VectorDivide.hpp"
#include "Containers/VectorDotProduct.hpp"
#include "Containers/Matrix.hpp"
#include "Containers/MatrixAdd.hpp"
#include "Containers/MatrixSubtract.hpp"
#include "Containers/MatrixDivide.hpp"
#include "Containers/MatrixMultiply.hpp"
#include "Containers/Math/MatrixMath.hpp"
#include "Containers/ConcurrentOperators.hpp"


#include "Distributions/Distributions.hpp"
#include "SpecialFunctions/Functions.hpp"
#include "Statistics/Descriptive.hpp"
#include "Utilities/IO/Console.hpp"
#include "Utilities/IO/IOStream.hpp"


#include "AutoDiff/AutoDiff.hpp"

//#include "Optimization/GradientBased/FunctionMinimizer.hpp"



//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
//! Namespace of the \b Analytics Template Library.
namespace atl {
}
//*************************************************************************************************


//**Mainpage***************************************************************************************
/*!\mainpage Main Page
 * 
 *This is the API documentation for the <b> Analytics Template Library (ATL)</b>. It gives an overview of the  
 *individual modules included in \b ATL. \b ATL is a high performance scientific computing and parameter 
 *estimation library written in C++.  
 * 
 * \section table_of_contents Table of Contents
 * <ul>
 *      <li> \ref installation </li>  
 *      <li> \ref tutorials </li> 
 *      <ul>
 *          <li>\ref containers</li>
 *          <ul>
 *               <li> \ref arrays </li>
 *               <li> \ref vectors </li>
 *               <li> \ref matrices </li>
 *               <li> \ref mixed_ops </li>
 *         </ul>
 *          <li> \ref auto_diff</li>
 *          <ul>
 *               <li> \ref variable</li>
 *               <li> \ref forward_mode</li>
 *               <li> \ref reverse_mode </li>
 *         </ul>
 *          <li> \ref optimization</li>
 *          <ul>
 *               <li> \ref gradient_based</li>
 *               <li> \ref stochastic </li>
 *               <li> \ref evolutionary </li>
 *         </ul>
 *          <li> \ref distributions_special_functions</li>
 *          <ul>
 *               <li> \ref distributions_objects</li>
 *               <li> \ref distributions_functions </li>
 *         </ul>
 *          <li> \ref statistics</li>
 *          <ul>
 *               <li> \ref descriptive</li>
 *         </ul>
 *          <li> \ref concurrency</li>
 *          <ul>
 *               <li> \ref threads</li>
 *               <li> \ref mpi</li>
 *         </ul>
 *         <li> \ref utilities</li>
 *          <ul>
 *               <li> \ref io</li>
 *               <li> \ref serialization</li>
 *         </ul>
 *         <li> \ref benchmarks</li>
 *         <li> \ref futurework</li>
 *      </ul>
 * </ul>
 */

//**Installation*****************************************************************
/*!\page installation Installation
 * Previous: \ref index <br><br>
 * \b ATL is designed to work in an multi-threaded environment and requires that the
 * <a href="https://www.threadingbuildingblocks.org">Intel Thread Building Blocks</a> 
 * library to be installed on your computer. Please note,  \b ATL will still compile and work without 
 * the <a href="https://www.threadingbuildingblocks.org">Intel Thread Building Blocks</a> , 
 * however multi-threading with \ref reverse_mode Automatic Differentiation will not be 
 * supported.   
 * 
 * <br>In addition, to make  and install \b ATL, you must also have 
 * <a href="http://www.cmake.org/">CMake</a> 
 * installed. <br>
 * 
 * Now that you have met the requirements for building and installing \b ATL, 
 * the rest is simple. From the command line:<br>
 * 
 * \code
 *  cd ATL
 *  mkdir build
 *  cd build
 *  cmake ..
 *  make
 *  make install
 * \endcode
 * 
 * <br>That's it! ATL should now be successfully installed on your computer.
 * <br><br>Next: \ref tutorials
 */

//**Tutorials*****************************************************************
/*!\page tutorials Tutorials
 * Previous: \ref installation <br><br>
 * Welcome to the tutorials section! These tutorials are designed to give you a 
 * basic introduction to \b ATL and it's modules. Each section will describe the 
 * module being introduced and show specific examples of its use.
 * 
 *  <br><br>Next: \ref containers
 * 
 */

//**Containers*****************************************************************
/*!\page containers Containers
 * Previous: \ref tutorials <br><br>
 * The \b ATL containers library is an expression template based library.
 * There are three container types in \b ATL, Array, Vector, and Matrix. All
 * are dynamic containers, and all support basic arithmetic operations. Also,
 * mixed precision types are also supported via type promotion. Furthermore,
 * all containers provide overloaded functions from the cmath library.
 * <br><br>
 * As mentioned above, the \b ATL Containers package relies on expression templates for 
 * efficiency. Expression templates is a template meta-programming technique in
 * which templates are used to represent portions of a mathematical expression. 
 * This becomes easy to implement using the operator overloading features of the 
 * C++ programming language. Expression templates are beneficial because they 
 * eliminate the need for temporary variables during the expression evaluation.
 * <br><br>A more detailed description of expression templates can be found  
 * <a href="http://en.wikipedia.org/wiki/Expression_templates">here.</a>
 * 
 * 
 *   <br><br>Next: \ref arrays
 * 
 */

//**Arrays*****************************************************************
/*!\page arrays Arrays
 * Previous: \ref containers <br><br>
 *  Arrays are simple, dynamic, multi-dimensional structures used to hold primitives and \b ATL 
 * \ref variable types. Arrays can have up to seven dimensions. 
 * 
 * <br>Instantiation of an Array structure:
 * \code
    atl::Array<double> a(10);                     //a 1d array of length 10 of type double.
    a = 3.1459;                                   //set all elements of a to 3.1459

    atl::Array<atl::Variable<double> > b(10, 10); //a 2d array of type atl::Variable<double>
    b = atl::Variable<double>(10.0);              //set all elements of b to 10.0
 * \endcode
 * <b>*</b>atl::Variable are part of the \ref auto_diff package.
 *
 * <br>Some Array operations:
 * \code
    atl::Array<atl::Variable<double> > c = a + b; // returns a 1d array of type atl::Variable<double>
    std::cout << "c = " << c << std::endl;

    atl::Array<atl::Variable<double> > d = a - b; // returns a 1d array of type atl::Variable<double>
    std::cout << "d = " << d << std::endl;

    atl::Array<atl::Variable<double> > e = a + b; // returns a 1d array of type atl::Variable<double>
    std::cout << "e = " << e << std::endl;

    atl::Array<atl::Variable<double> > f = a*b;   // returns a 1d array of type atl::Variable<double>
    std::cout << "f = " << f << std::endl;

    atl::Array<atl::Variable<double> > g = a / b; // returns a 1d array of type atl::Variable<double>
    std::cout << "g = " << g << std::endl;
 * \endcode
 * <br>\b Output:
 * \code 
c = [0]:13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 

d = [0]:-6.8541 -6.8541 -6.8541 -6.8541 -6.8541 -6.8541 -6.8541 -6.8541 -6.8541 -6.8541 

e = [0]:13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 13.1459 

f = [0]:31.459 31.459 31.459 31.459 31.459 31.459 31.459 31.459 31.459 31.459 

g = [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 

 *  \endcode
 *<br>Arrays can also be nested via template parameters:
 *\code
    atl::Array< atl::Array<atl::Variable<double> > > aa(10);
    aa = g;                                                  //aa is now a 10 x 10 array with each row equal to g.
    std::cout<<"aa = "<<aa;
 * \endcode
 * <br>\b Output:
 * \code
 aa = [0]:[0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 
 [0]:0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 0.31459 

 *\endcode
 * 
 * For more, see  atl::Array
 * <br><br>
 * <br><br>Next: \ref vectors
 */

//**Vectors*****************************************************************
/*!\page vectors Vectors
 * Previous: \ref arrays <br><br>
 * Vectors are similar to Arrays, but have only one dimension. They are dynamic 
 * structures and are used to hold primitive or \b ATL \ref variable types.
 * <br><br>Instantiation of a Vector:
 * \code
    atl::Vector<long double> a(10);                 //create a vector of type long double of length 10
    for (int i = 0; i < a.Size(0); i++) {
        a(i) = atl::runif<long double>(0.0, 1.0);   //set elements of a to a random value between 1-0
    }

    atl::Vector<atl::Variable<long double> > b(10); //create a vector of type atl::Variable<long double> of length 10
    b = atl::Variable<long double>(1.0);            //set all elements of b to 1.0
 * \endcode
 *  <b>*</b>atl::Variable are part of the \ref auto_diff package.
 *  <br><b>**</b>atl::runif is part of the \ref distributions_functions package.
 * <br><br>Vectors support the same operations as Arrays:
 * \code
    atl::Vector<atl::Variable<long double> > c = a + b;
    std::cout<<c<<"\n";
 * \endcode
 * <br>\b Output:
 * \code
 [0]:1.00001 1.13154 1.75561 1.45865 1.53277 1.21896 1.04704 1.67886 1.6793 1.93469 
 * \endcode
 * <br><br>
 * For more, see  atl::Vector
 * <br><br>
 * <br><br>Next: \ref matrices <br><br>
 */

//**Matrices*****************************************************************
/*!\page matrices Matrices
 */

//**MixedOps*****************************************************************
/*!\page mixed_ops Mixed Operations
 */


//**AutomaticDifferentiation****************************************************
/*!\page auto_diff Automatic Differentiation
 */

//**AutomaticDifferentiation****************************************************
/*!\page variable Variable
 */
//**AutomaticDifferentiation****************************************************
/*!\page forward_mode Forward Mode
 */

//**AutomaticDifferentiation****************************************************
/*!\page reverse_mode Reverse Mode
 */

//**Optimization****************************************************************
/*!\page optimization Optimization
 */

//**Optimization****************************************************************
/*!\page gradient_based Gradient Based 
 */

//**Optimization****************************************************************
/*!\page stochastic Stochastic
 */

//**Optimization****************************************************************
/*!\page evolutionary Evolutionary
 */

//**DistributionsAndSpecialFunctions********************************************
/*!\page distributions_special_functions Distributions and Special Functions
 *Every distribution in \b ATL has four functions and a class derived from
 * class atl::DistributionBase.<br> Distribution function root names are preceeded 
 * by a letter indicating the function type:<br><br>
 * <ul>
 * <li><b>p</b> for "probability", the cumulative distribution function(cdf)</li>
 * <li><b>q</b> for "quantile", the inverse of the the cumulative distribution function(cdf)</li>
 * <li><b>d</b> for "density", the probability density function</li>
 * <li><b>r</b> for "random", a random variable from the specified distribution</li>
 * </ul>
 * 
 *<br>Available Functions And Classes:<br><br>
 * <table style="width:75%">
 *<tr><td><b>Distribution Name</b></td><td><b>Cumulative Distribution Functions </b></td><td><b>Inverse Cumulative Distribution Functions </b></td><td><b>Probability Density Functions </b></td><td><b>Random Functions </b></td><td><b>Class Name</b></td>
 *<tr><td>Beta</td><td>pbeta</td><td>qbeta</td><td>dbeta</td><td>rbeta</td><td>BetaDistribution</td>
 *<tr><td>Binomial</td><td>pbinom</td><td>qbinom</td><td>dbinom</td><td>rbinom</td><td>BinomialDistribution</td>
 *<tr><td>Cauchy</td><td>pcauchy</td><td>qcauchy</td><td>dcauchy</td><td>rcauchy</td><td>CauchyDistribution</td>
 *<tr><td>Chi-Square</td><td>pchisq</td><td>qchisq</td><td>dchisq</td><td>rchisq</td><td>ChiSquaredDistribution</td>
 *<tr><td>Exponential</td><td>pexp</td><td>qexp</td><td>dexp</td><td>rexp</td><td>ExponentialDistribution</td>
 *<tr><td>F</td><td>pf</td><td>qf</td><td>df</td><td>rf</td><td>FDistribution</td>
 *<tr><td>Gamma</td><td>pgamma</td><td>qgamma</td><td>dgamma</td><td>rgamma</td><td>GammaDistribution</td>
 *<tr><td>Geometric</td><td>pgeom</td><td>qgeom</td><td>dgeom</td><td>rgeom</td><td>GeometricDistribution</td>
 *<tr><td>Hypergeometric</td><td>phyper</td><td>qhyper</td><td>dhyper</td><td>rhyper</td><td>HypergeometricDistribution</td>
 *<tr><td>Logistic</td><td>plogis</td><td>qlogis</td><td>dlogis</td><td>rlogis</td><td>LogisticDistribution</td>        
 *<tr><td>Log Normal</td><td>plnorm</td><td>qlnorm</td><td>dlnorm</td><td>rlnorm</td><td>LogNormalDistribution</td>
 *<tr><td>Negative Binomial</td><td>pnbinom</td><td>qnbinom</td><td>dnbinom</td><td>rnbinom</td><td>NegativeBinomialDistribution</td>
 *<tr><td>Normal</td><td>pnorm</td><td>qnorm</td><td>dnorm</td><td>rnorm</td><td>NormalDistribution</td>
 *<tr><td>Poisson</td><td>ppois</td><td>qpois</td><td>dpois</td><td>rpois</td><td>PoissonDistribution</td>     
 *<tr><td>Student's t</td><td>pt</td><td>qt</td><td>dt</td><td>rt</td><td>StudentsTDistribution</td>
 *<tr><td>Uniform</td><td>punif</td><td>qunif</td><td>dunif</td><td>runif</td><td>UniformDistribution</td>
 *<tr><td>Weibull</td><td>pweibull</td><td>qweibull</td><td>dweibull</td><td>rweibull</td><td>WeibullDistribution</td>
 * </table>
 */

//**DistributionsAndSpecialFunctions********************************************
/*!\page distributions_objects Distributions Objects
 */

//**DistributionsAndSpecialFunctions********************************************
/*!\page distributions_functions Distributions Functions
 */

//**Statistics********************************************
/*!\page statistics Statistics
 */

//**Statistics********************************************
/*!\page descriptive Descriptive
 */

//**Concurrency********************************************
/*!\page concurrency Concurrency
 */

//**Concurrency********************************************
/*!\page threads Threads
 */

//**Concurrency********************************************
/*!\page mpi MPI
 */

//**Utilities********************************************
/*!\page utilities Utilities
 */


//**Utilities********************************************
/*!\page io IO
 */

//**Utilities********************************************
/*!\page serialization Serialization
 */



//**FutureWork********************************************
/*!\page futurework Future Work
 */


#endif	/* ATL_HPP */

