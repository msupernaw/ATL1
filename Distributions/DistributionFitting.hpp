/* 
 * File:   DistributionFitting.hpp
 * Author: matthewsupernaw
 *
 * Created on December 9, 2014, 1:20 PM
 */

#ifndef DISTRIBUTIONFITTING_HPP
#define	DISTRIBUTIONFITTING_HPP

#include "Distributions.hpp"


namespace atl {
    
    
    template<class T>
    struct FitResults{
        atl::DistributionBase<T>* distribution;
        
    };

}


#endif	/* DISTRIBUTIONFITTING_HPP */

