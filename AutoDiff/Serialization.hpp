/* 
 * File:   Serialization.hpp
 * Author: matthewsupernaw
 *
 * Created on June 23, 2014, 3:33 PM
 */

#ifndef SERIALIZATION_HPP
#define	SERIALIZATION_HPP

#include "GradientStructure.hpp"


namespace atl{
  
    
    
    /**
     * 
     * Serialize a GradientStack to a machine independent string representation.
     * Base types must be primitives.
     * 
     * @param gradient_stack
     * @return 
     */
    template<class REAL_T>
    const std::string Serialize(const atl::GradientStack<REAL_T>& gradient_stack){
        
    }
    
    /**
     * Deserialize a machine independent string representation of a 
     * GradientStack. 
     * 
     * @param gradient_stack
     * @return 
     */
    template<class REAL_T>
    const atl::GradientStack<REAL_T> Deserialize(const std::string& gradient_stack){
        
    }
    
    /**
     * Serialize a Variable to a machine independent string representation.
     * Base types must be primitives.
     * @param variable
     * @return 
     */
    template<class REAL_T>
    const std::string Serialize(const atl::Variable<REAL_T>& variable){
        
    } 
    
    /**
     * 
     * Deserialize a machine independent string representation of a 
     * Variable. 
     *  
     * @param variable
     * @return 
     */
    template<class REAL_T>
    const atl::Variable<REAL_T> Deserialize(const std::string& variable){
        
    } 
}


#endif	/* SERIALIZATION_HPP */

