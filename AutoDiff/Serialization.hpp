/* 
 * File:   Serialization.hpp
 * Author: matthewsupernaw
 *
 * Created on June 23, 2014, 3:33 PM
 */

/**
 *
 * @author  Matthew R. Supernaw
 *
 * Public Domain Notice
 * National Oceanic And Atmospheric Administration
 *
 * This software is a "United States Government Work" under the terms of the
 * United States Copyright Act.  It was written as part of the author's official
 * duties as a United States Government employee and thus cannot be copyrighted.
 * This software is freely available to the public for use. The National Oceanic
 * And Atmospheric Administration and the U.S. Government have not placed any
 * restriction on its use or reproduction.  Although all reasonable efforts have
 * been taken to ensure the accuracy and reliability of the software and data,
 * the National Oceanic And Atmospheric Administration and the U.S. Government
 * do not and cannot warrant the performance warrant the performance or results
 * that may be obtained by using this  software or data. The National Oceanic
 * And Atmospheric Administration and the U.S. Government disclaim all
 * warranties, express or implied, including warranties of performance,
 * merchantability or fitness for any particular purpose.
 *
 * Please cite the author(s) in any work or product based on this material.
 *
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

