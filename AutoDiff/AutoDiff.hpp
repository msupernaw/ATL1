/* 
 * File:   ET4AD.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 11:57 AM
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


#ifndef ET4AD_HPP
#define	ET4AD_HPP
#include "Variable.hpp"
#include "Scalar.hpp"
#include "Add.hpp"
#include "Subtract.hpp"
#include "Multiply.hpp"
#include "Divide.hpp"
#include "Sin.hpp"
#include "Cos.hpp"
#include "Tan.hpp"
#include "ASin.hpp"
#include "ACos.hpp"
#include "ATan.hpp"
#include "Sqrt.hpp"
#include "Pow.hpp"
#include "Log.hpp"
#include "Log10.hpp"
#include "Exp.hpp"
#include "Sinh.hpp"
#include "Cosh.hpp"
#include "Tanh.hpp"
#include "Fabs.hpp"
#include "Floor.hpp"
#include "Ceil.hpp"


namespace atl{
  
 /**
 * Returns the maximum between a and b in a continuous manner using:
 * 
 * (a + b + |a - b|) / 2.0;
 * 
 * @param a
 * @param b
 * @return 
 */
template <typename T>
inline const atl::Variable<T> max(const atl::Variable<T>& a, const atl::Variable<T>& b) {
    return (a + b +  atl::fabs(a - b)) / 2.0;
}

/**
 * Returns the minimum between a and b in a continuous manner using:
 * 
 * (a + b - |a - b|) / 2.0;
 * 
 * @param a
 * @param b
 * @return 
 */
template <typename T>
inline const atl::Variable<T> min(const atl::Variable<T>& a, const atl::Variable<T>& b) {
    return (a + b - atl::fabs(a - b)) / 2.0;
} 
  
  
  /**
 * Returns the maximum between a and b in a continuous manner using:
 * 
 * (a + b + |a.GetValue() - b.GetValue()|) / 2.0;
 * 
 * @param a
 * @param b
 * @return 
 */
template <typename T>
inline const atl::Variable<T> max2(const atl::Variable<T>& a, const atl::Variable<T>& b) {
    return (a + b +  std::fabs(a.GetValue() - b.GetValue())) / 2.0;
}

/**
 * Returns the minimum between a and b in a continuous manner using:
 * 
 * (a + b - |a.GetValue() - b.GetValue()|) / 2.0;
 * 
 * @param a
 * @param b
 * @return 
 */
template <typename T>
inline const atl::Variable<T> min2(const atl::Variable<T>& a, const atl::Variable<T>& b) {
    return (a + b - std::fabs(a.GetValue() - b.GetValue())) / 2.0;
}
  
}

//
//typedef atl::Variable<double> variable;
//typedef atl::Variable<long double> variablel;
//typedef atl::Variable<float> variablef;








#endif	/* ET4AD_HPP */

