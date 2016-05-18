

/* 
 * File:   Config.hpp
 * Author: matthewsupernaw
 *
 * Created on February 10, 2016, 9:16 AM
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


#ifndef CONFIG_HPP
#define CONFIG_HPP

//#define USE_BOOST
//#define USE_GOOGLE_SET
#ifdef USE_BOOST
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#elif defined(USE_GOOGLE_SET)
#include "third_party/sparsehash/src/google/dense_hash_set"

template<class T>
class dense_set_wrapper : public google::dense_hash_set<T> {
public:

    dense_set_wrapper() {
        this->set_empty_key(NULL);
    }
    void reserve(size_t v){
        
    }
};

#else
#include <set>
#include "../Utilities/flat_map.hpp"
#include "../Utilities/flat_set.hpp"
#endif


#ifdef USE_BOOST
#define IDSet boost::container::flat_set
#elif defined(USE_GOOGLE_SET)
#define IDSet dense_set_wrapper
#else
#define IDSet flat_set
#endif

#endif /* CONFIG_HPP */

