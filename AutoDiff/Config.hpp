/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Config.hpp
 * Author: matthewsupernaw
 *
 * Created on February 10, 2016, 9:16 AM
 */

#ifndef CONFIG_HPP
#define CONFIG_HPP

//#define USE_BOOST
#define USE_GOOGLE_SET
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

