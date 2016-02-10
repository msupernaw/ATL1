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

#ifdef USE_BOOST
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#else
#include <set>
#include "../Utilities/flat_map.hpp"
#include "../Utilities/flat_set.hpp"
#endif

#ifdef USE_BOOST
#define IDSet boost::container::flat_set
#else
#define IDSet flat_set
#endif

#endif /* CONFIG_HPP */

