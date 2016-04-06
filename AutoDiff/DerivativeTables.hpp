/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   DerivativeTables.hpp
 * Author: matthewsupernaw
 *
 * Created on February 18, 2016, 10:50 AM
 */

#ifndef DERIVATIVETABLES_HPP
#define DERIVATIVETABLES_HPP

#include "third_party/sparsehash/src/google/dense_hash_map"
#include "third_party/sparsehash/src/google/sparse_hash_map"
namespace atl {

    template<class T>
    class dense_map_wrapper : public google::dense_hash_map<uint32_t, T> {
    public:

        dense_map_wrapper() {
            this->set_empty_key(0);
        }
    };
    
    template<class T>
    class dense_map_wrapper_long_key : public google::dense_hash_map<size_t, T> {
    public:

        dense_map_wrapper_long_key() {
            this->set_empty_key(0);
        }
    };
    

    class index_exception : public std::exception {

        virtual const char* what() const throw () {
            return "index out of bounds.";
        }
    };

    template<typename T>
    class FirstOrderTable {
        uint32_t min_key;
        uint32_t max_key;
        uint32_t size;

        std::vector<T> data;
    public:

        FirstOrderTable() {

        }

        FirstOrderTable(uint32_t min_index, uint32_t max_index) :
        min_key(min_index), max_key(max_index) {
            size = max_index - min_index;
            data.resize(size);
        }

        FirstOrderTable(const FirstOrderTable<T>& other) :
        min_key(other.min_key), max_key(other.max_key), size(other.size), data(other.data) {
        }

        T& operator[](uint32_t key) {
            return data[key - min_key];
        }

        T& at(uint32_t key) {
            if (key < min_key || key > max_key) {
                throw index_exception();
            }
            return data[key - min_key];
        }

        T get(uint32_t key) {

        }

        void clear() {
            data.clear();
        }

        inline void Reset(uint32_t min_index, uint32_t max_index) {
            this->max_key = max_index;
            this->min_key = min_index;
            size = max_index - min_index;
            data.resize(size + 1);
            for (int i = 0; i < data.size(); i++) {
                data[i] = 0.0;
            }
        }

    };

    template<typename T>
    class SecondOrderTable {
        uint32_t min_key;
        uint32_t max_key;
        uint32_t size;

        std::vector<dense_map_wrapper<T> > data;

    public:

        SecondOrderTable() {
        }

        SecondOrderTable(uint32_t min_index, uint32_t max_index) :
        min_key(min_index), max_key(max_index) {
            size = max_index - min_index;
            data.resize(size + 1);

            //            for (int i = 0; i < data.size(); i++) {
            //                data[i].set_empty_key(0);
            //                set_empty[i] = true;
            //            }
            std::cout << data.size() << std::endl;
        }

        SecondOrderTable(const SecondOrderTable<T>& other) :
        min_key(other.min_key), max_key(other.max_key), size(other.size), data(other.data) {
        }

        T& operator()(uint32_t key1, uint32_t key2) {
            return data[key1 - min_key][key2];
        }

        inline T& at(uint32_t key1, uint32_t key2) {
            if (key1 < min_key || key1 > max_key || key2 < min_key || key2 > max_key) {
                throw index_exception();
            }
            return data[key1 - min_key][key2];
        }
        //

        inline T get(uint32_t key1, uint32_t key2) {
            if (key1 < min_key || key1 > max_key || key2 < min_key || key2 > max_key) {
                throw index_exception();
            }
            typename dense_map_wrapper<T>::iterator it;
            it = data[key1 - min_key].find(key2);
            if (it != data[key1 - min_key].end()) {
                return (*it).second;
            } else {
                return 0;
            }
        }

        void clear() {
            data.clear();
        }

        inline void Reset(uint32_t min_index, uint32_t max_index) {
            this->max_key = max_index;
            this->min_key = min_index;
            size_t nsize = max_index - min_index;


            //            this->clear();
            data.resize(nsize + 1);
            for (int i = 0; i < data.size(); i++) {

                data[i].clear_no_resize();
            }
            size = nsize;
        }
    };

    template<typename T>
    class ThirdOrderTable {
    public:
        uint32_t min_key;
        uint32_t max_key;
        uint32_t size;
        std::vector<dense_map_wrapper<dense_map_wrapper<T> > > data;
        typedef typename std::vector<dense_map_wrapper<dense_map_wrapper<T> > >::iterator I_iterator;
        typedef typename dense_map_wrapper<dense_map_wrapper<T> >::iterator J_iterator;
        typedef typename dense_map_wrapper<T>::iterator K_iterator;
        I_iterator i_iter;
        J_iterator j_iter;

        ThirdOrderTable() {
        }

        ThirdOrderTable(uint32_t min_index, uint32_t max_index) :
        min_key(min_index), max_key(max_index) {
            size = max_index - min_index;
            data.resize(size + 1);
        }

        ThirdOrderTable(const ThirdOrderTable<T>& other) :
        min_key(other.min_key), max_key(other.max_key), size(other.size), data(other.data) {
        }

        inline T& operator()(uint32_t key1, uint32_t key2, uint32_t key3) {
            return data[key1 - min_key][key2][key3];
        }

        inline T& at(uint32_t key1, uint32_t key2, uint32_t key3) {
            if (key1 < min_key || key1 > max_key ||
                    key2 < min_key || key2 > max_key ||
                    key3 < min_key || key3 > max_key) {
                throw index_exception();
            }
            return data[key1 - min_key][key2][key3];
        }
        //

        inline T get(uint32_t key1, uint32_t key2, uint32_t key3) {

            typename dense_map_wrapper<dense_map_wrapper<T> >::iterator it;
            typename dense_map_wrapper< T>::iterator jt;
            it = data[key1 - min_key].find(key2);
            if (it != data[key1 - min_key].end()) {
                jt = (*it).second.find(key3);
                if (jt != (*it).second.end()) {
                    return (*jt).second;
                } else {
                    return 0;
                }
            } else {
                return 0;
            }
        }

        I_iterator getI(uint32_t i) {
            //            return &data[i - min_key];
            I_iterator iter = data.begin();
            std::advance(iter, (i - min_key));
            return iter;
        }

        J_iterator getJ(uint32_t i, uint32_t j) {
            typename dense_map_wrapper<dense_map_wrapper<T> >::iterator it;
            return data[i - min_key].find(j);
        }

        I_iterator endi() {
            return data.end();
        }

        J_iterator endj(uint32_t i) {
            return data[i - min_key].end();
        }

        inline void clear() {
            data.clear();
        }

        inline void Reset(uint32_t min_index, uint32_t max_index) {
            this->max_key = max_index;
            this->min_key = min_index;
            size_t nsize = max_index - min_index;


            //            this->clear();
            data.resize(nsize + 1);
            for (int i = 0; i < data.size(); i++) {
                typename dense_map_wrapper<dense_map_wrapper<T> >::iterator it;
                //                it = data[i].begin();
                for (it = data[i].begin(); it != data[i].end(); ++it) {
                    (*it).second.clear_no_resize();
                }
                data[i].clear_no_resize();
            }
            size = nsize;
        }
    };

}


#endif /* DERIVATIVETABLES_HPP */

