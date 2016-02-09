/* 
 * File:   DerivativePointer.hpp
 * Author: matthewsupernaw
 *
 * Created on April 8, 2015, 12:12 PM
 */

#ifndef DERIVATIVEPOINTER_HPP
#define DERIVATIVEPOINTER_HPP
#include <atomic>
#include <mutex>
#include <stack>
#include <memory>
#include "../Utilities/Platform.hpp"
#include "../Utilities/flat_set.hpp"
#include "../Utilities/MemoryPool.hpp"
#include "../Utilities/compressed_vector.hpp"
#include "PoolAllocator.hpp"
#include "third_party/dlmalloc/malloc.h"
#include <google/sparse_hash_map>
#include <cassert>
#include <unordered_map>
//#include "third_party/cache-table-0.2/mm/cache_map.hpp"
namespace atl {

    //    template <class T>
    //    class PoolAllocator {
    //    public:
    //
    //        static void* operator new(size_t size) {
    //            return s_memPool.malloc();
    ////                        return malloc(size);
    //        }
    //
    //        static void operator delete(void* deletable, size_t size) {
    //            //don't delete null pointers
    //            if (deletable)
    //                s_memPool.free(deletable);
    ////                        free(deletable);
    //        }
    //
    //    protected:
    //
    //        ~PoolAllocator() {
    //        }
    //
    //    private:
    //        //each FastAllocator specialization has it's own memory pool
    //        static util::MemoryPool<T> s_memPool;
    //    };
    //
    //    //the static variable s_memPool is defined here. It's constructor is passed the object size.
    //    template <class T>
    //    util::MemoryPool<T> PoolAllocator<T>::s_memPool(500000);

    /*!
     * Creates a unique identifier for variables. Identifiers are recyclable.
     * @return 
     */
    class VariableIdGenerator {
        std::stack<uint32_t> available;
        std::atomic<uint32_t> available_size;
        static std::mutex mutex_g;

    public:
        static std::shared_ptr<VariableIdGenerator> instance();

        const uint32_t next() {
            mutex_g.lock();
            uint32_t ret;
            if (!available.empty() > 0) {
                ret = available.top();
                available.pop();

            } else {
                ret = ++_id;
            }



            mutex_g.unlock();
            return ret; //(++_id);
        }

        void release(const uint32_t& id) {
            mutex_g.lock();
            available.push(id);
            available_size++;
            mutex_g.unlock();
        }

        const uint32_t current() {
            return _id;
        }

        //    private:

        VariableIdGenerator() : _id(0), available_size(0) {
        }

        std::atomic<uint32_t> _id;
    };




    std::mutex VariableIdGenerator::mutex_g;
    static std::shared_ptr<VariableIdGenerator> only_copy;

    inline std::shared_ptr<VariableIdGenerator>
    VariableIdGenerator::instance() {

        if (!only_copy) {
            only_copy = std::make_shared<VariableIdGenerator>();
        }

        return only_copy;
    }

    /*!
     * Creates a unique identifier for variables. Identifiers are recyclable.
     * @return 
     */
    class IndependentVariableIdGenerator {
        std::stack<uint32_t> available;
        std::atomic<uint32_t> available_size;
        static std::mutex mutex_g;

    public:
        static std::shared_ptr<IndependentVariableIdGenerator> instance();

        const uint32_t next() {
            mutex_g.lock();
            uint32_t ret;
            if (!available.empty() > 0) {
                ret = available.top();
                available.pop();

            } else {
                ret = ++_id;
            }



            mutex_g.unlock();
            return ret; //(++_id);
        }

        void release(const uint32_t& id) {
            mutex_g.lock();
            available.push(id);
            available_size++;
            mutex_g.unlock();
        }

        const uint32_t current() {
            return _id;
        }

        //    private:

        IndependentVariableIdGenerator() : _id(0), available_size(0) {
        }

        std::atomic<uint32_t> _id;
    };




    std::mutex IndependentVariableIdGenerator::mutex_g;
    static std::shared_ptr<IndependentVariableIdGenerator> only_copy2;

    inline std::shared_ptr<IndependentVariableIdGenerator>
    IndependentVariableIdGenerator::instance() {

        if (!only_copy2) {
            only_copy2 = std::make_shared<IndependentVariableIdGenerator>();
        }

        return only_copy2;
    }

    template<typename T>
    class HashMapWrapper {
        typedef google::sparse_hash_map<uint32_t, T> wrapped_map;
        wrapped_map data;

    public:
        typedef typename wrapped_map::iterator iterator;
        typedef typename wrapped_map::const_iterator const_iterator;

        HashMapWrapper() : data(1) {
            //            data.set_empty_key(0);
        }

        inline iterator find(uint32_t key) {
            return data.find(key);
        }

        inline const_iterator find(uint32_t key)const {
            return data.find(key);
        }

        inline iterator end() {
            return data.end();
        }

        inline const_iterator end() const {
            return data.end();
        }

        inline T& operator[](uint32_t key) {
            return data[key];
        }

        inline T& operator[](uint32_t key)const {
            return data[key];
        }

        void clear() {
            data.clear_no_resize();
        }


    };

    template<typename REAL_T>
    class VariableInfo : public atl::PoolAllocator<VariableInfo<REAL_T> > {
    public:
        //    static util::MemoryPool<VariableInfo<REAL_T>* > pool_g;
        //        static std::mutex mutex_g;
        static std::mutex vinfo_mutex_g;
        static std::vector<VariableInfo<REAL_T>* > freed;
        REAL_T dvalue;
        REAL_T d2value;
        REAL_T d3value;
        REAL_T vvalue;
        int count;
        std::atomic<int> occurences;

        typedef util::compressed_vector<REAL_T> HessianInfo;
        typedef util::compressed_vector<HessianInfo> ThirdOrderMixed;
        //                typedef google::sparse_hash_map<uint32_t, REAL_T> HessianInfo;
        //                typedef google::sparse_hash_map<uint32_t, HessianInfo> ThirdOrderMixed;
        //                        typedef mm::cache_map<uint32_t, REAL_T> HessianInfo;
        //                typedef mm::cache_map<uint32_t, HessianInfo> ThirdOrderMixed;
//        typedef std::unordered_map<uint32_t, REAL_T> HessianInfo;
//        typedef std::unordered_map<uint32_t, HessianInfo> ThirdOrderMixed;
        //        typedef HashMapWrapper<REAL_T> HessianInfo;
        //        typedef HashMapWrapper< HessianInfo> ThirdOrderMixed;
//        HessianInfo hessian_row;
//        ThirdOrderMixed third_order_mixed;
        typedef typename HessianInfo::iterator row_iterator;
        typedef typename ThirdOrderMixed::iterator h_iterator;
        uint32_t id;

        VariableInfo() : vvalue(0.0), dvalue(0.0), count(1), id(VariableIdGenerator::instance()->next()) {
            //            hessian_row.set_empty_key(0);
            //            third_order_mixed.set_empty_key(0);
            //            hessian_row.bucket_size(100);

        }

        inline void Aquire() {
            count++;
        }

        virtual ~VariableInfo() {
        }

//        inline REAL_T GetHessianRowValue(VariableInfo<REAL_T>* var) {
//            row_iterator it = hessian_row.find(var->id);
//            if (it != hessian_row.end()) {
//                return (*it).second;
//            } else {
//                return static_cast<REAL_T> (0.0);
//            }
//        }
//
//        inline REAL_T GetThirdOrderValue(VariableInfo<REAL_T>* a, VariableInfo<REAL_T>* b) {
//
//            h_iterator it = this->third_order_mixed.find(a->id);
//            if (it != this->third_order_mixed.end()) {
//                row_iterator jt = (*it).second.find(b->id);
//                if (jt != (*it).second.end()) {
//                    return (*jt).second;
//                } else {
//                    return static_cast<REAL_T> (0.0);
//                }
//            } else {
//                return static_cast<REAL_T> (0.0);
//            }
//        }

        inline void Release() {
            count--;

            if ((count) == 0) {
                //store this pointer in the freed list and delete when the gradient 
                //structure resets.
#ifdef ATL_THREAD_SAFE
                VariableInfo<REAL_T>::vinfo_mutex_g.lock();
                freed.push_back(this);
                VariableInfo<REAL_T>::vinfo_mutex_g.unlock();
#else
                freed.push_back(this);
#endif

            }
        }

        inline void Reset() {
//            this->hessian_row.clear();
//            this->third_order_mixed.clear();
            this->dvalue = 0;
            this->d2value = 0;
            this->d3value = 0;
            occurences = 0;
        }

        static void FreeAll() {
#pragma unroll
            for (int i = 0; i < freed.size(); i++) {
                VariableIdGenerator::instance()->release(freed[i]->id);
                VariableInfo<REAL_T>::operator delete(freed[i], sizeof (VariableInfo<REAL_T>));
                //                delete freed[i];
            }
            freed.resize(0);
        }
    };

    template<typename REAL_T>
    std::vector<VariableInfo<REAL_T>* > VariableInfo<REAL_T>::freed = [] {
        std::vector<VariableInfo<REAL_T>*> v;
        v.reserve(100000);
        return v;
    }();

    template<typename REAL_T>
    std::mutex VariableInfo<REAL_T>::vinfo_mutex_g;



}


#endif /* DERIVATIVEPOINTER_HPP */

