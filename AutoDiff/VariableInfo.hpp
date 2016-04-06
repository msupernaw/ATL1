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

#ifdef ATL_VARIABLE_INFO_USE_MEMORY_POOL
#include "PoolAllocator.hpp"
#endif



#include <cassert>
#include <unordered_map>

namespace atl {

    
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

        VariableIdGenerator() : _id(1), available_size(0) {
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

    template<typename REAL_T>
    class VariableInfo
#ifdef ATL_VARIABLE_INFO_USE_MEMORY_POOL
    : public atl::PoolAllocator<VariableInfo<REAL_T> > {
#else 
    {
#endif
        public:
     
        static std::mutex vinfo_mutex_g;
        static std::vector<VariableInfo<REAL_T>* > freed;
        REAL_T dvalue;
        REAL_T vvalue;
        std::atomic<int> count;
        std::atomic<int> dependence_level;
        std::atomic<int> is_dependent;
//        IDSet<atl::VariableInfo<REAL_T>* > dependencies;
        uint32_t id;
        uint32_t push_start = 0;//the beginning of nonlinear interaction
        std::string name;
        int push_count =0;
        bool has_nl_interaction =false;
        bool is_nl = false;
        
        VariableInfo() : dvalue(0.0), vvalue(0.0), count(1), dependence_level(1), is_dependent(0), id(VariableIdGenerator::instance()->next()) {


        }

        inline void Aquire() {
            count++;
        }

        virtual ~VariableInfo() {
        }


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
            this->dvalue = 0;
            this->dependence_level = 1;
            this->is_dependent = 0;
            this->push_count = 0;
            push_start = 0;
            this->has_nl_interaction = false;
            this->is_nl = false;
        }

        static void FreeAll() {
#pragma unroll
            for (int i = 0; i < freed.size(); i++) {
                VariableIdGenerator::instance()->release(freed[i]->id);
                freed[i]->vvalue = 0;
                delete freed[i];
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

