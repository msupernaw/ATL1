/*
 * File:   DerivativePointer.hpp
 * Author: matthewsupernaw
 *
 * Created on April 8, 2015, 12:12 PM
 */

#ifndef DERIVATIVEPOINTER_HPP
#define	DERIVATIVEPOINTER_HPP
#include <atomic>
#include <mutex>
#include <stack>

#include "../Utilities/Platform.hpp"

#ifndef ATL_WINDOWS
#include "third_party/clfmalloc.h"
#endif

#include "../Utilities/flat_map.hpp"


namespace atl {
    
    /*!
     * Creates a unique identifier for variables. Identifiers are recyclable.
     * @return
     */
    class VariableIdGenerator {
#ifdef USE_TBB
        tbb::concurrent_vector<uint32_t> available;
#else
        std::stack<uint32_t> available;
#endif
        std::atomic<uint32_t> available_size;
        static std::mutex mutex_g;
        
    public:
        static  std::shared_ptr<VariableIdGenerator> instance();
        
        const uint32_t next() {
            mutex_g.lock();
            uint32_t ret;
            if (!available.empty()> 0) {
                ret = available.top();
                available.pop();
                
            } else {
                ret = ++_id;
            }
            
            
            
            mutex_g.unlock();
            return ret;//(++_id);
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
    
    inline  std::shared_ptr<VariableIdGenerator>
    VariableIdGenerator::instance() {
        
        if (!only_copy) {
            only_copy = std::make_shared<VariableIdGenerator>();
        }
        
        return only_copy;
    }
    
    template<typename REAL_T>
    class VariableInfo {
    public:
        static std::mutex mutex_g;
        std::atomic<uint32_t> count;
        
        static std::mutex vinfo_mutex_g;
        static std::vector<VariableInfo<REAL_T>* > freed;
        REAL_T dvalue;
        REAL_T vvalue;
        
        typedef flat_map<VariableInfo<REAL_T>*, REAL_T> HessianInfo;
        HessianInfo hessian_row;
        typedef typename HessianInfo::iterator row_iterator;
        uint32_t id;
        
        VariableInfo() : vvalue(0.0),dvalue(0.0), count(1), id(VariableIdGenerator::instance()->next()) {
            
        }
        
        inline void Aquire() {
            count++;
        }
        
        inline REAL_T GetHessianRowValue(VariableInfo<REAL_T>* var){
            row_iterator it = hessian_row.find(var);
            if(it != hessian_row.end()){
                return (*it).second;
            }else{
                return static_cast<REAL_T>(0.0);
            }
        }
        
        inline void Release() {
            count--;
            if ((count) == 0) {
                VariableInfo<REAL_T>::vinfo_mutex_g.lock();
                freed.push_back(this);
                VariableInfo<REAL_T>::vinfo_mutex_g.unlock();
            }
        }
        
        static void FreeAll() {
            for (int i = 0; i < freed.size(); i++) {
                VariableIdGenerator::instance()->release(freed[i]->id);
                delete freed[i];
            }
            freed.resize(0);
        }
        
        void* operator new(size_t size) {
            return malloc(size);
        }
        
        void operator delete(void* ptr) {
            free(ptr);
        }
    };
    
    
    template<typename REAL_T>
    std::vector<VariableInfo<REAL_T>* > VariableInfo<REAL_T>::freed;
    
    template<typename REAL_T>
    std::mutex VariableInfo<REAL_T>::vinfo_mutex_g;
    
    template<typename REAL_T>
    std::mutex VariableInfo<REAL_T>::mutex_g;
    
    
}


#endif	/* DERIVATIVEPOINTER_HPP */
