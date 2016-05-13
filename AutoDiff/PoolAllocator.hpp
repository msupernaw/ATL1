/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   PoolAllocator.hpp
 * Author: matthewsupernaw
 *
 * Created on January 4, 2016, 8:50 AM
 */

#ifndef POOLALLOCATOR_HPP
#define POOLALLOCATOR_HPP

#include "../Utilities/MemoryPool.hpp"

namespace atl{
    template <class T>
    class PoolAllocator {
    public:

        static void* operator new(size_t size) {
            return s_memPool.malloc();
//                        return malloc(size);
        }

        static void operator delete(void* deletable, size_t size) {
            //don't delete null pointers
            if (deletable)
                s_memPool.free(deletable);
//                        free(deletable);
        }

    protected:

        ~PoolAllocator() {
        }

    private:
        //each FastAllocator specialization has it's own memory pool
        static util::MemoryPool<T> s_memPool;
    };

    //the static variable s_memPool is defined here. It's constructor is passed the object size.
    template <class T>
    util::MemoryPool<T> PoolAllocator<T>::s_memPool(500000);



template <class T>
    class DynamicExpressionPoolAllocator {
    public:

        static void* operator new(size_t size) {
            return s_memPool.malloc();
//                        return malloc(size);
        }

        static void operator delete(void* deletable, size_t size) {
            //don't delete null pointers
            if (deletable)
                s_memPool.free(deletable);
//                        free(deletable);
        }

    protected:

        ~DynamicExpressionPoolAllocator() {
        }

    private:
        //each FastAllocator specialization has it's own memory pool
        static util::MemoryPool<T> s_memPool;
    };

    //the static variable s_memPool is defined here. It's constructor is passed the object size.
    template <class T>
    util::MemoryPool<T> DynamicExpressionPoolAllocator<T>::s_memPool(500000);

}



#endif /* POOLALLOCATOR_HPP */

