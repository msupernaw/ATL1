
/* 
 * File:   PoolAllocator.hpp
 * Author: matthewsupernaw
 *
 * Created on January 4, 2016, 8:50 AM
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
    util::MemoryPool<T> PoolAllocator<T>::s_memPool(50000);



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

