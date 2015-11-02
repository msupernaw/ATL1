/*
 * File:   GradientStructure.hpp
 * Author: matthewsupernaw
 *
 * Created on April 7, 2015, 4:05 PM
 */

#ifndef GRADIENTSTRUCTURE_HPP
#define	GRADIENTSTRUCTURE_HPP
#include <memory>
#include <vector>
#include <algorithm>
#include <functional>
#include <atomic>
#include <iostream>
#include "VariableInfo.hpp"
#include "AlignedAllocator.hpp"
#ifdef USE_TBB
#include "third_party/tbb42_20140601oss/include/tbb/concurrent_vector.h"
#endif
//#include "third_party/clfmalloc.h"
#include <fstream>
#include <cmath>
//#include <unordered_map>

//#define USE_BOOST
#ifdef USE_BOOST
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#else
#include <set>
#include "../Utilities/flat_map.hpp"
#include "../Utilities/flat_set.hpp"
#endif

#include "Variable.hpp"

//#warning add jacobian matrix calculations

template<typename T>
class Allocator {
public:
    //    typedefs
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    
public:
    //    convert an allocator<T> to allocator<U>
    
    template<typename U>
    struct rebind {
        typedef Allocator<U> other;
    };
    
public:
    
    inline explicit Allocator() {
    }
    
    inline ~Allocator() {
    }
    
    inline explicit Allocator(Allocator const&) {
    }
    
    template<typename U>
    inline explicit Allocator(Allocator<U> const&) {
    }
    
    //    address
    
    inline pointer address(reference r) {
        return &r;
    }
    
    inline const_pointer address(const_reference r) {
        return &r;
    }
    
    //    memory allocation
    
    inline pointer allocate(size_type cnt,
                            typename std::allocator<void>::const_pointer = 0) {
        return reinterpret_cast<pointer> (malloc(cnt * sizeof (T)));
    }
    
    inline void deallocate(pointer p, size_type) {
        free(p);
    }
    
    //    size
    
    inline size_type max_size() const {
        return std::numeric_limits<size_type>::max() / sizeof (T);
    }
    
    //    construction/destruction
    
    inline void construct(pointer p, const T& t) {
        new(p) T(t);
    }
    
    inline void destroy(pointer p) {
        p->~T();
    }
    
    inline bool operator==(Allocator const&) {
        return true;
    }
    
    inline bool operator!=(Allocator const& a) {
        return !operator==(a);
    }
};

template<class T>
class SetWrapper {
#ifdef USE_BOOST
    boost::container::flat_set<T, std::less<T>, Allocator<T> > s;
public:
    
    typedef typename boost::container::flat_set<T, std::less<T>, Allocator<T> >::value_type value_type;
    typedef typename boost::container::flat_set<T, std::less<T>, Allocator<T> >::iterator iterator;
    typedef typename boost::container::flat_set<T, std::less<T>, Allocator<T> >::const_iterator const_iterator;
#else
    std::set<T, std::less<T>, Allocator<T> > s;
public:
    
    typedef typename std::set<T, std::less<T>, Allocator<T> >::value_type value_type;
    typedef typename std::set<T, std::less<T>, Allocator<T> >::iterator iterator;
    typedef typename std::set<T, std::less<T>, Allocator<T> >::const_iterator const_iterator;
#endif
    
    inline std::pair<iterator, bool> insert(const value_type& value) {
        return s.insert(value);
    }
    
    inline std::pair<iterator, bool> insert(value_type&& value) {
        return s.insert(value);
    }
    
    inline iterator end() {
        return s.end();
    }
    
    inline const_iterator end() const {
        return s.end();
    }
    
    inline iterator begin() {
        return s.begin();
    }
    
    inline const_iterator begin() const {
        return s.begin();
    }
    
    inline size_t size() {
        return s.size();
    }
    
    inline void clear() {
        s.clear();
    }
    
    
};

template <class T, class Compare = std::less<T> >
struct sorted_vector {
    //    using std::vector;
    //    using std::lower_bound;
    std::vector<T> V;
    Compare cmp;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    
    iterator begin() {
        return V.begin();
    }
    
    iterator end() {
        return V.end();
    }
    
    const_iterator begin() const {
        return V.begin();
    }
    
    const_iterator end() const {
        return V.end();
    }
    
    sorted_vector(const Compare& c = Compare())
    : V(), cmp(c) {
    }
    
    template <class InputIterator>
    sorted_vector(InputIterator first, InputIterator last,
                  const Compare& c = Compare())
    : V(first, last), cmp(c) {
        std::sort(begin(), end(), cmp);
    }
    
    iterator insert(const T& t) {
        iterator i = std::lower_bound(begin(), end(), t, cmp);
        if (i == end() || cmp(t, *i))
            V.insert(i, t);
        return i;
    }
    
    size_t size() {
        return V.size();
    }
    
    void clear() {
        //        V.clear();
        V.resize(0);
    }
    
    const_iterator find(const T& t) const {
        const_iterator i = std::lower_bound(begin(), end(), t, cmp);
        return i == end() || cmp(t, *i) ? end() : i;
    }
};
#ifdef USE_BOOST
#define IDSet boost::container::flat_set
#else
#define IDSet flat_set
#endif
#define Entry Adjoint<REAL_T>



namespace atl {
    
    template<typename REAL_T>
    class AdjointDerivative {
        //        typedef std::vector<std::vector<REAL_T > > Mat;
        //        Mat hessian;
    public:
        VariableInfo<REAL_T>* dependent;
        REAL_T forward;
        REAL_T forward2;
        //        bool separable;
        
        AdjointDerivative() : dependent(NULL), forward(0.0), forward2(0.0) {
        }
        
        AdjointDerivative(const AdjointDerivative<REAL_T>& o)
        : dependent(std::move(o.dependent)),
        forward(std::move(o.forward)),
        forward2(std::move(o.forward2)) {
            
            
        }
        
        //        AdjointDerivative(AdjointDerivative<REAL_T>&& o) : dependent(std::move(o.dependent)), forward(std::move(o.forward)) {
        //            o.forward = 0.0;
        //            o.forward2 = 0.0;
        //            o.dependent = NULL;
        //        }
        
        AdjointDerivative& operator=(const AdjointDerivative<REAL_T>& other) {
            this->dependent = other.dependent;
            forward = other.forward;
            forward2 = other.forward2;
            return *this;
        }
        
        AdjointDerivative(VariableInfo<REAL_T>* d, REAL_T f) :
        dependent(d), forward(f), forward2(1) {
        }
        
        AdjointDerivative(VariableInfo<REAL_T>* d, REAL_T f, REAL_T f2) :
        dependent(d), forward(f), forward2(f2) {
            //            std::cout << "ad " << forward << " " << forward2 << "\n";
        }
        
        
    };
    
    template<typename REAL_T>
    struct Adjoint {
        VariableInfo<REAL_T>* w; //function or dependent variable.
        std::vector<AdjointDerivative<REAL_T> > entries;
        typedef typename flat_map<VariableInfo<REAL_T>*, std::vector<REAL_T> >::iterator entry2_iterator;
        //mixed second order partials used for
        //the exact hessian calculations.
        
        IDSet<VariableInfo<REAL_T>* > ids;
        std::vector<REAL_T> second_order_partials;
        
        Adjoint() : w(NULL) {
            entries.reserve(4);
        }
        
        Adjoint(const Adjoint<REAL_T>& other) : w(other.w), entries((other.entries)) {
            
        }
        
        //                Adjoint(Adjoint<REAL_T>&& other) : w(other.w), entries(std::move(other.entries)) {
        //
        //                }
        
        inline void Reset() {
            w->dvalue = 0;
            w->hessian_row.clear();
            w = NULL;
            
#pragma unroll
            for (int i = 0; i < entries.size(); i++) {
                if (NULL != entries[i].dependent) {
                    entries[i].dependent->dvalue = 0;
                    entries[i].dependent->hessian_row.reset();
                    //                    std::cout<<"count "<<entries[i].dependent->count<<"\n";
                }
            }
            entries.resize(0);
            second_order_partials.resize(0);
            ids.clear();
        }
        
    };
    
    enum DerivativeTraceLevel {
        GRADIENT = 0,
        GRADIENT_AND_HESSIAN
    };
    
    /**
     * Class to record operations. Often refered to as a "Tape". Holds a stack of
     * first and second order partial derivatives used in adjoint accumulation of
     * gradients and Hessian matrices.
     */
    template<typename REAL_T>
    class GradientStructure {
    public:
        DerivativeTraceLevel derivative_trace_level;
        std::vector<Adjoint<REAL_T> > gradient_stack;
        //        Adjoint<REAL_T>* gradient_stack;
        std::atomic<size_t> stack_current;
        bool recording;
        size_t max_stack_size;
        size_t max_initialized_size;
        
        bool gradient_computed;
        
        GradientStructure(uint32_t size = 1000000) : recording(true), stack_current(0), gradient_computed(false), derivative_trace_level(GRADIENT_AND_HESSIAN) {
            gradient_stack.resize(size);
            //            gradient_stack = (Adjoint<REAL_T>*)malloc(size);
            max_stack_size = size;
            max_initialized_size = 0;
        }
        
        /**
         * Sets the size of the stack.
         * @param size
         */
        void SetSize(size_t size) {
            //            gradient_stack.reserve(size);
        }
        
        virtual ~GradientStructure() {
            //            if (this->stack_current != 0) {
            //                this->Reset();
            //            }
            //            for (int i = 0; i < max_initialized_size; i++) {
            //                (&gradient_stack[i])->~Adjoint<REAL_T>();
            //            }
            //            //            }
            //            free(this->gradient_stack);
        }
        
        inline const uint32_t GetStartIndex(uint32_t count) {
            return stack_current.fetch_add(count, std::memory_order_relaxed);
        }
        
        inline void SetRecording(bool recording) {
            this->recording = recording;
        }
        
        /**
         * Atomic operation. Gets the next available index in the stack.
         *
         * @return
         */
        inline size_t NextIndex() {
            return stack_current.fetch_add(1, std::memory_order_relaxed);
        }
        
        inline Adjoint<REAL_T>& NextEntry() {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            if ((this->stack_current + 1) >= this->max_stack_size) {
                std::cout << "Current stack index exceeds gradient stack limits.\n" << std::flush;
                exit(0);
            }
#endif
            return this->gradient_stack[this->NextIndex()];
        }
        
        /**
         * Computes the gradient via reverse mode accumulation.
         * \image html gradient.png
         */
        inline void Accumulate() {
            gradient_computed = true;
            
            if (recording) {
                
                REAL_T w = 0.0;
                this->gradient_stack[stack_current - 1].w->dvalue = 1.0;
#pragma unroll
                for (int i = (stack_current - 1); i >= 0; i--) {
                    w = this->gradient_stack[i].w->dvalue;
                    this->gradient_stack[i].w->dvalue = 0;
                    if (w != static_cast<REAL_T> (0)) {
#pragma unroll
                        for (int j = 0; j < this->gradient_stack[i].entries.size(); j++) {
                            this->gradient_stack[i].entries[j].dependent->dvalue += w * this->gradient_stack[i].entries[j].forward;
                        }
                    }
                }
                
            }
        }
        
        /**
         * Computes the gradient and Hessian matrix via reverse mode accumulation.
         *
         * \image html hessian.png
         */
        void HessianAndGradientAccumulate() {
            typedef typename VariableInfo<REAL_T>::HessianInfo HessianInfo;
            
            if (recording) {
                
                
                if (this->derivative_trace_level == GRADIENT) {
                    this->Accumulate();
                } else {
                    REAL_T w, hw;
                    this->gradient_stack[stack_current - 1].w->dvalue = 1.0;
                    
                    for (int i = (stack_current - 1); i >= 0; i--) {
                        typename HessianInfo::iterator vit;
                        atl::VariableInfo<REAL_T>* vi = this->gradient_stack[i].w;
                        w = this->gradient_stack[i].w->dvalue;
                        this->gradient_stack[i].w->dvalue = 0;
                        if (w != static_cast<REAL_T> (0)) {
                            
                            REAL_T hii = 0.0;
                            vit = vi->hessian_row.find(vi);
                            if (vit != vi->hessian_row.end()) {
                                hii = (*vit).second;
                            }
#pragma unroll
                            for (int j = 0; j < this->gradient_stack[i].entries.size(); j++) {
                                this->gradient_stack[i].entries[j].dependent->dvalue += w * this->gradient_stack[i].entries[j].forward;
                                atl::VariableInfo<REAL_T>* vj = this->gradient_stack[i].entries[j].dependent;
                                REAL_T hij = 0.0;
                                typename HessianInfo::iterator vijt;
                                vijt = vi->hessian_row.find(vj);
                                if (vijt != vi->hessian_row.end()) {
                                    hij = (*vijt).second;
                                }
#pragma unroll
                                for (int k = 0; k < this->gradient_stack[i].entries.size(); k++) {
                                    atl::VariableInfo<REAL_T>* vk = this->gradient_stack[i].entries[k].dependent;
                                    
                                    REAL_T hjk = 0.0;
                                    typename HessianInfo::iterator vjt;
                                    vjt = vj->hessian_row.find(vk);
                                    if (vjt != vj->hessian_row.end()) {
                                        hjk = (*vjt).second;
                                    }
                                    REAL_T hik = 0.0;
                                    typename HessianInfo::iterator vkt;
                                    vkt = vi->hessian_row.find(vk);
                                    if (vkt != vi->hessian_row.end()) {
                                        hik = (*vkt).second;
                                    }
                                    REAL_T entry = 0.0;
                                    //
                                    entry += (hik * this->gradient_stack[i].entries[j].forward);
                                    entry += (hij * this->gradient_stack[i].entries[k].forward);
                                    entry += hii * this->gradient_stack[i].entries[j].forward * this->gradient_stack[i].entries[k].forward;
                                    entry += w * this->gradient_stack[i].second_order_partials[j * this->gradient_stack[i].entries.size() + k];
                                    
                                    if (entry != 0.0/* && std::fabs(entry) != std::numeric_limits<REAL_T>::infinity()*/) {
                                        REAL_T hentry = entry + hjk;
                                        vj->hessian_row[vk] = hentry;
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }
        }
        
        //        /**
        //         * Creates and adds an entry into this gradient stack for a variable(the dependent)
        //         * from a gradient vector, Hessian matrix with respect to to a list of variables.
        //         * @param var
        //         * @param gradient
        //         * @param hessian
        //         */
        //        void AddAdjointEntry(const atl::Variable<REAL_T>& var,
        //                const std::vector<atl::VariableInfo<REAL_T> >& variables,
        //                const std::vector<double>& gradient,
        //                const std::vector<const std::vector<double> >& hessian) {
        //            Adjoint<REAL_T>& adjoint = this->gradient_stack[this->NextIndex()];
        //            adjoint.w = var.info;
        //            adjoint.entries.resize(variables.size());
        //            adjoint.second_order_partials.resize(variables.size()*variables.size());
        //            for (int i = 0; i < variables.size(); i++) {
        //                adjoint.entries[i] = (std::move(AdjointDerivative<REAL_T>(variables[i].info, gradient[i])));
        //                for (int j = 0; j < variables.size(); j++) {
        //
        //                }
        //            }
        //
        //        }
        
        /**
         * Resets this stack and makes it available for a new recording.
         * 
         * @param empty_trash
         */
        inline void Reset(bool empty_trash = true) {
            
            if (max_initialized_size < stack_current) {
                max_initialized_size = stack_current;
            }
            
            if (this->recording) {
#pragma unroll
                for (int i = (stack_current - 1); i >= 0; i--) {
                    this->gradient_stack[i].Reset();
                }
                if (empty_trash) {
                    VariableInfo<REAL_T>::FreeAll();
                }
                stack_current = 0;
                
                gradient_computed = false;
            }
        }
        
        
    };
    
    
    
}

#endif	/* GRADIENTSTRUCTURE_HPP */
