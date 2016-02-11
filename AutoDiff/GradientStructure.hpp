/*
 * File:   GradientStructure.hpp
 * Author: matthewsupernaw
 *
 * Created on April 7, 2015, 4:05 PM
 */

#ifndef GRADIENTSTRUCTURE_HPP
#define GRADIENTSTRUCTURE_HPP
#include "Config.hpp"
#include <memory>
#include <vector>
#include <algorithm>
#include <functional>
#include <atomic>
#include <iostream>
#include "VariableInfo.hpp"
#include <fstream>
#include <cmath>
#include "AlignedAllocator.hpp"
#include "../Utilities/Combinations.hpp"
#include <unordered_set>
//#include <unordered_set>
//#define USE_BOOST

#include "DynamicExpression.hpp"


//#define ATL_USE_SMID

#ifdef ATL_USE_SMID
#include "../Utilities/SIMD.hpp"
#endif

#define ATL_ENABLE_BOUNDS_CHECKING
//#ifdef ATL_USE_SMID
//#include "../Utilities/SMID.hpp"
//#endif

#include "Variable.hpp"

//#define HESSIAN_TRACE

//#warning add jacobian matrix calculations

//#ifdef USE_BOOST
//#define IDSet boost::container::flat_set
//#else
//#define IDSet flat_set
////std::unordered_set
////flat_set
////
////flat_set
//#endif
#define Entry StackEntry<REAL_T>

//#define ATL_THREAD_SAFE


namespace atl {
    
    template<typename REAL_T>
    struct StackEntry {
        VariableInfo<REAL_T>* w; //function or dependent variable.
        atl::DynamicExpression<REAL_T>* exp;
        IDSet<atl::VariableInfo<REAL_T>* > ids;
        std::vector<VariableInfo<REAL_T>* > live_ids; //live variables used in reverse accumulation
        std::vector<atl::VariableInfo<REAL_T>* >id_list;
        std::vector<REAL_T> first;
        std::vector<REAL_T> second;
        std::vector<REAL_T> second_mixed;
        
        std::vector<REAL_T> third;
        std::vector<REAL_T> third_mixed;
        
        StackEntry() : w(NULL), exp(NULL) {
            first.reserve(10);
            ids.reserve(5);
        }
        
        StackEntry(const StackEntry<REAL_T>& orig) {
            this->w = orig.w;
            this->exp = orig.exp->Clone();
            typename IDSet<atl::VariableInfo<REAL_T>* >::const_iterator it;
            for (it = orig.ids.begin(); it != orig.ids.end(); ++it) {
                this->ids.insert((*it));
            }
            
            this->first.insert(this->first.begin(), orig.first.begin(), orig.first.end());
            this->second.insert(this->second.begin(), orig.second.begin(), orig.second.end());
        }
        
        inline void PushVariable(VariableInfo<REAL_T>* v) {
            if (v != w) {
                //                v->occurences--;
                //                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it =
                //                        ids.find(v);
                if (!ids.contains(v)) {
                    //                if (it == ids.end()) {
                    live_ids.push_back(v);
                }
            }
        }
        
        inline void PushVariables(const std::vector<VariableInfo<REAL_T>* >& v) {
            for (size_t i = 0; i < v.size(); i++) {
                if (v[i] != w) {
                    if (!ids.contains(v)) {
                        live_ids.push_back(v[i]);
                    }
                }
            }
        }
        
        inline void Prepare() {
            id_list.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator e;
            e = ids.end();
            for (it = ids.begin(); it != e; ++it) {
                id_list.push_back((*it));
            }
            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator ee;
            ee = live_ids.end();
            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator jt;
            for (jt = live_ids.begin(); jt != ee; ++jt) {
                id_list.push_back((*jt));
            }
            
        }
        
        //        inline void Prepare(atl::VariableInfo<REAL_T>* w) {
        //            id_list.resize(0);
        //            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
        //            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator e;
        //            e = ids.end();
        //            for (it = ids.begin(); it != e; ++it) {
        //                id_list.push_back((*it));
        //            }
        //
        //            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator jt;
        //            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator ee;
        //            ee = live_ids.end();
        //            for (jt = live_ids.begin(); jt != ee; ++jt) {
        //                REAL_T hij = 0.0;
        //                typename HessianInfo::iterator vijt;
        //                vijt = w->hessian_row.find((*jt)->id);
        //                if (vijt != w->hessian_row.end()) {
        //                    if ((*vijt) != 0) {
        //                        id_list.push_back((*jt));
        //                    }
        //                }
        //
        //            }
        //
        //        }
        
        inline void Reset() {
            
            first.resize(0);
            second_mixed.resize(0);
            third_mixed.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator end = ids.end();
            for (it = ids.begin(); it != end; ++it) {
                (*it)->Reset();
            }
            
            w->Reset();
            w = NULL;
            //            local_size = 0;
            live_ids.resize(0);
            if (exp) {
                delete exp;
                exp = NULL;
            }
            ids.clear();
            
        }
        
        //        inline void SoftReset() {
        //
        //
        //            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
        //            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator end = ids.end();
        //            for (it = ids.begin(); it != end; ++it) {
        //                (*it)->dvalue = 0;
        //            }
        //
        //            w->dvalue = 0;
        //
        //        }
        
    };
    
    enum DerivativeTraceLevel {
        FIRST_ORDER = 0, //SAME AS GRADIENT
        SECOND_ORDER, // SECOND ORDER PER VARIABLE ONLY
        THIRD_ORDER, // THIRD ORDER PER VARIABLE ONLY
        SECOND_ORDER_MIXED_PARTIALS, //SAME AS HESSIAN_AND_GRADIENT
        THIRD_ORDER_MIXED_PARTIALS,
        GRADIENT,
        GRADIENT_AND_HESSIAN,
        DYNAMIC_RECORD,
    };
    
    /**
     * Class to record operations. Often refered to as a "Tape". Holds a stack of
     * first and second order partial derivatives used in adjoint accumulation of
     * gradients and Hessian matrices.
     */
    template<typename REAL_T>
    class GradientStructure {
        typedef std::unordered_map<uint32_t, REAL_T> FirstOrder;
        typedef std::unordered_map<uint32_t, FirstOrder> SecondOrderMixed;
        typedef std::unordered_map<uint32_t, SecondOrderMixed> ThirdOrderMixed;
        
        FirstOrder first_order;
        SecondOrderMixed second_order_mixed;
        ThirdOrderMixed third_order_mixed;
        typedef typename FirstOrder::iterator fo_iterator;
        typedef typename SecondOrderMixed::iterator so_iterator;
        typedef typename ThirdOrderMixed::iterator to_iterator;
        
    public:
        DerivativeTraceLevel derivative_trace_level;
        std::vector<StackEntry<REAL_T> > gradient_stack;
#ifdef ATL_THREAD_SAFE
        std::atomic<size_t> stack_current;
#else
        size_t stack_current;
#endif
        bool recording;
        size_t max_stack_size;
        size_t max_initialized_size;
        
        bool gradient_computed;
        
        GradientStructure(uint32_t size = 1000000)
        : recording(true), stack_current(0),
        gradient_computed(false),
        derivative_trace_level(GRADIENT_AND_HESSIAN) {
            gradient_stack.resize(size);
            max_stack_size = size;
            max_initialized_size = 0;
        }
        
        GradientStructure(const GradientStructure<REAL_T>& other) :
        derivative_trace_level(other.derivative_trace_level),
        stack_current(other.stack_current),
        recording(other.recording),
        max_stack_size(other.max_stack_size),
        max_initialized_size(other.max_initialized_size),
        gradient_computed(other.gradient_computed) {
            
            //            this->gradient_stack.resize(other.stack_current);
            for (int i = 0; i < other.stack_current; i++) {
                this->gradient_stack.push_back(other.gradient_stack[i]);
            }
        }
        
        /**
         * Sets the size of the stack.
         * @param size
         */
        void SetSize(size_t size) {
            gradient_stack.resize(size);
        }
        
        virtual ~GradientStructure() {
            
        }
        
        inline void SetRecording(bool recording) {
            this->recording = recording;
        }
        
        //        inline REAL_T& Reference(uint32_t i) {
        //            return first_order[i];
        //        }
        
        inline REAL_T& Reference(uint32_t i, uint32_t j) {
            std::vector<uint32_t> indexes = {i, j};
            std::sort(indexes.begin(), indexes.end());
            return second_order_mixed[indexes[0]][indexes[1]];
        }
        
        inline REAL_T& Reference_No_Sort(uint32_t i, uint32_t j) {
            return second_order_mixed[i][j];
        }
        
        inline REAL_T& Reference(uint32_t i, uint32_t j, uint32_t k) {
            std::vector<uint32_t> indexes = {i, j, k};
            std::sort(indexes.begin(), indexes.end());
            return third_order_mixed[indexes[0]][indexes[1]][indexes[2]];
        }
        
        //        inline const REAL_T Value(uint32_t i) {
        //            fo_iterator it = first_order.find(i);
        //            if (it != first_order.end()) {
        //                return (*it).second;
        //            } else {
        //                return 0;
        //            }
        //            //            return this->first_order.get(i);
        //        }
        
        inline const REAL_T Value(uint32_t i, uint32_t j) {
            std::vector<uint32_t> indexes = {i, j};
            std::sort(indexes.begin(), indexes.end());
            so_iterator it = this->second_order_mixed.find(indexes[0]);
            
            if (it != this->second_order_mixed.end()) {
                fo_iterator itt = (*it).second.find(indexes[1]);
                if (itt != (*it).second.end()) {
                    return (*itt).second;
                } else {
                    return 0;
                }
            } else {
                return 0.0;
            }
            
            //            return this->second_order_mixed.get(indexes[0]).get(indexes[1]);
        }
        
        inline const REAL_T Value_No_Sort(uint32_t i, uint32_t j) {
            //            return this->second_order_mixed.get(i).get(j);
            so_iterator it = this->second_order_mixed.find(i);
            
            if (it != this->second_order_mixed.end()) {
                fo_iterator itt = (*it).second.find(j);
                if (itt != (*it).second.end()) {
                    return (*itt).second;
                } else {
                    return 0;
                }
            } else {
                return 0.0;
            }
        }
        
        inline const REAL_T Value(uint32_t i, uint32_t j, uint32_t k) {
            std::vector<uint32_t> indexes = {i, j, k};
            std::sort(indexes.begin(), indexes.end());
            
            to_iterator it = this->third_order_mixed.find(indexes[0]);
            
            if (it != this->third_order_mixed.end()) {
                so_iterator jt = (*it).second.find(indexes[1]);
                
                if (jt != (*it).second.end()) {
                    
                    fo_iterator kt = (*jt).second.find(indexes[2]);
                    if (kt != (*jt).second.end()) {
                        return (*kt).second;
                    } else {
                        return 0;
                    }
                    
                } else {
                    return 0.0;
                }
                
            } else {
                return 0.0;
            }
            
            //            return this->third_order_mixed.get(indexes[0]).get(indexes[1]).get(indexes[2]);
        }
        
        /**
         * Atomic operation. Gets the next available index in the stack.
         *
         * @return
         */
        inline const size_t NextIndex() {
#ifdef ATL_THREAD_SAFE
            return stack_current.fetch_add(1, std::memory_order_relaxed);
#else
            return stack_current++;
#endif
        }
        
        inline StackEntry<REAL_T>& NextEntry() {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            if ((this->stack_current) >= this->max_stack_size - 5) {
                std::cout << "Current derivative stack index exceeds stack limits.\n" << std::flush;
                
                gradient_stack.resize(this->max_stack_size * 2);
                this->max_stack_size = gradient_stack.size();
            }
#endif
            return this->gradient_stack[this->NextIndex()];
        }
        
        /**
         * Accumulates derivatives in reverse mode according to the member <i>derivative_trace</i>.
         *<br><br> <b>Reverse mode accumulation equations for each <i>derivative_trace</i> flag:</b><br>
         * <br><br>For <b><i>GRADIENT</i></b> or <b><i>FIRST_ORDER</i></b>
         * \image html gradient.png
         * For <b><i>HESSIAN</i></b> or <b><i>SECOND_ORDER_MIXED_PARTIALS</i></b>
         * \image html hessian.png
         * For <b><i>THIRD_ORDER_MIXED_PARTIALS</i></b>
         * \image html third_order.png
         *
         */
        inline void Accumulate() {
            gradient_computed = true;
            
            if (recording) {
                REAL_T w = 0.0;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator end;
                
                int j = 0;
                switch (this->derivative_trace_level) {
                        
                    case GRADIENT:
                        this->AccumulateFirstOrder();
                        //                        this->Reference(gradient_stack[stack_current - 1].w->id) = 1.0;
                        //#pragma unroll
                        //                        for (int i = (stack_current - 1); i >= 0; i--) {
                        //                            w = this->Value(gradient_stack[i].w->id); //gradient_stack[i].w->dvalue; //set w
                        //                            this->Reference(gradient_stack[i].w->id) = 0;
                        //                            if (w != static_cast<REAL_T> (0)) {
                        //                                j = 0;
                        //                                for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        //                                    REAL_T dv = w * gradient_stack[i].first[j];
                        //                                    //                            vj->dvalue += w * gradient_stack[i].first[j];
                        //                                    if (dv != 0.0) {
                        //                                        this->Reference((*it)->id) += dv;
                        //                                    }
                        //                                    j++;
                        //                                }
                        //                            }
                        //                        }
                        break;
                    case DYNAMIC_RECORD:
                        //                        gradient_stack[stack_current - 1].w->dvalue = 1.0;
                        //#pragma unroll
                        //                        for (int i = (stack_current - 1); i >= 0; i--) {
                        //                            w = this->gradient_stack[i].w->dvalue;
                        //                            gradient_stack[i].w->dvalue = seed;
                        //                            if (w != static_cast<REAL_T> (0)) {
                        //                                j = 0;
                        //                                for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        //                                    (*it)->dvalue += w * gradient_stack[i].exp->EvaluateDerivative((*it)->id);
                        //                                    j++;
                        //                                }
                        //                            }
                        //                        }
                        break;
                    case FIRST_ORDER:
                        this->AccumulateFirstOrder();
                        //
                        //                        this->Reference(gradient_stack[stack_current - 1].w->id) = 1.0;
                        //#pragma unroll
                        //                        for (int i = (stack_current - 1); i >= 0; i--) {
                        //                            w = this->Value(gradient_stack[i].w->id); //gradient_stack[i].w->dvalue; //set w
                        //                            this->Reference(gradient_stack[i].w->id) = 0;
                        //                            if (w != static_cast<REAL_T> (0)) {
                        //
                        //
                        //
                        //                                for (int j = 0; j < gradient_stack[i].ids.data().size(); j++) {
                        //
                        //                                }
                        //
                        //                                j = 0;
                        //
                        //                                for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        //                                    REAL_T dv = w * gradient_stack[i].first[j];
                        //                                    //                            vj->dvalue += w * gradient_stack[i].first[j];
                        //                                    if (dv != 0.0) {
                        //                                        this->Reference((*it)->id) += dv;
                        //                                    }
                        //                                    j++;
                        //                                }
                        //                            }
                        //
                        //                        }
                        break;
                        
                    case SECOND_ORDER:
                        //                        gradient_stack[stack_current - 1].w->dvalue = 1.0;
                        //#pragma unroll
                        //                        for (int i = (stack_current - 1); i >= 0; i--) {
                        //                            w = this->gradient_stack[i].w->dvalue;
                        //                            gradient_stack[i].w->dvalue = seed;
                        //                            if (w != static_cast<REAL_T> (0)) {
                        //                                j = 0;
                        //                                for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        //                                    (*it)->dvalue += w * gradient_stack[i].first[j];
                        //                                    (*it)->d2value += w * gradient_stack[i].second[j];
                        //                                    j++;
                        //                                }
                        //                            }
                        //                        }
                        break;
                    case THIRD_ORDER:
                        //                        gradient_stack[stack_current - 1].w->dvalue = 1.0;
                        //#pragma unroll
                        //                        for (int i = (stack_current - 1); i >= 0; i--) {
                        //                            w = this->gradient_stack[i].w->dvalue;
                        //                            gradient_stack[i].w->dvalue = seed;
                        //                            if (w != static_cast<REAL_T> (0)) {
                        //                                j = 0;
                        //                                for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        //                                    (*it)->dvalue += w * gradient_stack[i].first[j];
                        //                                    (*it)->d2value += w * gradient_stack[i].second[j];
                        //                                    (*it)->d3value += w * gradient_stack[i].third[j];
                        //                                    j++;
                        //                                }
                        //                            }
                        //                        }
                        break;
                    case GRADIENT_AND_HESSIAN:
                        this->AccumulateSecondOrderMixed();
                        //                        this->HessianAndGradientAccumulate();
                        break;
                    case SECOND_ORDER_MIXED_PARTIALS:
                        this->AccumulateSecondOrderMixed();
                        //                        this->HessianAndGradientAccumulate();
                        break;
                    case THIRD_ORDER_MIXED_PARTIALS:
                        this->AccumulateThirdOrderMixed();
                        //                        this->ThirdOrderMixedAccumulate();
                        break;
                    default:
                        std::cout << __func__ << "unknown trace level...\n";
                        //                        std::cout<<"\n";
                }
            }
        }
        
        void AccumulateFirstOrder() {
            
            
#ifdef ATL_USE_SMID
            REAL_T w;
            REAL_T adj[2];
            typedef typename simd::simd_traits<REAL_T>::type sse_t;
            sse_t sse_w;
            sse_t sse_d;
            sse_t sse_result;
            size_t sse_size = simd::simd_traits<REAL_T>::size;
            
#else
            REAL_T w = 0.0;
            int j = 0;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
#endif
            gradient_stack[stack_current - 1].w->dvalue = 1.0;
#pragma unroll
            for (int i = (stack_current - 1); i >= 0; i--) {
#ifdef ATL_USE_SMID
                
                w = gradient_stack[i].w->dvalue; //gradient_stack[i].w->dvalue; //set w
                //                w[1] = w[0];
                sse_w.load1(&w);
                
                
                
                if (w != static_cast<REAL_T> (0)) {
                    gradient_stack[i].w->dvalue = 0;
                    
                    int j = 0;
                    size_t size = gradient_stack[i].ids.data().size();
                    
                    for (; (j + sse_size) < size; j += sse_size) {
                        sse_d.load_u(&gradient_stack[i].first[j]);
                        sse_result = sse_d*sse_w;
                        sse_result.store_u(adj);
                        
                        if (adj[0] != static_cast<REAL_T> (0)) {
                            gradient_stack[i].ids.data()[j]->dvalue += adj[0];
                        }
                        
                        if (adj[1] != static_cast<REAL_T> (0)) {
                            gradient_stack[i].ids.data()[j + 1]->dvalue += adj[1];
                        }
                    }
                    
                    for (; j < size; j++) {
                        gradient_stack[i].ids.data()[j]->dvalue += w * gradient_stack[i].first[j];
                    }
                    
                }
#else
                w = gradient_stack[i].w->dvalue;
                
                if (w != static_cast<REAL_T> (0.0)) {
                    gradient_stack[i].w->dvalue = 0.0;
                    j = 0;
                    for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        (*it)->dvalue += w * gradient_stack[i].first[j];
                        j++;
                    }
                }
#endif
            }
        }
        
        /**
         * Computes the gradient and Hessian matrix via reverse mode accumulation.
         * \f[
         \begin{equation}
         \begin{split}
         \frac{\hat{\partial}}{\hat{\partial} v_b}\left[\frac{\hat{\partial} f_i}{\hat{\partial} v_c}\right] &=
         \frac{\partial^2 f_{i+1}}{\partial v_b \partial v_c} +
         \left(\frac{\partial^2 \phi_i}{\partial v_b \partial v_c} * \frac{\partial f_{i+1}}{\partial v_i} \right) +
         \left(\frac{\partial\phi_i}{\partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_b \partial v_i}\right)
         \\
         &+ \left(\frac{\partial\phi_i}{\partial v_b} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_c}\right) +
         \left(\frac{\partial\phi_i}{\partial v_b} * \frac{\partial\phi_i}{\partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_i}\right)
         \end{split}
         \end{equation}
         
         \begin{equation}
         \begin{split}
         \frac{\hat{\partial}}{\hat{\partial}  v_a}\left[\frac{\hat{\partial}}{\hat{\partial} v_b}\left(\frac{\hat{\partial} f_i}{\hat{\partial} v_c}\right)\right] &=
         \frac{\partial^3 f_{i+1}}{\partial v_a \partial v_b \partial v_c} +
         \left(\frac{\partial^3 \phi_i}{\partial v_a \partial v_b \partial v_c} * \frac{\partial f_{i+1}}{\partial v_i}\right) +
         \left(\frac{\partial^2 \phi_i}{\partial v_b \partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_a \partial v_i}\right)
         \\
         &+ \left(\frac{\partial^2 \phi_i}{\partial v_a \partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_b \partial v_i}\right) +
         \left(\frac{\partial \phi_i}{\partial v_c} * \frac{\partial^3 f_{i+1}}{\partial v_a \partial v_b \partial v_i}\right)
         \\
         &+ \left(\frac{\partial^2 \phi_i}{\partial v_a \partial v_b} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_c}\right) +
         \left(\frac{\partial \phi_i}{\partial v_b} * \frac{\partial^3 f_{i+1}}{\partial v_a \partial v_i \partial v_c}\right)
         \\
         &+ \left(\frac{\partial^2 \phi_i}{\partial v_a \partial v_b} * \frac{\partial \phi_i}{\partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_i}\right) +
         \left(\frac{\partial \phi_i}{\partial v_b} * \frac{\partial^2 \phi_i}{\partial v_a \partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_i}\right) +
         \\
         &+ \left(\frac{\partial \phi_i}{\partial v_b} * \frac{\partial \phi_i}{\partial v_c} * \frac{\partial^3 f_{i+1}}{\partial v_a \partial v_i \partial v_i}\right)
         \\
         &+ \frac{\partial \phi_i}{\partial v_a} *
         \left[\frac{\partial^3 f_{i+1}}{\partial v_i \partial v_b \partial v_c} +
         \left(\frac{\partial^3 \phi_i}{\partial v_i \partial v_b \partial v_c} * \frac{\partial f_{i+1}}{\partial v_i}\right) +
         \left(\frac{\partial^2 \phi_i}{\partial v_b \partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_i}\right)
         \right.
         \\
         &\left.
         + \left(\frac{\partial^2 \phi_i}{\partial v_i \partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_b \partial v_i}\right) +
         \left(\frac{\partial \phi_i}{\partial v_c} * \frac{\partial^3 f_{i+1}}{\partial v_i \partial v_b \partial v_i}\right)
         \right.
         \\
         &\left.
         + \left(\frac{\partial^2 \phi_i}{\partial v_i \partial v_b} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_c}\right) +
         \left(\frac{\partial \phi_i}{\partial v_b} * \frac{\partial^3 f_{i+1}}{\partial v_i \partial v_i \partial v_c}\right)
         \right.
         \\
         &\left.
         + \left(\frac{\partial^2 \phi_i}{\partial v_i \partial v_b} * \frac{\partial \phi_i}{\partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_i}\right) +
         \left(\frac{\partial \phi_i}{\partial v_b} * \frac{\partial^2 \phi_i}{\partial v_i \partial v_c} * \frac{\partial^2 f_{i+1}}{\partial v_i \partial v_i}\right)
         \right.
         \\
         &\left.
         + \left(\frac{\partial \phi_i}{\partial v_b} * \frac{\partial \phi_i}{\partial v_c} * \frac{\partial^3 f_{i+1}}{\partial v_i \partial v_i \partial v_i}\right)\right]
         \end{split}
         * ]
         *
         * \image html hessian.png
         */
        
        //        void HessianAndGradientAccumulate() {
        //
        //            if (recording) {
        //
        //
        //                if (this->derivative_trace_level == GRADIENT) {
        //                    this->Accumulate();
        //                } else {
        //                    REAL_T w;
        //                    REAL_T w2;
        //                    typedef typename VariableInfo<REAL_T>::HessianInfo HessianInfo;
        //
        //                    flat_set<VariableInfo<REAL_T>*> live;
        //
        //                    //initialize w
        //                    this->gradient_stack[stack_current - 1].w->dvalue = 1.0;
        //                    this->gradient_stack[stack_current - 1].w->d2value = 1.0;
        //
        //                    unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation
        //
        //                    std::vector<REAL_T> vij; //holds current second order derivative for i wrt j
        //
        //                    //Hessian iterators
        //                    typename HessianInfo::iterator vit;
        //                    typename HessianInfo::iterator vijt;
        //                    typename HessianInfo::iterator iend;
        //                    typename HessianInfo::iterator jend;
        //                    typename HessianInfo::iterator vjt;
        //
        //                    REAL_T hii = 0.0;
        //                    REAL_T hij = 0.0;
        //                    REAL_T hjk = 0;
        //                    REAL_T dj = 0;
        //                    REAL_T dk = 0;
        //
        //                    for (int i = (stack_current - 1); i >= 0; i--) {
        //
        //#ifdef HESSIAN_TRACE
        //                        std::cout << "\n\n\n";
        //                        std::stringstream ss;
        //                        int pushed_count = 0;
        //#endif
        //                        atl::VariableInfo<REAL_T>* vi = gradient_stack[i].w; //variable info for i
        //                        w = gradient_stack[i].w->dvalue; //set w
        //                        gradient_stack[i].w->dvalue = 0; //cancel out derivative for i
        //
        //                        rows = gradient_stack[i].first.size();
        //
        //                        //get h[i][i]
        //                        hii = 0.0;
        //
        //                        iend = vi->hessian_row.end();
        //                        vit = vi->hessian_row.find(vi->id);
        //                        if (vit != iend) {
        //                            hii = (*vit).second;
        //                            if (hii != REAL_T(0.0)) {
        //                                (*vit).second = 0.0;
        //                            }
        //                        }
        //
        //                        //prepare for the hessian calculation.
        //                        //builds a list of variables to use, statement level variables come first,
        //                        //then any pushed variables are after.
        //                        gradient_stack[i].Prepare();
        //
        //
        //
        //                        //resize second order derivative for i wrt j
        //                        vij.resize(gradient_stack[i].id_list.size());
        //
        //#pragma unroll
        //                        for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {
        //                            atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
        //
        //                            //compute gradient
        //                            if (j < rows && w != REAL_T(0.0)) {
        //                                vj->dvalue += w * gradient_stack[i].first[j];
        //                                vj->d2value += w * gradient_stack[i].second_mixed[j * rows + j];
        //                            }
        //
        //                            //load second order partial derivative for i wrt j and k
        //                            hij = 0.0;
        //
        //                            vijt = vi->hessian_row.find(vj->id);
        //                            if (vijt != iend) {
        //                                hij = (*vijt).second;
        //                                (*vijt).second = 0;
        //                            }
        //                            vij[j] = (hij);
        //
        //                        }
        //
        //                        std::vector<int> pushed_js(gradient_stack[i].id_list.size(), 0);
        //
        //#pragma unroll
        //                        //start the Hessian calculations
        //                        for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {
        //
        //                            atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j]; //vj
        //                            jend = vj->hessian_row.end();
        //                            REAL_T hij = vij[j]; //h[i][j]
        //
        //                            dj = 0;
        //                            bool j_in_local = false;
        //                            if (j < rows) {
        //                                dj = gradient_stack[i].first[j];
        //                                j_in_local = true;
        //                            }
        //
        //
        //                            //                            bool j_pushed = false;
        //#pragma unroll
        //                            //use symmetry
        //                            for (unsigned k = 0; k <= j; k++) {
        //                                REAL_T entry = 0.0; //the entry value for h[j][k]
        //                                atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k]; //vk
        //
        //                                dk = 0;
        //                                bool k_in_local = false;
        //                                if (k < rows) {
        //                                    dk = gradient_stack[i].first[k];
        //                                    k_in_local = true;
        //                                }
        //
        //                                if (!j_in_local && !k_in_local) {
        //                                    if (dk == REAL_T(0.0) && dj == REAL_T(0.0)) {
        //                                        continue; //pushed variable doesn't matter
        //                                    }
        //                                }
        //
        //
        //
        //                                entry += vij[k] * dj + (hij * dk) + hii * dj*dk;
        //
        //#ifdef HESSIAN_TRACE
        //                                REAL_T s = 0;
        //#endif
        //
        //                                if (j_in_local && k_in_local) {
        //                                    entry += w * gradient_stack[i].second_mixed[j * rows + k];
        //                                }
        //
        //#ifdef HESSIAN_TRACE
        //                                std::cout << "h[" << vj->id << "][" << vk->id << "] +=" << "h[" << vi->id << "][" << vj->id << "]{" << hij << "} *" << dk << " + " <<
        //                                        "h[" << vi->id << "][" << vk->id << "]{" << vij[k] << "}*" << dj << " + " <<
        //                                        "h[" << vi->id << "][" << vi->id << "]{" << hii << "} *" << dj << "*" << dk << " + " << w << "*" << s;
        //#endif
        //
        //                                if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
        //
        //#ifdef HESSIAN_TRACE
        //                                    std::cout << " = " << (entry + hjk) << "\n";
        //                                    if (j > this->gradient_stack[i].first.size()) {
        //                                        pushed_count++;
        //                                        std::cout << "Push mattered. occurences " << vj->occurences << "\n";
        //                                        ;
        //                                    }
        //#endif
        //                                    //set h[j][k]
        //                                    vj->hessian_row[vk->id] += entry; // + hjk;
        //                                    if (j > 0 && j != k) {
        //                                        //set h[k][j]
        //                                        vk->hessian_row[vj->id] += entry; /// + hjk;
        //                                    }
        //
        //                                    if (i > 0) {
        //                                        if (!pushed_js[j]) {
        //                                            //this variable may be needed in the future, so push it to the next entry
        //                                            gradient_stack[i - 1].PushVariable(vj);
        //                                            pushed_js[j] = 1;
        //                                            //                                            std::cout<<"pushed...\n";
        //                                            //                                            live.insert(vj);
        //#ifdef HESSIAN_TRACE
        //                                            ss << vj->id << " ";
        //#endif
        //                                            //                                            j_pushed = true;
        //                                        }
        //                                    }
        //                                }
        //#ifdef HESSIAN_TRACE
        //                                else {
        //                                    std::cout << " = 0\n";
        //                                }
        //#endif
        //                            }
        //                        }
        //
        //#ifdef HESSIAN_TRACE
        //                        std::cout << ss.str() << "\n";
        //#endif
        //                    }
        //                }
        //            }
        //        }
        
        void AccumulateSecondOrderMixed() {
            
            if (recording) {
                
                std::vector<std::vector<int> > combinations;
                
                REAL_T w;
                REAL_T w2;
                
                flat_set<VariableInfo<REAL_T>*> live;
                
                //initialize w
                this->gradient_stack[stack_current - 1].w->dvalue = 1.0;
                
                
                unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation
                
                std::vector<REAL_T> vij; //holds current second order derivative for i wrt j
                
                
                
                REAL_T hii = 0.0;
                REAL_T hij = 0.0;
                REAL_T hjk = 0;
                REAL_T dj = 0;
                REAL_T dk = 0;
                
                
                
                util::CombinationsWithRepetition combos(10, 2);
                
                for (int i = (stack_current - 1); i >= 0; i--) {
#ifdef HESSIAN_TRACE
                    std::cout << "\n\n\n";
                    std::stringstream ss;
                    int pushed_count = 0;
#endif
                    atl::VariableInfo<REAL_T>* vi = gradient_stack[i].w; //variable info for i
                    w = gradient_stack[i].w->dvalue; //gradient_stack[i].w->dvalue; //set w
                    gradient_stack[i].w->dvalue = 0; //cancel out derivative for i
                    
                    rows = gradient_stack[i].first.size();
                    
                    //get h[i][i]
                    hii = this->Value_No_Sort(vi->id, vi->id);
                    if (hii != 0) {
                        this->Reference(vi->id, vi->id) = 0.0;
                    }
                    
                    
                    //prepare for the hessian calculation.
                    //builds a list of variables to use, statement level variables come first,
                    //then any pushed variables are after.
                    gradient_stack[i].Prepare();
                    
                    
                    
                    //resize second order derivative for i wrt j
                    vij.resize(gradient_stack[i].id_list.size());
                    
#pragma unroll
                    for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {
                        //                        std::cout << "push " << j << std::endl;
                        atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                        
                        //compute gradient
                        if (j < rows && w != REAL_T(0.0)) {
                            vj->dvalue += w * gradient_stack[i].first[j];
                            //                            vj->dvalue += w * gradient_stack[i].first[j];
                            //                            if (dv != 0.0) {
                            //                                this->Reference(vj->id) += dv;
                            //                            }
                        }
                        
                        //load second order partial derivative for i wrt j and k
                        hij = this->Value(vi->id, vj->id);
                        
                        if (hij != 0) {
                            this->Reference(vi->id, vj->id) = 0.0;
                        }
                        
                        
                        vij[j] = (hij);
                        
                    }
                    
                    
                    
                    std::vector<int> pushed_js(gradient_stack[i].id_list.size(), 0);
                    
                    combos.Reset(gradient_stack[i].id_list.size(), 2);
                    
                    do {
                        int j = combos[0];
                        int k = combos[1];
                        
                        atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                        atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k];
                        
                        
                        REAL_T hij = vij[j]; //h[i][j]
                        
                        dj = 0;
                        //                        bool j_in_local = false;
                        if (j < rows) {
                            dj = gradient_stack[i].first[j];
                            //                            j_in_local = true;
                        }
                        
                        REAL_T entry = 0.0; //the entry value for h[j][k]
                        
                        dk = 0;
                        if (k < rows) {
                            dk = gradient_stack[i].first[k];
                        }
                        
                        if ((j >= rows) && (k >= rows)) {
                            if (dk == REAL_T(0.0) && dj == REAL_T(0.0)) {
                                continue; //pushed variable doesn't matter
                            }
                        }
                        
                        
                        
                        entry += vij[k] * dj + (hij * dk) + hii * dj*dk;
                        
                        
                        //                            REAL_T s = 0;
                        if (j < rows && k < rows) {
                            entry += w * gradient_stack[i].second_mixed[j * rows + k];
                            //                                s = gradient_stack[i].second_mixed[j * rows + k];
                        }
                        //                            std::cout << "h[" << vj->id << "][" << vk->id << "] +=" << "h[" << vi->id << "][" << vj->id << "]{" << hij << "} *" << dk << " + " <<
                        //                                    "h[" << vi->id << "][" << vk->id << "]{" << vij[k] << "}*" << dj << " + " <<
                        //                                    "h[" << vi->id << "][" << vi->id << "]{" << hii << "} *" << dj << "*" << dk << " + " << w << "*" << s << " = " << entry << "\n";
                        //
                        
                        
                        if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                            
                            this->Reference(vj->id, vk->id) += entry;
                            
                            if (i > 0) {
                                if (!pushed_js[j]) {
                                    //this variable may be needed in the future, so push it to the next entry
                                    gradient_stack[i - 1].PushVariable(vj);
                                    //                                        std::cout << "Pushing " << vj->id << "\n";
                                    pushed_js[j] = 1;
                                }
                                if (!pushed_js[k]) {
                                    //this variable may be needed in the future, so push it to the next entry
                                    gradient_stack[i - 1].PushVariable(vk);
                                    //                                        std::cout << "Pushing " << vk->id << "\n";
                                    pushed_js[k] = 1;
                                }
                            }
                        }
                        
                        
                        
                        
                    } while (combos.Next());
                    //                    }
                }
            }
        }
        
        void AccumulateThirdOrderMixed() {
            
            if (recording) {
                
                std::vector<std::vector<int> > combinations;
                
                REAL_T w;
                
                //initialize w
                this->gradient_stack[stack_current - 1].w->dvalue = 1.0;
                unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation
                
                std::vector<REAL_T> vij; //holds current second order derivative for i wrt j
                std::vector<REAL_T> viij_;
                std::vector<REAL_T> vijk_;
                
                
                REAL_T hii = 0.0;
                REAL_T hij = 0.0;
                //                REAL_T hik = 0.0;
                //                REAL_T hjk = 0;
                REAL_T diii = 0.0;
                REAL_T dijl = 0.0;
                REAL_T dikl = 0.0;
                //                REAL_T djkl = 0.0;
                REAL_T dijk = 0.0;
                REAL_T diil = 0.0;
                REAL_T diik = 0.0;
                REAL_T diij = 0.0;
                REAL_T dj = 0.0;
                REAL_T dk = 0.0;
                REAL_T dl = 0.0;
                
                util::CombinationsWithRepetition combos(10, 3);
                
                for (int i = (stack_current - 1); i >= 0; i--) {
                    atl::VariableInfo<REAL_T>* vi = gradient_stack[i].w; //variable info for i
                    
                    w = this->gradient_stack[i].w->dvalue; //set w
                    this->gradient_stack[i].w->dvalue = 0; //cancel out derivative for i
                    
                    //                        std::cout<<"w = "<<w<<"\n";
                    rows = gradient_stack[i].first.size();
                    
                    //get h[i][i]
                    hii = this->Value(vi->id, vi->id);
                    
                    if (hii != 0.0) {
                        this->Reference(vi->id, vi->id) = 0.0;
                    }
                    
                    diii = this->Value(vi->id, vi->id, vi->id);
                    
                    if (diii != 0.0) {
                        this->Reference(vi->id, vi->id, vi->id) = 0.0;
                    }
                    
                    
                    //prepare for the hessian calculation.
                    //builds a list of variables to use, statement level variables come first,
                    //then any pushed "live" variables are after.
                    gradient_stack[i].Prepare();
                    
                    
                    //resize second order derivative for i wrt j
                    vij.resize(gradient_stack[i].id_list.size());
                    viij_.resize(gradient_stack[i].id_list.size());
                    vijk_.resize(gradient_stack[i].id_list.size() * gradient_stack[i].id_list.size());
#pragma unroll
                    for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {
                        atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                        
                        //compute gradient
                        if (j < rows && w != REAL_T(0.0)) {
                            vj->dvalue += w * gradient_stack[i].first[j];
                        }
                        
                        
                        
                        //load second order partial derivative for i wrt j and k
                        hij = this->Value(vi->id, vj->id);
                        
                        if (hij != 0.0) {
                            this->Reference(vi->id, vj->id) = 0.0;
                        }
                        
                        
                        vij[j] = hij;
                        
                        
                        
                        diil = this->Value(vi->id, vi->id, vj->id);
                        
                        if (diil != 0.0) {
                            this->Reference(vi->id, vi->id, vj->id) = 0.0;
                        }
                        
                        viij_[j] = diil;
                        
                        for (unsigned k = j; k < gradient_stack[i].id_list.size(); k++) {
                            
                            dijk = 0.0;
                            atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k];
                            
                            dijk = this->Value(vi->id, vj->id, vk->id);
                            
                            if (dijk != 0.0) {
                                this->Reference(vi->id, vj->id, vk->id) = 0.0;
                            }
                            
                            vijk_[(j * gradient_stack[i].id_list.size()) + k] = dijk;
                            vijk_[(k * gradient_stack[i].id_list.size()) + j] = dijk;
                        }
                    }
                    
                    
                    std::vector<int> pushed_js(gradient_stack[i].id_list.size(), 0);
                    
                    
                    combos.Reset(gradient_stack[i].id_list.size(), 3);
                    do {
                        
                        int j = combos[0]; //combinations[combos][0];
                        int k = combos[1]; //combinations[combos][1];
                        int l = combos[2]; //combinations[combos][2];
                        
                        
                        
                        atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                        atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k];
                        atl::VariableInfo<REAL_T>* vl = gradient_stack[i].id_list[l];
                        
                        dj = 0;
                        if (j < rows) {
                            dj = gradient_stack[i].first[j];
                        }
                        
                        REAL_T entry = 0.0; //the entry value for h[j][k]
                        
                        dk = 0;
                        if (k < rows) {
                            dk = gradient_stack[i].first[k];
                        }
                        
                        
                        if (j == 0) {
                            
                            REAL_T hdk = 0;
                            REAL_T hdj = 0;
                            
                            if (k < rows) {
                                hdj = gradient_stack[i].first[k];
                            }
                            
                            REAL_T entry = 0.0; //the entry value for h[j][k]
                            
                            if (l < rows) {
                                hdk = gradient_stack[i].first[l];
                            }
                            
                            
                            
                            entry += vij[l] * hdj + (vij[k] * hdk) + hii * hdj*hdk;
                            
                            
                            if (l < rows && k < rows) {
                                entry += w * gradient_stack[i].second_mixed[k * rows + l];
                            }
                            
                            
                            if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                                
                                //                            std::cout << "entry = " << entry << "\n";
                                this->Reference(vk->id, vl->id) += entry;
                                
                                if (i > 0) {
                                    //this variable may be needed in the future, so push it to the next entry
                                    if (!pushed_js[k]) {
                                        gradient_stack[i - 1].PushVariable(vk);
                                        pushed_js[k] = 1;
                                    }
                                    if (!pushed_js[l]) {
                                        gradient_stack[i - 1].PushVariable(vl);
                                        pushed_js[l] = 1;
                                    }
                                    //                                    if (!pushed_js[j]) {
                                    //                                        gradient_stack[i - 1].PushVariable(vk);
                                    //                                        pushed_js[j] = 1;
                                    //                                    }
                                    //                                    if (!pushed_js[k]) {
                                    //                                        gradient_stack[i - 1].PushVariable(vl);
                                    //                                        pushed_js[k] = 1;
                                    //                                    }
                                    
                                }
                            }
                        }
                        
                        REAL_T hij = vij[j]; //h[i][j]
                        REAL_T hik = vij[k];
                        REAL_T hil = vij[l];
                        diij = viij_[j];
                        diik = viij_[k];
                        diil = viij_[l];
                        
                        
                        //find dijk
                        dijk = vijk_[(j * gradient_stack[i].id_list.size() + k)];
                        dijl = vijk_[(j * gradient_stack[i].id_list.size() + l)];
                        dikl = vijk_[(k * gradient_stack[i].id_list.size() + l)];
                        
                        
                        
                        REAL_T d3 = 0.0;
                        
                        dl = 0.0;
                        
                        //
                        if (l < rows) {
                            dl = gradient_stack[i].first[l];
                            if (k < rows && j < rows) {
                                d3 = gradient_stack[i].third_mixed[(j * rows * rows) + (k * rows) + l];
                            }
                        }
                        
                        REAL_T pjk = 0.0;
                        REAL_T pjl = 0.0;
                        REAL_T pkl = 0.0;
                        
                        
                        if (k < rows && l < rows) {
                            pkl = gradient_stack[i].second_mixed[k * rows + l];
                        }
                        
                        if (j < rows && l < rows) {
                            pjl = gradient_stack[i].second_mixed[j * rows + l];
                        }
                        
                        if (j < rows && k < rows) {
                            pjk = gradient_stack[i].second_mixed[j * rows + k];
                        }
                        
                        
                        
                        
                        REAL_T entry_3 = 0;
                        entry_3 = (d3 * w) + (pkl * hij)
                        + (pjl * hik) + (dl * dijk)
                        +(pjk * hil) + (dk * dijl)
                        +(pjk * dl * hii) + (dk * pjl * hii)+(dk * dl * diij)
                        + dj * (dikl + (pkl * hii)+(dl * diik) + (dk * diil)
                                +(dk * dl * diii));
                        
                        //                        std::cout << "dj = " << dj << "," << vj->dvalue << "\n";
                        //                        std::cout << "dk = " << dk << "," << vk->dvalue << "\n";
                        //                        std::cout << "dl = " << dl << "," << vl->dvalue << "\n";
                        //                        std::cout << "diii = " << diii << "\n";
                        //                        std::cout << "dikl = " << dikl << "\n";
                        //                        std::cout << "dijl = " << dijl << "\n";
                        //                        std::cout << "dijk = " << dijk << "\n";
                        //                        std::cout << "diij = " << diij << "\n";
                        //                        std::cout << "d3 = " << d3 << "\n";
                        //                        std::cout << "pjk = " << pjk << "\n";
                        //                        std::cout << "pjl = " << pjl << "\n";
                        //                        std::cout << "pkl = " << pjl << "\n";
                        //                        std::cout << "hii = " << hii << "\n";
                        //                        std::cout << "hij = " << hij << "\n";
                        //                        std::cout << "hik = " << hik << "\n";
                        //                        std::cout << "hil = " << hil << "\n";
                        //                        std::cout << "w = " << w << "\n";
                        //
                        //
                        //
                        //                        std::cout << "entry_3 = " << (d3 * w) << " + " << (pkl * hij)
                        //                                << " + " << (pjl * hik) << " + " << (dl * dijk)
                        //                                << " + " << (pjk * hil) << " + " << (dk * dijl)
                        //                                << " + " << (pjk * dl * hii) << " + " << (dk * pjl * hii) << " + " << (dk * dl * diij)
                        //                                << " + " << dj * (dikl + (pkl * hii)+(dl * diik) + (dk * diil)
                        //                                +(dk * dl * diii));
                        //                        //
                        //                        std::cout << " = " << entry_3 << "\n";
                        
                        if (entry_3 != 0.0) {
                            this->Reference(vj->id, vk->id, vl->id) += entry_3;
                            if (i > 0) {
                                if (!pushed_js[j]) {
                                    gradient_stack[i - 1].PushVariable(vj);
                                    pushed_js[j] = 1;
                                }
                                if (!pushed_js[k]) {
                                    gradient_stack[i - 1].PushVariable(vk);
                                    pushed_js[k] = 1;
                                }
                                if (!pushed_js[l]) {
                                    gradient_stack[i - 1].PushVariable(vl);
                                    pushed_js[l] = 1;
                                }
                            }
                        }
                    } while (combos.Next());
                }
            }
        }
        
        void AccumulateThirdOrderMixed2() {
            
            if (recording) {
                
                
                
                REAL_T w;
                
                //initialize w
                this->Reference(this->gradient_stack[stack_current - 1].w->id) = 1.0;
                unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation
                
                std::vector<REAL_T> vij; //holds current second order derivative for i wrt j
                std::vector<REAL_T> viij_;
                std::vector<REAL_T> vijk_;
                
                
                REAL_T hii = 0.0;
                REAL_T hij = 0.0;
                REAL_T hik = 0.0;
                REAL_T hjk = 0;
                REAL_T diii = 0.0;
                REAL_T dijl = 0.0;
                REAL_T dikl = 0.0;
                REAL_T djkl = 0.0;
                REAL_T dijk = 0.0;
                REAL_T diil = 0.0;
                REAL_T diik = 0.0;
                REAL_T diij = 0.0;
                REAL_T dj = 0.0;
                REAL_T dk = 0.0;
                REAL_T dl = 0.0;
                
                std::vector<int> indexes;
                for (int i = (stack_current - 1); i >= 0; i--) {
                    atl::VariableInfo<REAL_T>* vi = gradient_stack[i].w; //variable info for i
                    
                    w = this->Value(this->gradient_stack[stack_current - 1].w->id); //set w
                    this->Reference(this->gradient_stack[stack_current - 1].w->id) = 0; //cancel out derivative for i
                    
                    //                        std::cout<<"w = "<<w<<"\n";
                    rows = gradient_stack[i].first.size();
                    
                    //get h[i][i]
                    hii = this->Value(vi->id, vi->id);
                    
                    if (hii != 0.0) {
                        this->Reference(vi->id, vi->id) = 0.0;
                    }
                    
                    diii = this->Value(vi->id, vi->id, vi->id);
                    
                    if (diii != 0.0) {
                        this->Reference(vi->id, vi->id, vi->id) = 0.0;
                    }
                    
                    
                    //prepare for the hessian calculation.
                    //builds a list of variables to use, statement level variables come first,
                    //then any pushed "live" variables are after.
                    gradient_stack[i].Prepare();
                    
                    
                    //resize second order derivative for i wrt j
                    vij.resize(gradient_stack[i].id_list.size());
                    viij_.resize(gradient_stack[i].id_list.size());
                    vijk_.resize(gradient_stack[i].id_list.size() * gradient_stack[i].id_list.size());
                    indexes.resize(0);
#pragma unroll
                    for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {
                        indexes.push_back(j);
                        atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                        
                        //compute gradient
                        REAL_T dv = w * gradient_stack[i].first[j];
                        //                            vj->dvalue += w * gradient_stack[i].first[j];
                        if (dv != 0.0) {
                            this->Reference(vj->id) += dv;
                        }
                        
                        //load second order partial derivative for i wrt j and k
                        hij = this->Value(vi->id, vj->id);
                        
                        if (hij != 0.0) {
                            this->Reference(vi->id, vj->id) = 0.0;
                        }
                        
                        
                        vij[j] = (hij);
                        
                        diil = this->Value(vi->id, vi->id, vj->id);
                        
                        if (diil != 0.0) {
                            this->Reference(vi->id, vi->id, vj->id) = 0.0;
                        }
                        
                        viij_[j] = diil;
                        
                        for (unsigned k = 0; k < gradient_stack[i].id_list.size(); k++) {
                            
                            dijk = 0.0;
                            atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k];
                            
                            dijk = this->Value(vi->id, vj->id, vk->id);
                            
                            if (dijk != 0.0) {
                                this->Reference(vi->id, vj->id, vk->id) = 0.0;
                            }
                            
                            vijk_[(j * gradient_stack[i].id_list.size()) + k] = dijk;
                            
                        }
                    }
                    
                    
                    std::vector<int> pushed_js(gradient_stack[i].id_list.size(), 0);
                    
                    
                    if (indexes.size() == 1) {
                        
                        int j = 0;
                        int k = 0;
                        int l = 0;
                        
                        atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                        atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k];
                        atl::VariableInfo<REAL_T>* vl = gradient_stack[i].id_list[l];
                        
                        
                        
                        dj = 0;
                        bool j_in_local = false;
                        if (j < rows) {
                            dj = gradient_stack[i].first[j];
                            j_in_local = true;
                        }
                        
                        REAL_T entry = 0.0; //the entry value for h[j][k]
                        
                        dk = 0;
                        bool k_in_local = false;
                        if (k < rows) {
                            dk = gradient_stack[i].first[k];
                            k_in_local = true;
                        }
                        
                        if (!j_in_local && !k_in_local) {
                            if (dk == REAL_T(0.0) && dj == REAL_T(0.0)) {
                                continue; //pushed variable doesn't matter
                            }
                        }
                        
                        
                        
                        
                        entry += vij[k] * dj + (vij[j] * dk) + hii * dj*dk;
                        
                        
                        //                        REAL_T s = 0;
                        if (j_in_local && k_in_local) {
                            entry += w * gradient_stack[i].second_mixed[j * rows + k];
                        }
                        
                        
                        if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                            
                            //                            std::cout << "entry = " << entry << "\n";
                            this->Reference(vj->id, vk->id) += entry;
                            
                            if (i > 0) {
                                //this variable may be needed in the future, so push it to the next entry
                                gradient_stack[i - 1].PushVariable(vj);
                                gradient_stack[i - 1].PushVariable(vk);
                            }
                        }
                        
                        diik = viij_[k];
                        
                        
                        //find dijk
                        dijk = vijk_[(j * gradient_stack[i].id_list.size() + k)];
                        
                        REAL_T hil = vij[l];
                        
                        
                        dijl = vijk_[(j * gradient_stack[i].id_list.size() + l)];
                        
                        dikl = vijk_[(k * gradient_stack[i].id_list.size() + l)];
                        
                        diil = viij_[l];
                        
                        
                        REAL_T d3 = 0.0;
                        
                        dl = 0.0;
                        
                        //
                        if (l < rows) {
                            dl = gradient_stack[i].first[l];
                            if (k < rows && j < rows) {
                                d3 = gradient_stack[i].third_mixed[(j * rows * rows) + (k * rows) + l];
                            }
                        }
                        
                        
                        REAL_T sj = 0.0;
                        REAL_T pjk = 0.0;
                        REAL_T pjl = 0.0;
                        REAL_T pkl = 0.0;
                        REAL_T sk = 0.0;
                        REAL_T sl = 0.0;
                        
                        
                        
                        if (k < rows && l < rows) {
                            pkl = gradient_stack[i].second_mixed[k * rows + l];
                            sj = dj * gradient_stack[i].second_mixed[k * rows + l];
                        }
                        
                        if (j < rows && l < rows) {
                            pjl = gradient_stack[i].second_mixed[j * rows + l];
                            sk = dk * gradient_stack[i].second_mixed[j * rows + l];
                        }
                        
                        if (j < rows && k < rows) {
                            pjk = gradient_stack[i].second_mixed[j * rows + k];
                            sl = dl * gradient_stack[i].second_mixed[j * rows + k];
                        }
                        
                        
                        
                        
                        REAL_T entry_3 = 0;
                        entry_3 = (d3 * w) + (pkl * hij)
                        + (pjl * hik) + (dl * dijk)
                        +(pjk * hil) + (dk * dijl)
                        +(pjk * dl * hii) + (dk * pjl * hii)+(dk * dl * diij)
                        + dj * (dikl + (pkl * hii)+(dl * diik) + (dk * diil)
                                +(dk * dl * diii));
                        
                        if (entry_3 != 0.0) {
                            this->Reference(vj->id, vk->id, vl->id) += entry_3;
                            if (i > 0) {
                                gradient_stack[i - 1].PushVariable(vj);
                                gradient_stack[i - 1].PushVariable(vk);
                                gradient_stack[i - 1].PushVariable(vl);
                            }
                        }
                        
                        
                        
                    } else {
                        
                        util::Combonations combinations(indexes, 3);
                        
                        for (int combos = 0; combos < combinations.GetCombinations().size(); combos++) {
                            
                            
                            int j = combinations.GetCombinations()[combos][0];
                            int k = combinations.GetCombinations()[combos][1];
                            int l = combinations.GetCombinations()[combos][2];
                            
                            
                            
                            atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                            atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k];
                            atl::VariableInfo<REAL_T>* vl = gradient_stack[i].id_list[l];
                            
                            
                            
                            dj = 0;
                            bool j_in_local = false;
                            if (j < rows) {
                                dj = gradient_stack[i].first[j];
                                j_in_local = true;
                            }
                            
                            REAL_T entry = 0.0; //the entry value for h[j][k]
                            
                            dk = 0;
                            bool k_in_local = false;
                            if (k < rows) {
                                dk = gradient_stack[i].first[k];
                                k_in_local = true;
                            }
                            
                            if (!j_in_local && !k_in_local) {
                                if (dk == REAL_T(0.0) && dj == REAL_T(0.0)) {
                                    continue; //pushed variable doesn't matter
                                }
                            }
                            
                            
                            
                            
                            entry += vij[k] * dj + (vij[j] * dk) + hii * dj*dk;
                            
                            
                            //                        REAL_T s = 0;
                            if (j_in_local && k_in_local) {
                                entry += w * gradient_stack[i].second_mixed[j * rows + k];
                            }
                            
                            
                            if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                                
                                //                            std::cout << "entry = " << entry << "\n";
                                this->Reference(vj->id, vk->id) += entry;
                                
                                if (i > 0) {
                                    //this variable may be needed in the future, so push it to the next entry
                                    gradient_stack[i - 1].PushVariable(vj);
                                    gradient_stack[i - 1].PushVariable(vk);
                                }
                            }
                            
                            diik = viij_[k];
                            
                            
                            //find dijk
                            dijk = vijk_[(j * gradient_stack[i].id_list.size() + k)];
                            
                            REAL_T hil = vij[l];
                            
                            
                            dijl = vijk_[(j * gradient_stack[i].id_list.size() + l)];
                            
                            dikl = vijk_[(k * gradient_stack[i].id_list.size() + l)];
                            
                            diil = viij_[l];
                            
                            
                            REAL_T d3 = 0.0;
                            dl = 0.0;
                            
                            //                                
                            if (l < rows) {
                                dl = gradient_stack[i].first[l];
                                if (k < rows && j < rows) {
                                    d3 = gradient_stack[i].third_mixed[(j * rows * rows) + (k * rows) + l];
                                }
                            }
                            
                            
                            REAL_T sj = 0.0;
                            REAL_T pjk = 0.0;
                            REAL_T pjl = 0.0;
                            REAL_T pkl = 0.0;
                            REAL_T sk = 0.0;
                            REAL_T sl = 0.0;
                            
                            
                            
                            if (k < rows && l < rows) {
                                pkl = gradient_stack[i].second_mixed[k * rows + l];
                                //                                sj = dj * gradient_stack[i].second_mixed[k * rows + l];
                            }
                            
                            if (j < rows && l < rows) {
                                pjl = gradient_stack[i].second_mixed[j * rows + l];
                                //                                sk = dk * gradient_stack[i].second_mixed[j * rows + l];
                            }
                            
                            if (j < rows && k < rows) {
                                pjk = gradient_stack[i].second_mixed[j * rows + k];
                                //                                sl = dl * gradient_stack[i].second_mixed[j * rows + k];
                            }
                            
                            
                            
                            
                            REAL_T entry_3 = 0;
                            entry_3 = (d3 * w) + (pkl * hij)
                            + (pjl * hik) + (dl * dijk)
                            +(pjk * hil) + (dk * dijl)
                            +(pjk * dl * hii) + (dk * pjl * hii)+(dk * dl * diij)
                            + dj * (dikl + (pkl * hii)+(dl * diik) + (dk * diil)
                                    +(dk * dl * diii));
                            
                            if (entry_3 != 0.0) {
                                this->Reference(vj->id, vk->id, vl->id) += entry_3;
                                if (i > 0) {
                                    gradient_stack[i - 1].PushVariable(vj);
                                    gradient_stack[i - 1].PushVariable(vk);
                                    gradient_stack[i - 1].PushVariable(vl);
                                    
                                }
                            }
                            
                            
                        }
                        
                    }
                    
                }
            }
        }
        
        /**
         * Resets this stack and makes it available for a new recording.
         *
         * @param empty_trash
         */
        inline void Reset(bool empty_trash = true) {
            
            this->first_order.clear();
            this->second_order_mixed.clear();
            this->third_order_mixed.clear();
            
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
        
        inline void SoftReset() {
            this->first_order.clear();
            this->second_order_mixed.clear();
            this->third_order_mixed.clear();
            
            for (int i = (stack_current - 1); i >= 0; i--) {
                this->gradient_stack[i].SoftReset();
            }
        }
        
        
        
    };
    
    
    
}

#endif /* GRADIENTSTRUCTURE_HPP */
