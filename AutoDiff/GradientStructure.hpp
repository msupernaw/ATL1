/*
 * File:   GradientStructure.hpp
 * Author: matthewsupernaw
 *
 * Created on April 7, 2015, 4:05 PM
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


#ifndef GRADIENTSTRUCTURE_HPP
#define GRADIENTSTRUCTURE_HPP
#define ATL_THREAD_SAFE

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
#include "../Utilities/Combinations.hpp"
#include "../Utilities/flat_map.hpp"
#include "DynamicExpression.hpp"

#ifdef ATL_USE_SMID
#include "../Utilities/SIMD.hpp"
#endif

#define ATL_ENABLE_BOUNDS_CHECKING


#include "Variable.hpp"


#define Entry StackEntry<REAL_T>



namespace atl {

    template<typename REAL_T>
    struct StackEntry {
        VariableInfo<REAL_T>* w; //function or dependent variable.
        atl::DynamicExpression<REAL_T>* exp;
        IDSet<atl::VariableInfo<REAL_T>* > ids;
        typedef typename IDSet<atl::VariableInfo<REAL_T>* >::iterator id_itereator;
        IDSet<VariableInfo<REAL_T>* > live_ids; //live variables used in reverse accumulation
        std::vector<atl::VariableInfo<REAL_T>* > id_list;
        std::vector<atl::VariableInfo<REAL_T>* > valid_id_list;
        std::vector<REAL_T> first;
        std::vector<REAL_T> second;
        std::vector<REAL_T> second_mixed;
        std::vector<REAL_T> third;
        std::vector<REAL_T> third_mixed;
        uint32_t max_id = std::numeric_limits<uint32_t>::min();
        uint32_t min_id = std::numeric_limits<uint32_t>::max();

        StackEntry() : w(NULL), exp(NULL) {

        }

        StackEntry(const StackEntry<REAL_T>& other) :
        w(other.w), exp(other.exp), ids(other.ids), live_ids(other.live_ids), id_list(other.id_list), first(other.first), second(other.second), second_mixed(other.second_mixed), third(other.third), third_mixed(other.third_mixed), max_id(other.max_id), min_id(other.min_id) {
        }
        //        StackEntry(const StackEntry<REAL_T>& orig) {
        //            this->w = orig.w;
        //            
        //            if (orig.exp)
        //                this->exp = orig.exp->Clone();
        //            typename IDSet<atl::VariableInfo<REAL_T>* >::const_iterator it;
        //            for (it = orig.ids.begin(); it != orig.ids.end(); ++it) {
        //                this->ids.insert((*it));
        //            }
        //
        //            this->first.insert(this->first.begin(), orig.first.begin(), orig.first.end());
        //            this->second.insert(this->second.begin(), orig.second.begin(), orig.second.end());
        //        }

        const std::pair<uint32_t, uint32_t> FindMinMax() {
            std::pair<uint32_t, uint32_t> p;
            p.first = this->w->id;
            p.second = this->w->id;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            for (it = this->ids.begin(); it != ids.end(); ++it) {
                if ((*it)->id < p.first) {
                    p.first = (*it)->id;
                }
                if ((*it)->id > p.second) {
                    p.second = (*it)->id;
                }
            }
            return p;
        }

        inline void PushVariable(VariableInfo<REAL_T>* v) {
            if (v != w) {

                id_itereator it = ids.find(v);

                if (it == ids.end()) {
                    live_ids.insert(v);
                }
            }
        }

        inline void PushVariables(const std::vector<VariableInfo<REAL_T>* >& v) {
            for (size_t i = 0; i < v.size(); i++) {
                if (v[i] != w && v[i]->is_dependent == 0) {
                    id_itereator it = ids.find(v[i]);

                    if (it == ids.end()) {
                        live_ids.insert(v[i]);
                    }
                }
            }
        }

        inline void Prepare() {
            id_list.resize(0);
            valid_id_list.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator dt;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator e;
            e = ids.end();
            for (it = ids.begin(); it != e; ++it) {
                id_list.push_back((*it));
                valid_id_list.push_back((*it));
            }

            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator ee;
            ee = live_ids.end();
            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator jt;

            for (jt = live_ids.begin(); jt != ee; ++jt) {

                valid_id_list.push_back((*jt));
                id_list.push_back((*jt));
            }

            for (jt = w->nldependencies.begin(); jt != w->nldependencies.end(); ++jt) {
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator dt;
                dt = ids.find((*jt));
                if (dt == ids.end()) {
                    dt = live_ids.find((*jt));
                    if (dt == live_ids.end()) {
                        if ((*jt) != this->w) {
                            valid_id_list.push_back((*jt));
                            id_list.push_back((*jt));
                        }
                    }
                }


            }


        }

        inline void IntermediateReset() {
            this->live_ids.clear();
        }

        inline void SoftReset() {
            first.resize(0);
            second_mixed.resize(0);
            third_mixed.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator end = ids.end();
            for (it = ids.begin(); it != end; ++it) {
                (*it)->Reset();
            }

            w->Reset();
            this->live_ids.clear();
        }

        inline void Reset() {

            max_id = std::numeric_limits<uint32_t>::min();
            min_id = std::numeric_limits<uint32_t>::max();
            first.resize(0);
            second_mixed.resize(0);
            third_mixed.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator end = ids.end();
            //if using the memory pool, these maybe deleted before reset is called
            for (it = ids.begin(); it != end; ++it) {
                if ((*it) != NULL) {
                    (*it)->Reset();
                }
            }
            if (w != NULL) {
                w->Reset();
                w = NULL;
            }
            this->live_ids.clear();

            //            if (exp) {
            //                delete exp;
            //                exp = NULL;
            //            }
            ids.clear_no_resize();

        }


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

    template<class T>
    class DerivativeMatrix {
        size_t rows;
        size_t columns;
        std::vector<flat_map<size_t, size_t> > row_indices;

        std::vector<T> data;


    public:

        DerivativeMatrix(size_t rows = 1, size_t columns = 1) :
        rows(rows), columns(columns) {
            row_indices.resize(rows);
            data.reserve(rows);

        }

        inline T& operator()(const size_t& r, const size_t& c) {

            flat_map<size_t, size_t>::iterator it = row_indices[r].find(c);

            if (it != row_indices[r].end()) {
                return data[(*it).second];
            } else {
                data.push_back(T());
                row_indices[r][c] = data.size() - 1;
                return data[data.size() - 1];
            }

        }

        inline const T value(const size_t& r, const size_t& c) {
            flat_map<size_t, size_t>::iterator it = row_indices[r].find(c);
            return it != row_indices[r].end() ? data[(*it).second] : T();
        }

        inline void remove(const size_t& r, const size_t& c) {
            flat_map<size_t, size_t>::iterator it = row_indices[r].find(c);
            data[(*it).second] = 0.0;
            row_indices[r].erase(c);
        }

        void zero() {
            std::fill(data.begin(), data.end(), static_cast<T> (0.0));
        }

        void resize(size_t r, size_t c) {
            row_indices.resize(r);
            rows = r;
            columns = c;
        }
    };


    template<typename T>
    class GradientStructure;

    template<typename T>
    void PrepareThread(const int& start, const int& end, std::vector<T>& vij, std::vector<T>& viij_,
            std::vector<T>& vijk_, const size_t& ID_LIST_SIZE, atl::VariableInfo<T>* vi,
            StackEntry<T>& entry, atl::GradientStructure<T>* gs) {

        atl::VariableInfo<T>* vk;
        atl::VariableInfo<T>* vj;
        atl::VariableInfo<T>* vl;
        for (unsigned j = start; j < end; j++) {

            vj = entry.id_list[j];

            //load second order partial derivative for i wrt j and k
            vij[j] = gs->Value(vi->id, vj->id);

            if (std::fabs(vij[j]) > 0.0) {
                gs->Reference(vi->id, vj->id) = 0.0;
            }

            viij_[j] = gs->Value(vi->id, vi->id, vj->id);

            if (std::fabs(viij_[j]) > 0.0) {
                gs->Reference(vi->id, vi->id, vj->id) = 0.0;
            }


#pragma unroll
            for (unsigned k = j; k < ID_LIST_SIZE; k++) {
                vk = entry.id_list[k];

                vijk_[(j * ID_LIST_SIZE) + k] = gs->Value(vi->id, vj->id, vk->id);

                if (vijk_[(j * ID_LIST_SIZE) + k] != 0.0) {
                    gs->Reference(vi->id, vj->id, vk->id) = 0.0;

                }
                vijk_[(k * ID_LIST_SIZE) + j] = vijk_[(j * ID_LIST_SIZE) + k]; // dijk;
            }
        }


    }

    /**
     * Class to record operations. Often refered to as a "Tape". Holds a stack of
     * first and second order partial derivatives used in adjoint accumulation of
     * gradients and Hessian matrices.
     */
    template<typename REAL_T>
    class GradientStructure {
        typedef typename flat_map<uint32_t, REAL_T>::iterator derivative_iterator;

        flat_map<uint32_t, flat_map<uint32_t, REAL_T> > second;
        typedef typename flat_map<uint32_t, flat_map<uint32_t, REAL_T> >::iterator second_order_iterator;

        flat_map<uint32_t, flat_map<uint32_t, flat_map<uint32_t, REAL_T> > > third;
        typedef typename flat_map<uint32_t, flat_map<uint32_t, flat_map<uint32_t, REAL_T> > >::iterator third_order_iterator;

    public:
        uint32_t max_id = std::numeric_limits<uint32_t>::min();
        uint32_t min_id = std::numeric_limits<uint32_t>::max();
        size_t range;
        DerivativeTraceLevel derivative_trace_level;
        std::vector<StackEntry<REAL_T> > gradient_stack;
        std::map<uint32_t, StackEntry<REAL_T> > initialized_variables;
        typedef typename std::map<uint32_t, StackEntry<REAL_T> >::iterator initialized_variables_iterator;
#ifdef ATL_THREAD_SAFE
        std::mutex stack_lock;
#endif
        std::atomic<size_t> stack_current;
        size_t stack_begin;

        bool recording;
        size_t max_stack_size;
        size_t max_initialized_size;

        bool gradient_computed;
//        std::mutex stack_lock;

        GradientStructure(uint32_t size = 10000)
        : recording(true), stack_current(0), stack_begin(0),
        gradient_computed(false),
        derivative_trace_level(FIRST_ORDER) {
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

            for (int i = 0; i < other.stack_current; i++) {
                this->gradient_stack.push_back(other.gradient_stack[i]);
            }
        }

        void Begin() {
            this->stack_begin = stack_current + 1;
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

        inline REAL_T& Reference(uint32_t i) {
            //            return first_order[i];
        }

        inline REAL_T& Reference(uint32_t i, uint32_t j) {

            if (j < i) {
                std::swap(i, j);
            }
            return this->second[i][j];
        }

        inline REAL_T& ReferenceNoSort(uint32_t i, uint32_t j) {


            return this->second[i][j];
        }

        inline void MakeZero(uint32_t i, uint32_t j) {
            if (j < i) {
                std::swap(i, j);
            }
            this->second[i][j] = 0.0;
        }

        inline void MakeZeroNoSort(uint32_t i, uint32_t j) {
            this->second[i][j] = 0.0;
        }

        inline REAL_T& Reference(uint32_t i, uint32_t j, uint32_t k) {
            if (i < j) {
                if (k < i) std::swap(i, k);
            } else {
                if (j < k) std::swap(i, j);
                else std::swap(i, k);
            }
            if (k < j) std::swap(j, k);
            return third[i][j][k];
        }

        inline REAL_T& ReferenceNoSort(uint32_t i, uint32_t j, uint32_t k) {
            return third[i][j][k];
        }

        inline void MakeZero(uint32_t i, uint32_t j, uint32_t k) {
            if (i < j) {
                if (k < i) std::swap(i, k);
            } else {
                if (j < k) std::swap(i, j);
                else std::swap(i, k);
            }
            if (k < j) std::swap(j, k);
            third[i][j][k] = 0.0;
        }

        inline void MakeZeroNoSort(uint32_t i, uint32_t j, uint32_t k) {
            third[i][j][k] = 0.0;
        }

        inline const REAL_T Value(uint32_t i) {
//                        return this->first_order.get(i);
        }

        inline const REAL_T Value(uint32_t i, uint32_t j) {
            if (j < i) {
                std::swap(i, j);
            }

            second_order_iterator it = second.find(i);
            if (it != second.end()) {
                derivative_iterator jt = (*it).second.find(j);
                if (jt != (*it).second.end()) {
                    return (*jt).second;
                } else {
                    return static_cast<REAL_T> (0.0);
                }

            }
            return 0.0;

        }

        inline const REAL_T ValueNoSort(uint32_t i, uint32_t j) {

            second_order_iterator it = second.find(i);
            if (it != second.end()) {
                derivative_iterator jt = (*it).second.find(j);
                if (jt != (*it).second.end()) {
                    return (*jt).second;
                } else {
                    return static_cast<REAL_T> (0.0);
                }

            }
            return 0.0;

        }

        inline const REAL_T Value(uint32_t i, uint32_t j, uint32_t k) {

            if (i < j) {
                if (k < i) std::swap(i, k);
            } else {
                if (j < k) std::swap(i, j);
                else std::swap(i, k);
            }

            if (k < j) std::swap(j, k);

            third_order_iterator it = third.find(i);
            if (it != third.end()) {
                second_order_iterator jt = (*it).second.find(j);
                if (jt != (*it).second.end()) {
                    derivative_iterator kt = (*jt).second.find(k);
                    if (kt != (*jt).second.end()) {
                        return (*kt).second;
                    } else {
                        return static_cast<REAL_T> (0.0);
                    }
                } else {
                    return static_cast<REAL_T> (0.0);
                }
            }

            return static_cast<REAL_T> (0.0);

        }

        inline const REAL_T ValueNoSort(uint32_t i, uint32_t j, uint32_t k) {

            third_order_iterator it = third.find(i);
            if (it != third.end()) {
                second_order_iterator jt = (*it).second.find(j);
                if (jt != (*it).second.end()) {
                    derivative_iterator kt = (*jt).second.find(k);
                    if (kt != (*jt).second.end()) {
                        return (*kt).second;
                    } else {
                        return static_cast<REAL_T> (0.0);
                    }
                } else {
                    return static_cast<REAL_T> (0.0);
                }
            }

            return static_cast<REAL_T> (0.0);

        }

        /**
         * Atomic operation. Gets the next available index in the stack.
         *
         * @return
         */
        inline const size_t NextIndex() {
#ifdef ATL_THREAD_SAFE
            stack_lock.lock();
#endif
            if (stack_current + 1 >= this->gradient_stack.size()) {
                //                std::cout<<"Resizing Tape structure...\n";
                this->gradient_stack.resize(this->gradient_stack.size() + 100);
            }

#ifdef ATL_THREAD_SAFE
            stack_lock.unlock();
#endif       
            return this->stack_current++;
        }

        inline StackEntry<REAL_T>& NextEntry() {
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
                        break;
                    case DYNAMIC_RECORD:
                        break;
                    case FIRST_ORDER:
                        this->AccumulateFirstOrder();
                        break;

                    case SECOND_ORDER:
                        std::cout << "SECOND_ORDER not yet implemented!";
                        exit(0);
                        break;
                    case THIRD_ORDER:
                        std::cout << "THIRD_ORDER not yet implemented!";
                        exit(0);
                        break;
                    case GRADIENT_AND_HESSIAN:
                        this->AccumulateSecondOrderMixed();
                        break;
                    case SECOND_ORDER_MIXED_PARTIALS:
                        this->AccumulateSecondOrderMixed();
                        break;
                    case THIRD_ORDER_MIXED_PARTIALS:
                        this->AccumulateThirdOrderMixed();
                        break;
                    default:
                        std::cout << __func__ << "unknown trace level...\n";
                }
            }
        }

        void AccumulateFirstOrder() {

            this->PrepareDerivativeTables(1);
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
                //                std::cout<<"I = "<<i<<"\n"<<gradient_stack[i].w->dependence_level<<"\n";
                if (w != static_cast<REAL_T> (0.0)) {
                    gradient_stack[i].w->dvalue = 0.0;
                    j = 0;
                    for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        //                        std::cout<< (*it)->dvalue<<"+="<<w<<"*"<<gradient_stack[i].first[j]<<"\n";
                        (*it)->dvalue += w * gradient_stack[i].first[j];

                        j++;
                    }
                }
#endif
            }
        }

        void AccumulateFirstOrderDynamic() {


            REAL_T w = 0.0;
            int j = 0;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;

            gradient_stack[stack_current - 1].w->dvalue = 1.0;
#pragma unroll
            for (int i = (stack_current - 1); i >= 0; i--) {
                gradient_stack[i].w->vvalue = gradient_stack[i].exp->Evaluate();
                w = gradient_stack[i].w->dvalue;

                if (w != static_cast<REAL_T> (0.0)) {
                    gradient_stack[i].w->dvalue = 0.0;
                    j = 0;
                    for (it = gradient_stack[i].ids.begin(); it != gradient_stack[i].ids.end(); ++it) {
                        (*it)->dvalue += w * gradient_stack[i].exp->EvaluateDerivative((*it)->id);
                        j++;
                    }
                }

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

        void AccumulateSecondOrderMixed() {

            if (recording) {
                //                this->PrepareDerivativeTables(2);           

                REAL_T w;
                REAL_T w2;


                //initialize w
                this->gradient_stack[stack_current - 1].w->dvalue = 1.0;


                unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation

                std::vector<REAL_T> vij; //holds current second order derivative for i wrt j



                REAL_T hii = 0.0;
                REAL_T hij = 0.0;
                REAL_T hjk = 0;
                REAL_T dj = 0;
                REAL_T dk = 0;

                atl::VariableInfo<REAL_T>* vi;
                atl::VariableInfo<REAL_T>* vj;
                atl::VariableInfo<REAL_T>* vk;

                for (int i = (stack_current - 1); i >= 0; i--) {
#ifdef HESSIAN_TRACE
                    std::cout << "\n\n\n";
                    std::stringstream ss;
                    int pushed_count = 0;
#endif
                    vi = gradient_stack[i].w; //variable info for i
                    w = gradient_stack[i].w->dvalue; //gradient_stack[i].w->dvalue; //set w
                    gradient_stack[i].w->dvalue = 0; //cancel out derivative for i

                    rows = gradient_stack[i].first.size();

                    //get h[i][i]
                    hii = 0.0;
                    if (vi->so_mark) {
                        hii = this->Value(vi->id, vi->id);
                        if (hii != 0) {
                            this->MakeZero(vi->id, vi->id); // = 0.0;
                        }
                    }

                    //prepare for the hessian calculation. 
                    //builds a list of variables to use, statement level variables come first,
                    //then any pushed variables are after.
                    gradient_stack[i].Prepare();

                    size_t ID_LIST_SIZE = gradient_stack[i].id_list.size();

                    //resize second order derivative for i wrt j
                    vij.resize(ID_LIST_SIZE);
                    std::vector<bool> needs_push(ID_LIST_SIZE, false);



                    if (w != 0.0) {
                        for (unsigned j = 0; j < rows; j++) {
                            atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];
                            vj->dvalue += w * gradient_stack[i].first[j];
                        }
                    }

#pragma unroll
                    for (unsigned j = 0; j < ID_LIST_SIZE; j++) {
                        //                        std::cout << "push " << j << std::endl;
                        vj = gradient_stack[i].id_list[j];

                        //load second order partial derivative for i wrt j and k
                        hij = 0.0;
                        if (vi->so_mark && vj->so_mark) {
                            hij = this->Value(vi->id, vj->id);

                            if (hij != 0) {
                                this->MakeZero(vi->id, vj->id); // = 0.0;
                            }
                        }
                        vij[j] = (hij);
                    }

                    REAL_T entry;
#pragma unroll
                    for (int j = 0; j < rows; j++) {
                        vj = gradient_stack[i].id_list[j];
                        dj = gradient_stack[i].first[j];
                        REAL_T hij = vij[j]; //h[i][j]
#pragma unroll
                        for (int k = j; k < rows; k++) {

                            vk = gradient_stack[i].id_list[k];

                            entry = 0.0; //the entry value for h[j][k]


                            dk = gradient_stack[i].first[k];






                            entry += vij[k] * dj + (hij * dk) + hii * dj*dk;

                            //                            std::cout<<"w = "<<w<<std::endl;
                            //                            std::cout<<"i = "<<i<<" of "<<stack_current<<" rows = "<<rows<<" "<<gradient_stack[i].second_mixed.size()<<" "<<j<<"-"<<k<<" "<<(j*rows+k)<<std::endl;
                            entry += w * gradient_stack[i].second_mixed[j * rows + k];


                            //                            if (gradient_stack[i].second_mixed[j * rows + k] != 0.0) {
                            //                                vj->push_count = 1;
                            //                                vk->push_count = 1;
                            //                            }



                            if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                                if (entry != entry) {
                                    std::cout << "Derivative signaling NaN\n";
                                    //                                    exit(0);
                                }
                                this->Reference(vj->id, vk->id) += entry;
                                vj->so_mark = true;
                                vk->so_mark = true;
                                needs_push[k] = true;
                                //                                needs_push[j] = true;
                            }

                        }

                        for (int k = rows; k < ID_LIST_SIZE; k++) {

                            vk = gradient_stack[i].id_list[k];

                            entry = 0.0; //the entry value for h[j][k]

                            dk = 0;





                            entry += vij[k] * dj; // + (hij * dk) + hii * dj*dk;


                            if (j < rows && k < rows) {

                                entry += w * gradient_stack[i].second_mixed[j * rows + k];
                                //                                if (gradient_stack[i].second_mixed[j * rows + k] != 0.0) {
                                //                                    vj->push_count = 1;
                                //                                    vk->push_count = 1;
                                //                                }
                            }


                            if (/*std::fabs(entry)*/entry != REAL_T(0.0) && entry == entry) {//h[j][k] needs to be updated
                                if (entry != entry) {
                                    std::cout << "Derivative signaling NaN\n";
                                    //                                    exit(0);
                                }
                                this->Reference(vj->id, vk->id) += entry;
                                vj->so_mark = true;
                                vk->so_mark = true;
                                needs_push[k] = true;
                                //                                needs_push[j] = true;
                            }

                        }
                    }
                    if (gradient_stack[i].w->dependence_level > 0) {//this was a compound assignment and its dependencies must be pushed
                        if (i > 0) {
#pragma unroll
                            for (int ii = 0; ii < rows; ii++) {
                                gradient_stack[i - 1].PushVariable(gradient_stack[i].id_list[ii]);

                            }
#pragma unroll
                            for (int ii = rows; ii < ID_LIST_SIZE; ii++) {

                                if ((!gradient_stack[i].id_list[ii]->is_dependent) && needs_push[ii]) {
                                    gradient_stack[i - 1].PushVariable(gradient_stack[i].id_list[ii]);
                                }

                            }


                        }
                        gradient_stack[i].w->dependence_level--;
                    }
                }
            }
        }

        void AccumulateThirdOrderMixed() {


            if (recording) {

                //                this->PrepareDerivativeTables(3);


                REAL_T w;

                //initialize w
                gradient_stack[stack_current - 1].w->dvalue = 1.0;
                unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation

                std::vector<REAL_T> vij; //holds current second order derivative for i wrt j
                std::vector<REAL_T> viij_;
                std::vector<REAL_T> vijk_;


                REAL_T hii = 0.0;
                REAL_T diii = 0.0;
                REAL_T dijk = 0.0;
                REAL_T diil = 0.0;
                REAL_T dj = 0.0;
                REAL_T dk = 0.0;
                REAL_T dl = 0.0;
                REAL_T pjl = 0.0;
                REAL_T pkl = 0.0;
                REAL_T entry_3 = 0;
                REAL_T d3 = 0.0;
                REAL_T hij = 0.0;
                REAL_T pjk = 0.0;




                for (int i = (stack_current - 1); i >= 0; i--) {

                    atl::VariableInfo<REAL_T>* vi = gradient_stack[i].w; //variable info for i
                    w = gradient_stack[i].w->dvalue; //set w
                    gradient_stack[i].w->dvalue = 0; //cancel out derivative for i

                    rows = gradient_stack[i].first.size();

                    hii = 0.0;

                    //get h[i][i]
                    if (vi->so_mark) {
                        hii = Value(vi->id, vi->id);

                        if (hii != 0.0) {
                            Reference(vi->id, vi->id) = 0.0;
                        }
                    }

                    diii = 0.0;
                    if (vi->to_mark) {
                        diii = Value(vi->id, vi->id, vi->id);

                        if (diii != 0.0) {
                            Reference(vi->id, vi->id, vi->id) = 0.0;
                        }
                    }

                    //prepare for the hessian calculation. 
                    //builds a list of variables to use, statement level variables come first,
                    //then any pushed "independent" variables are after.
                    gradient_stack[i].Prepare();

                    size_t ID_LIST_SIZE = gradient_stack[i].valid_id_list.size();
                    //                    std::cout << "Rows = " << rows << "\n";
                    //                    std::cout << "ID_LIST_SIZE = " << ID_LIST_SIZE << "\n";
                    //resize second order derivative for i wrt j
                    vij.resize(ID_LIST_SIZE);
                    viij_.resize(ID_LIST_SIZE);
                    vijk_.resize(ID_LIST_SIZE * ID_LIST_SIZE);

                    std::vector<bool> needs_push(ID_LIST_SIZE, false);


                    std::map<uint32_t, IDSet<atl::VariableInfo<REAL_T>* > > mattered_with;

                    //compute gradient
                    if (w != REAL_T(0.0)) {
                        for (unsigned j = 0; j < rows; j++) {
                            gradient_stack[i].valid_id_list[j]->dvalue += w * gradient_stack[i].first[j];
                        }
                    }
                    atl::VariableInfo<REAL_T>* vk;
                    atl::VariableInfo<REAL_T>* vj;
                    atl::VariableInfo<REAL_T>* vl;
                    //prepare higher order stuff
#pragma unroll

                    for (unsigned j = 0; j < ID_LIST_SIZE; j++) {

                        vj = gradient_stack[i].valid_id_list[j];

                        //load second order partial derivative for i wrt j and k
                        vij[j] = 0.0;
                        if (vi->so_mark && vj->so_mark) {
                            vij[j] = Value(vi->id, vj->id);

                            if (std::fabs(vij[j]) > 0.0) {
                                MakeZero(vi->id, vj->id);
                            }
                        }

                        viij_[j] = 0.0;
                        if (vi->to_mark && vj->to_mark) {
                            viij_[j] = Value(vi->id, vi->id, vj->id);

                            if (std::fabs(viij_[j]) > 0.0) {
                                MakeZero(vi->id, vi->id, vj->id);
                            }
                        }

#pragma unroll
                        for (unsigned k = j; k < ID_LIST_SIZE; k++) {
                            vk = gradient_stack[i].valid_id_list[k];

                            vijk_[(j * ID_LIST_SIZE) + k] = 0.0;

                            if (vi->to_mark && vj->to_mark && vk->to_mark) {
                                vijk_[(j * ID_LIST_SIZE) + k] = Value(vi->id, vj->id, vk->id);

                                if (vijk_[(j * ID_LIST_SIZE) + k] != 0.0) {

                                    MakeZero(vi->id, vj->id, vk->id);

                                }
                            }
                            vijk_[(k * ID_LIST_SIZE) + j] = vijk_[(j * ID_LIST_SIZE) + k]; // dijk;
                        }
                    }


                    REAL_T entry;
                    REAL_T hdk;
                    REAL_T hdj;
                    int mcount = 0;
#pragma unroll
                    for (int j = 0; j < rows; j++) {
                        vj = gradient_stack[i].valid_id_list[j];
                        dj = gradient_stack[i].first[j];

                        if (j == 0) {
#pragma unroll
                            for (int k = j; k < rows; k++) {
                                hdj = 0;
                                hdj = gradient_stack[i].first[k];
                                atl::VariableInfo<REAL_T>* vk = gradient_stack[i].valid_id_list[k];
#pragma unroll

                                for (int l = k; l < rows; l++) {

                                    vl = gradient_stack[i].valid_id_list[l];
                                    hdk = 0;

                                    entry = 0.0; //the entry value for h[j][k]

                                    hdk = gradient_stack[i].first[l];


                                    entry += vij[l] * hdj + (vij[k] * hdk) + hii * hdj*hdk;


                                    entry += w * gradient_stack[i].second_mixed[k * rows + l];


                                    if (entry != entry) {
                                        std::cout << "Derivative signaling NaN\n";
                                        //                                        exit(0);
                                    }
                                    if (/*std::fabs(entry)*/entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                                        this->Reference(vk->id, vl->id) += entry;
                                        vk->so_mark = true;
                                        vl->so_mark = true;
                                        needs_push[l] = true;
                                    }
                                }

                                for (int l = rows; l < ID_LIST_SIZE; l++) {

                                    vl = gradient_stack[i].valid_id_list[l];
                                    
                                    hdk = 0;

                                    entry = 0.0; //the entry value for h[j][k]

                                    entry += vij[l] * hdj; // + (vij[k] * hdk) + hii * hdj*hdk;

                                    if (entry_3 != entry_3) {
                                        std::cout << "Derivative signaling NaN\n";
                                        //                                        exit(0);
                                    }
                                    if (entry != REAL_T(0.0)) {//h[j][k] needs to be updated
                                        Reference(vk->id, vl->id) += entry;
                                        vk->so_mark = true;
                                        vl->so_mark = true;
                                        needs_push[l] = true;
                                        needs_push[k] = true;


                                    }
                                }
                            }
                        }

#pragma unroll
                        for (int k = j; k < rows; k++) {
                            vk = gradient_stack[i].valid_id_list[k];

                            dk = gradient_stack[i].first[k];
                            pjk = gradient_stack[i].second_mixed[j * rows + k];

                            for (int l = k; l < rows; l++) {
                                vl = gradient_stack[i].valid_id_list[l];
                                entry_3 = 0;

                                dl = gradient_stack[i].first[l];
                                pjl = gradient_stack[i].second_mixed[j * rows + l];
                                pkl = gradient_stack[i].second_mixed[k * rows + l];

                                d3 = gradient_stack[i].third_mixed[(j * rows * rows) + (k * rows) + l];

                                entry_3 += (d3 * w)
                                        +(pjl * vij[k])
                                        + (dk * pjl * hii)
                                        + (dl * vijk_[(j * ID_LIST_SIZE + k)])
                                        + (pkl * vij[j])
                                        + (dk * dl * viij_[j])
                                        + (pjk * dl * hii);


                                entry_3 += (pjk * vij[l])
                                        + (dk * vijk_[(j * ID_LIST_SIZE + l)]);





                                entry_3 += dj * (vijk_[(k * ID_LIST_SIZE + l)] + (pkl * hii)+(dl * viij_[k]) + (dk * viij_[l])
                                        +(dk * dl * diii));

                                if (entry_3 != 0.0) {
                                    Reference(vj->id, vk->id, vl->id) += entry_3;
                                    vj->to_mark = true;
                                    vk->to_mark = true;
                                    vl->to_mark = true;
                                }

                            }

                            for (int l = rows; l < ID_LIST_SIZE; l++) {

                                vl = gradient_stack[i].valid_id_list[l];
                                entry_3 = 0;

                                dl = 0.0;
                                pjl = 0.0;
                                pkl = 0.0;


                                entry_3 += (pjk * vij[l])
                                        + (dk * vijk_[(j * ID_LIST_SIZE + l)]);





                                entry_3 += dj * (vijk_[(k * ID_LIST_SIZE + l)] + (dk * viij_[l]));

                                if (entry_3 != entry_3) {
                                    std::cout << "Derivative signaling NaN\n";
                                }
                                if (entry_3 != 0.0) {
                                    Reference(vj->id, vk->id, vl->id) += entry_3;
                                    vj->to_mark = true;
                                    vk->to_mark = true;
                                    vl->to_mark = true;
                                    needs_push[l] = true;
                                    needs_push[k] = true;
                                }

                            }
                        }

                        if (dj != 0.0) {
#pragma unroll
                            for (int k = rows; k < ID_LIST_SIZE; k++) {
                                atl::VariableInfo<REAL_T>* vk = gradient_stack[i].valid_id_list[k];
#pragma unroll
                                for (int l = k; l < ID_LIST_SIZE; l++) {
                                    atl::VariableInfo<REAL_T>* vl = gradient_stack[i].valid_id_list[l];
                                    entry_3 = dj * (vijk_[(k * ID_LIST_SIZE + l)]);

                                    if (entry_3 != entry_3) {
                                        std::cout << "Derivative signaling NaN\n";
                                    }
                                    if (entry_3 != 0.0) {
                                        Reference(vj->id, vk->id, vl->id) += entry_3;
                                        vj->to_mark = true;
                                        vk->to_mark = true;
                                        vl->to_mark = true;
                                        needs_push[l] = true;
                                        needs_push[k] = true;
                                    }

                                }
                            }
                        }
                    }





                    if (gradient_stack[i].w->dependence_level > 0) {//this was a compound assignment and its dependencies must be pushed
                        if (i > 0) {

                            for (int ii = 0; ii < rows; ii++) {
                                gradient_stack[i - 1].PushVariable(gradient_stack[i].id_list[ii]);

                            }

                            for (int ii = rows; ii < ID_LIST_SIZE; ii++) {
                                if ((!gradient_stack[i].id_list[ii]->is_dependent) && needs_push[ii]) {
                                    gradient_stack[i - 1].PushVariable(gradient_stack[i].id_list[ii]);
                                }

                            }

                        }
                        gradient_stack[i].w->dependence_level--;
                    }
                }


            }

        }

        void AccumulateThirdOrderMixedDynamic() {

            std::cout << __func__ << "Not yet implemented!\n";
            exit(0);

        }

        /**
         * Resets this stack and makes it available for a new recording.
         *
         * @param empty_trash
         */
        inline void Reset(bool empty_trash = true) {
            max_id = std::numeric_limits<uint32_t>::min();
            min_id = std::numeric_limits<uint32_t>::max();

            if (max_initialized_size < stack_current) {
                max_initialized_size = stack_current;
            }
            if (this->second.size()) {
                this->second.clear();
                this->third.clear();
            }
#pragma unroll
            for (int i = (stack_current - 1); i >= 0; i--) {
                this->gradient_stack[i].Reset();
            }
            //
            //            initialized_variables_iterator it;
            //            for (it = this->initialized_variables.begin(); it != this->initialized_variables.end();) {
            //                if ((*it).second.w->count == 0) {
            //                    this->initialized_variables.erase(it++);
            //                } else {
            //                    ++it;
            //                }
            //            }

            if (empty_trash) {
                VariableInfo<REAL_T>::FreeAll();
            }
            stack_current = 0;

            gradient_computed = false;
        }

        void PrepareDerivativeTables(int level = 2) {
            if (level > 1) {



                max_id = 0;
                min_id = std::numeric_limits<uint32_t>::max();

                for (int i = 0; i < this->stack_current; i++) {
                    std::pair<size_t, size_t> p = this->gradient_stack[i].FindMinMax();
                    if (p.first< this->min_id) {
                        this->min_id = p.first;
                    }

                    if (p.second > this->max_id) {
                        this->max_id = p.second;
                    }
                }
                //                range = (max_id - min_id) + 1;
                //                if (level == 2) {
                ////                    somp.zero();
                ////                    this->somp.resize(range + 1, range + 1);
                ////                    for (int i = 0; i < tomp.size(); i++) {
                ////                        tomp[i].zero();
                ////                    }
                //                }
                //                if (level == 3) {
                //                    this->second.clear();
                //                    this->third.clear();
                //                    //                    somp.zero();
                //                    //                    this->somp.resize(range + 1, range + 1);
                //                    //                    this->tomp.resize(range + 1);
                //                    //                    for (int i = 0; i < tomp.size(); i++) {
                //                    //                        tomp[i].zero();
                //                    //                        tomp[i].resize(range + 1, range + 1);
                //                    //                    }
                //                }
            }

        }

        inline void SoftReset() {

            //            somp.zero();
            //            this->somp.resize(range + 1, range + 1);
            //            this->tomp.resize(range + 1);
            //            for (int i = 0; i < tomp.size(); i++) {
            //                tomp[i].zero();
            //                tomp[i].resize(range + 1, range + 1);
            //            }
            //
            this->second.clear();
            this->third.clear();
            for (int i = (stack_current - 1); i >= 0; i--) {
                this->gradient_stack[i].SoftReset();
            }
        }



    };



}

#endif /* GRADIENTSTRUCTURE_HPP */
