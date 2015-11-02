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
#include <fstream>
#include <cmath>

#include <unordered_set>
//#define USE_BOOST
#ifdef USE_BOOST
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#else
#include <set>
#include "../Utilities/flat_map.hpp"
#include "../Utilities/flat_set.hpp"
#endif

#define ATL_USE_SMID

//#ifdef ATL_USE_SMID
//#include "../Utilities/SMID.hpp"
//#endif

#include "Variable.hpp"

//#define HESSIAN_TRACE 

//#warning add jacobian matrix calculations

#ifdef USE_BOOST
#define IDSet boost::container::flat_set
#else
#define IDSet flat_set
//
//flat_set
#endif
#define Entry StackEntry<REAL_T>



namespace atl {

    template<typename REAL_T>
    struct StackEntry {
        typedef typename VariableInfo<REAL_T>::HessianInfo HessianInfo;

        std::stack<size_t> entry_indexes;
        VariableInfo<REAL_T>* w; //function or dependent variable.
        IDSet<atl::VariableInfo<REAL_T>* > ids;
        //        flat_map<VariableInfo<REAL_T>*, REAL_T > first_order;
        //        flat_map<VariableInfo<REAL_T>*, flat_map<VariableInfo<REAL_T>*, REAL_T> > second_order;
        std::vector<VariableInfo<REAL_T>* > live_ids; //live variables used in reverse accumulation
        std::vector<atl::VariableInfo<REAL_T>* >id_list;
        std::vector<REAL_T> first;
        std::vector<REAL_T> second;
        size_t local_size;
        //        bool is_aliased;

        StackEntry() : w(NULL) {
        }

        inline void PushVariable(VariableInfo<REAL_T>* v) {
            if (v != w) {
                if (!ids.contains(v)) {
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

        inline void Prepare(atl::VariableInfo<REAL_T>* w) {
            id_list.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator e;
            e = ids.end();
            for (it = ids.begin(); it != e; ++it) {
                id_list.push_back((*it));
            }

            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator jt;
            typename std::vector<atl::VariableInfo<REAL_T>* >::iterator ee;
            ee = live_ids.end();
            for (jt = live_ids.begin(); jt != ee; ++jt) {
                REAL_T hij = 0.0;
                typename HessianInfo::iterator vijt;
                vijt = w->hessian_row.find((*jt)->id);
                if (vijt != w->hessian_row.end()) {
                    if ((*vijt) != 0) {
                        id_list.push_back((*jt));
                    }
                }

            }

        }

        inline void Reset() {
            first.resize(0);
            second.resize(0);
            typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            for (it = ids.begin(); it != ids.end(); ++it) {
                (*it)->Reset();
            }
            w->Reset();
            w = NULL;
            local_size = 0;
            live_ids.resize(0);
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
        std::vector<StackEntry<REAL_T> > gradient_stack;
        //        Adjoint<REAL_T>* gradient_stack;
        std::atomic<size_t> stack_current;
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

        /**
         * Sets the size of the stack.
         * @param size
         */
        void SetSize(size_t size) {
            gradient_stack.resize(size);
        }

        virtual ~GradientStructure() {

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

        inline StackEntry<REAL_T>& NextEntry() {
#ifdef ATL_ENABLE_BOUNDS_CHECKING
            if ((this->stack_current + 1) >= this->max_stack_size) {
                std::cout << "Current derivative stack index exceeds stack limits.\n" << std::flush;
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
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename flat_map<VariableInfo<REAL_T>*, REAL_T>::iterator git;

                this->gradient_stack[stack_current - 1].w->dvalue = 1.0;
#pragma unroll
                for (int i = (stack_current - 1); i >= 0; i--) {
                    w = this->gradient_stack[i].w->dvalue;
                    this->gradient_stack[i].w->dvalue = 0;
                    if (w != static_cast<REAL_T> (0)) {
                        int j = 0;
                        for (it = this->gradient_stack[i].ids.begin(); it != this->gradient_stack[i].ids.end(); ++it) {
                            (*it)->dvalue += w * this->gradient_stack[i].first[j];
                            j++;
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

            if (recording) {


                if (this->derivative_trace_level == GRADIENT) {
                    this->Accumulate();
                } else {
                    REAL_T w;
                    typedef typename VariableInfo<REAL_T>::HessianInfo HessianInfo;


                    //initialize w
                    this->gradient_stack[stack_current - 1].w->dvalue = 1.0;


                    unsigned rows = 0; //the size of the local derivatives, anything higher was pushed from previous calculation

                    std::vector<REAL_T> vij; //holds current second order derivative for i wrt j

                    //Hessian iterators
                    typename HessianInfo::iterator vit;
                    typename HessianInfo::iterator vijt;
                    typename HessianInfo::iterator iend;
                    typename HessianInfo::iterator jend;
                    typename HessianInfo::iterator vjt;

                    REAL_T hii = 0.0;
                    REAL_T hij = 0.0;
                    REAL_T hjk = 0;
                    REAL_T dj = 0;
                    REAL_T dk = 0;
                    for (int i = (stack_current - 1); i >= 0; i--) {

#ifdef HESSIAN_TRACE
                        std::cout << "\n\n\n";
                        std::stringstream ss;
                        int pushed_count = 0;
#endif
                        atl::VariableInfo<REAL_T>* vi = gradient_stack[i].w; //variable info for i
                        w = gradient_stack[i].w->dvalue; //set w
                        gradient_stack[i].w->dvalue = 0; //cancel out derivative for i

                        rows = gradient_stack[i].first.size();

                        //get h[i][i]
                        hii = 0.0;

                        iend = vi->hessian_row.end();
                        vit = vi->hessian_row.find(vi->id);
                        if (vit != iend) {
                            hii = (*vit);
                            if (hii != 0) {
                                (*vit) = 0.0;
                            }
                        }

                        //prepare for the hessian calculation. 
                        //builds a list of variables to use, statement level variables come first,
                        //then any pushed variables are after.
                        gradient_stack[i].Prepare();

                        //resize second order derivative for i wrt j
                        vij.resize(gradient_stack[i].id_list.size());

#pragma unroll
                        for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {
                            atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j];

                            //compute gradient
                            if (j < rows && w != 0) {
                                vj->dvalue += w * gradient_stack[i].first[j];
                            }

                            //load second order partial derivative for i wrt j and k
                            hij = 0.0;

                            vijt = vi->hessian_row.find(vj->id);
                            if (vijt != iend) {
                                hij = (*vijt);
                                (*vijt) = 0;
                            }
                            vij[j] = (hij);

                        }

#pragma unroll
                        //start the Hessian calculations    
                        for (unsigned j = 0; j < gradient_stack[i].id_list.size(); j++) {

                            atl::VariableInfo<REAL_T>* vj = gradient_stack[i].id_list[j]; //vj
                            jend = vj->hessian_row.end();
                            REAL_T hij = vij[j]; //h[i][j]

                            dj = 0;
                            bool j_in_local = false;
                            if (j < rows) {
                                dj = gradient_stack[i].first[j];
                                j_in_local = true;
                            }


                            bool j_pushed = false;
#pragma unroll
                            //use symmetry
                            for (unsigned k = 0; k <= j; k++) {
                                REAL_T entry = 0.0; //the entry value for h[j][k]
                                atl::VariableInfo<REAL_T>* vk = gradient_stack[i].id_list[k]; //vk

                                dk = 0;
                                bool k_in_local = false;
                                if (k < rows) {
                                    dk = gradient_stack[i].first[k];
                                    k_in_local = true;
                                }

                                if (!j_in_local && !k_in_local) {
                                    if (dk == 0 && dj == 0) {
                                        continue; //pushed variable doesn't matter
                                    }
                                }

                              

                                entry += vij[k] * dj + (hij * dk)+ hii * dj*dk;

#ifdef HESSIAN_TRACE
                                REAL_T s = 0;
#endif
                                //                                if (w != 0) {
                                if (j_in_local && k_in_local) {
                                    entry += w * gradient_stack[i].second[j * rows + k];
                                }
                                //                                    if (unsigned(j) < unsigned(rows)) {
                                //                                        if (unsigned(k) < unsigned(rows)) {
                                //                                            entry += w * gradient_stack[i].second[j * rows + k];
                                //#ifdef HESSIAN_TRACE
                                //                                            s = this->gradient_stack[i].second[j * rows + k];
                                //#endif
                                //                                        }
                                //                                    }
                                //                                }
#ifdef HESSIAN_TRACE
                                std::cout << "h[" << vj->id << "][" << vk->id << "] +=" << "h[" << vi->id << "][" << vj->id << "]{" << hij << "} *" << dk << " + " <<
                                        "h[" << vi->id << "][" << vk->id << "]{" << vij[k] << "}*" << dj << " + " <<
                                        "h[" << vi->id << "][" << vi->id << "]{" << hii << "} *" << dj << "*" << dk << " + " << w << "*" << s;
#endif

                                if (std::fabs(entry) > 0.0) {//h[j][k] needs to be updated

                                    //                                    hjk = 0.0;

                                    //                                    vjt = vj->hessian_row.find(vk->id);
                                    //                                    if (vjt != jend) {
                                    //                                        hjk = (*vjt);
                                    //                                    }
#ifdef HESSIAN_TRACE
                                    std::cout << " = " << (entry + hjk) << "\n";
                                    if (j > this->gradient_stack[i].first.size()) {
                                        pushed_count++;
                                        std::cout << "Push mattered!!!!\n";
                                    }
#endif
                                    //set h[j][k]
                                    vj->hessian_row[vk->id] += entry; // + hjk;
                                    if (j > 0 && j != k) {
                                        //set h[k][j]
                                        vk->hessian_row[vj->id] += entry; /// + hjk;
                                    }

                                    if (i > 0) {
                                        if (!j_pushed) {
                                            //this variable may be needed in the future, so push it to the next entry
                                            gradient_stack[i - 1].PushVariable(vj);
#ifdef HESSIAN_TRACE
                                            ss << vj->id << " ";
#endif
                                            j_pushed = true;
                                        }
                                    }
                                }
#ifdef HESSIAN_TRACE
                                else {
                                    std::cout << " = 0\n";
                                }
#endif
                            }
                        }

#ifdef HESSIAN_TRACE
                        std::cout << ss.str() << "\n";
#endif
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
