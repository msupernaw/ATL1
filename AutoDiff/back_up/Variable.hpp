/*
 * File:   Variable.hpp
 * Author: matthewsupernaw
 *
 * Created on September 24, 2014, 1:09 PM
 */

#ifndef ET4AD_VARIABLE_HPP
#define	ET4AD_VARIABLE_HPP

#include <cmath>
#include <stack>
#include <vector>
#include <valarray>
#include <thread>
#include "AlignedAllocator.hpp"
#include "Expression.hpp"
#include "GradientStructure.hpp"
#include "Add.hpp"
#include "Subtract.hpp"
#include "Multiply.hpp"
#include "Divide.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846	/* pi */
#endif



#ifdef ATL_USE_THREAD_LOCAL_STORAGE

#if defined (__GNUC__)
#define THREAD_LOCAL_STORAGE __thread
#elif defined (_MSC_VER)
#define THREAD_LOCAL_STORAGE __declspec(thread)
#else // !__GNUC__ && !_MSC_VER
#error "Define a thread local storage qualifier for your compiler/platform!"
#endif

#endif


namespace atl {
    
    /**
     * Base class for parameter transformations. Used in optimization
     * problems involving bounded parameters.
     * @param val
     * @param min
     * @param max
     * @return
     */
    template<typename REAL_T>
    class ParameterTransformation {
    public:
        /**
         * Convert the external value to its internal representation.
         *
         * @param val
         * @param min
         * @param max
         * @return
         */
        virtual REAL_T External2Internal(REAL_T val, REAL_T min_, REAL_T max_) const = 0;
        /**
         * Convert a variables internal value to its external representation.
         * @param val
         * @param min
         * @param max
         * @return
         */
        virtual REAL_T Internal2External(REAL_T val, REAL_T min_, REAL_T max_) const = 0;
        /**
         * The derivative of Internal2External.
         * @param val
         * @param min
         * @param max
         * @return
         */
        virtual REAL_T DerivativeInternal2External(REAL_T val, REAL_T min_, REAL_T max_) const = 0;
        
    };
    
    /**
     * Sine transformation for a variable.
     */
    template<typename REAL_T>
    class TanhParameterTransformation : public ParameterTransformation<REAL_T> {
    public:
        
        virtual REAL_T External2Internal(REAL_T val, REAL_T min_, REAL_T max_) const {
            return min_ + .5 * (max_ - min_)*(1.0 + std::tanh(val));
        }
        
        virtual REAL_T Internal2External(REAL_T val, REAL_T min_, REAL_T max_) const {
            return std::atanh(2.0 * (val - min_) / (max_ - min_) - 1.0);
        }
        
        virtual REAL_T DerivativeInternal2External(REAL_T val, REAL_T min_, REAL_T max_)const {
            return 2.0 / ((max_ - min_) * std::pow((1.0 - ((2.0 * (val - min_)) / max_ - min_ - 1.0)), 2.0));
        }
    };
    
    /**
     * Sine transformation for a variable.
     */
    template<typename REAL_T>
    class SinParameterTransformation : public ParameterTransformation<REAL_T> {
    public:
        
        virtual REAL_T External2Internal(REAL_T val, REAL_T min_, REAL_T max_) const {
            //                        if (val < min_) {
            //                            return -M_PI / 2.0;
            //                        } else if (val > max_) {
            //                            return M_PI / 2.0;
            //                        } else {
            //                            return 2.0 * M_PI + std::asin(std::max(-1.0, std::min(1.0, (2.0 * (val - min_) / (max_ - min_) - 1.0))));
            //                        }
            return  std::asin((2.0*val)/(max_-min_)-min_/(max_-min_)-max_/(max_-min_));
            ////
            //            REAL_T piby2 = 2. * std::atan(1.);
            //            REAL_T distnn = 8. * sqrt(std::numeric_limits<REAL_T>::epsilon());
            //            REAL_T vlimhi = piby2 - distnn;
            //            REAL_T vlimlo = -piby2 + distnn;
            //
            //            REAL_T yy = 2. * (val - min_) / (max_ - min_) - 1.;
            //            REAL_T yy2 = yy*yy;
            //            if (yy2 > (1. - std::numeric_limits<REAL_T>::epsilon())) {
            //                if (yy < 0.) {
            //                    // Lower limit
            //                    //       std::cout<<"SinParameterTransformation warning: is at its Lower allowed limit. "<<Value<<std::endl;
            //                    return vlimlo;
            //                } else {
            //                    // Upper limit
            //                    //       std::cout<<"SinParameterTransformation warning: is at its Upper allowed limit."<<std::endl;
            //                    return vlimhi;
            //                }
            //
            //            } else {
            //                return std::asin(yy);
            //            }
            
        }
        
        virtual REAL_T Internal2External(REAL_T val, REAL_T min_, REAL_T max_) const {
            
            //            return ((std::sin(val) + 1.0) / 2.0)*(max - min) + min;
            //                        return std::max(min_, std::min(max_, ((std::sin(val) + static_cast<REAL_T> (1.0)) / static_cast<REAL_T> (2.0))*(max_ - min_) + min_));
            return min_ + 0.5 * (max_ - min_)*(std::sin(val) + 1.);
        }
        
        virtual REAL_T DerivativeInternal2External(REAL_T val, REAL_T min_, REAL_T max_)const {
            return 0.5 * ((max_ - min_) * std::cos(val));
        }
    };
    
    template<typename REAL_T>
    class LogitParameterTransformation : public ParameterTransformation<REAL_T> {
    public:
        
        virtual REAL_T External2Internal(REAL_T val, REAL_T min_, REAL_T max_)const {
            REAL_T p = (val - min_) / (max_ - min_);
            return std::log(p / (1.0 - p));
        }
        
        virtual REAL_T Internal2External(REAL_T val, REAL_T min_, REAL_T max_) const {
            REAL_T p = ::exp(val) / (1.0 + ::exp(val));
            return p * (max_ - min_) + min_;
        }
        
        virtual REAL_T DerivativeInternal2External(REAL_T val, REAL_T min_, REAL_T max_)const {
            return (::exp(val) * ::log(M_E)*(max_ - min_)) / (::exp(val) + 1.0)-
            (::exp(2.0 * val) * ::log(M_E)*(max_ - min_)) / ::pow((::exp(val) + 1), 2.0);
        }
    };
    
    template<typename REAL_T, //base type
    int group = 0 > //group identifier
    class Variable : public atl::ExpressionBase<REAL_T, Variable<REAL_T, group > > {
        static SinParameterTransformation<REAL_T> default_transformation;
        IDSet<atl::VariableInfo<REAL_T>* > ids;
        VariableInfo<REAL_T>* mapped_info;
        ParameterTransformation<REAL_T>* transformation;
    public:
        typedef REAL_T BASE_TYPE;
        mutable VariableInfo<REAL_T>* info;
        REAL_T min_boundary_m;
        REAL_T max_boundary_m;
        
        
        std::string name_m;
        bool bounded_m;
        
        
        
        
        
        static /* ATTRIBUTE_TLS */ GradientStructure<REAL_T> gradient_structure_g;
        
        static bool IsRecording() {
            return Variable<REAL_T, group>::gradient_structure_g.recording;
        }
        
        static void SetRecording(bool record) {
            Variable<REAL_T, group>::gradient_structure_g.recording = record;
        }
        
        Variable(REAL_T val = 0.0) :
        info(new VariableInfo<REAL_T>),
        bounded_m(false),
        min_boundary_m(std::numeric_limits<REAL_T>::min()),
        max_boundary_m(std::numeric_limits<REAL_T>::max()),
        transformation(&default_transformation) {
            info->vvalue = (val);
            mapped_info = (this->info);
        }
        
        Variable(const Variable& other)
        : info(other.info),
        min_boundary_m(other.min_boundary_m),
        max_boundary_m(other.max_boundary_m),
        bounded_m(other.bounded_m),
        transformation(&default_transformation) {
            info->count++;
            mapped_info = (other.mapped_info);
        }
        
        //        Variable(Variable&& other)
        //        : info(other.info),
        //        min_boundary_m(other.min_boundary_m),
        //        max_boundary_m(other.max_boundary_m),
        //        bounded_m(other.bounded_m),
        //        transformation(other.transformation) {
        //            //                        info->count++;
        //            mapped_info = (other.mapped_info);
        //            //            other.info->Release();
        //            other.info = new atl::VariableInfo<REAL_T>();
        //            other.min_boundary_m = std::numeric_limits<REAL_T>::min();
        //            other.max_boundary_m = std::numeric_limits<REAL_T>::max();
        //            other.bounded_m = false;
        //            other.transformation = &default_transformation;
        //            other.mapped_info = NULL;
        //        }
        
        template<typename A>
        Variable(const ExpressionBase<REAL_T, A>& exp) :
        info(new VariableInfo<REAL_T>),
        bounded_m(false),
        min_boundary_m(std::numeric_limits<REAL_T>::min()),
        max_boundary_m(std::numeric_limits<REAL_T>::max()),
        transformation(&default_transformation) {
            mapped_info = (this->info);
            if (Variable<REAL_T>::gradient_structure_g.recording) {
                
                Adjoint<REAL_T>& entry = Variable<REAL_T>::gradient_structure_g.gradient_stack[Variable<REAL_T>::gradient_structure_g.NextIndex()];
                //                new(&entry)Adjoint<REAL_T>();
                entry.w = info;
                
                ids.clear();
                exp.PushIds(ids);
                size_t isize = ids.size();
                entry.entries.resize(isize);
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator jt;
                
                int i = 0;
                entry.second_order_partials.resize(ids.size() * ids.size());
                for (it = ids.begin(); it != ids.end(); ++it) {
                    REAL_T fx = exp.EvaluateDerivative((*it)->id);
                    int j = 0;
                    entry.entries[i] = (std::move(AdjointDerivative<REAL_T>((*it), fx)));
                    
                    if (Variable<REAL_T>::gradient_structure_g.derivative_trace_level == GRADIENT_AND_HESSIAN) {
                        for (jt = ids.begin(); jt != ids.end(); ++jt) {
                            REAL_T fxx = exp.EvaluateDerivative((*it)->id, (*jt)->id);
                            entry.second_order_partials[i * isize + j] = fxx;
                            j++;
                        }
                        
                    }
                    i++;
                }
            }
            info->vvalue = exp.GetValue();
        }
        
        virtual ~Variable() {
            info->Release();
        }
        
        /**
         * Allows for computations and gradient calculations separate from the
         * global gradient structure. Helpful for multi-threaded functions. It
         * allows for adjoint derivatives to be computed and set in the global
         * gradient structure with respect to the desired variables.
         *
         * @param gs
         * @param var
         * @param exp
         */
        template<typename A>
        inline void Assign(atl::GradientStructure<REAL_T>& gs, const atl::ExpressionBase<REAL_T, A>& exp) {
            if (gs.recording) {
                
                Adjoint<REAL_T>& entry = gs.gradient_stack[gs.NextIndex()];
                // new(&entry)Adjoint<REAL_T>();
                entry.w = info;
                
                ids.clear();
                exp.PushIds(ids);
                int isize = ids.size();
                entry.entries.resize(isize);
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator jt;
                
                int i = 0;
                entry.second_order_partials.resize(ids.size() * ids.size());
                for (it = ids.begin(); it != ids.end(); ++it) {
                    REAL_T fx = exp.EvaluateDerivative((*it)->id);
                    int j = 0;
                    entry.entries[i] = (std::move(AdjointDerivative<REAL_T>((*it), fx)));
                    
                    if (gs.derivative_trace_level == GRADIENT_AND_HESSIAN) {
                        for (jt = ids.begin(); jt != ids.end(); ++jt) {
                            REAL_T fxx = exp.EvaluateDerivative((*it)->id, (*jt)->id);
                            entry.second_order_partials[i * isize + j] = fxx;
                            j++;
                        }
                        
                    }
                    i++;
                }
            }
            this->SetValue(exp.GetValue());
        }
        
        inline Variable<REAL_T>& operator=(const REAL_T & value) {
            
            this->SetValue(value);
            return *this;
        }
        
        inline Variable<REAL_T>& operator=(const Variable<REAL_T> & other) {
            if (Variable<REAL_T>::gradient_structure_g.recording) {
                Adjoint<REAL_T>& entry = Variable<REAL_T>::gradient_structure_g.gradient_stack[Variable<REAL_T>::gradient_structure_g.NextIndex()];
                entry.w = info;
                entry.entries.push_back(AdjointDerivative<REAL_T>((other.info), 1.0));
                entry.second_order_partials.resize(1, 0.0);
                
                
            }
            this->SetValue(other.GetValue());
            return *this;
        }
        
        //        Variable<REAL_T>& operator=(Variable&& other) {
        //
        //            info = (other.info);
        //            min_boundary_m = (other.min_boundary_m);
        //            max_boundary_m = (other.max_boundary_m);
        //            bounded_m = (other.bounded_m);
        //            transformation = (other.transformation);
        //            //                                    info->count++;
        //            mapped_info = (other.mapped_info);
        //            other.info = new atl::VariableInfo<REAL_T>();
        //            other.min_boundary_m = std::numeric_limits<REAL_T>::min();
        //            other.max_boundary_m = std::numeric_limits<REAL_T>::max();
        //            other.bounded_m = false;
        //            other.transformation = &default_transformation;
        //            other.mapped_info = NULL;
        //            return *this;
        //        }
        
        //        operator REAL_T() {
        //            return this->GetValue();
        //        }
        //
        //        operator REAL_T()const {
        //            return this->GetValue();
        //        }
        
        template<class A>
        inline Variable& operator=(const ExpressionBase<REAL_T, A>& exp) {
            if (Variable<REAL_T>::gradient_structure_g.recording) {
                
                Adjoint<REAL_T>& entry = Variable<REAL_T>::gradient_structure_g.gradient_stack[Variable<REAL_T>::gradient_structure_g.NextIndex()];
                //                new(&entry)Adjoint<REAL_T>();
                entry.w = info;
                
                ids.clear();
                exp.PushIds(ids);
                size_t isize = 0;
                isize = ids.size();
                entry.entries.resize(isize);
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator jt;
                
                int i = 0;
                entry.second_order_partials.resize(ids.size() * ids.size());
                for (it = ids.begin(); it != ids.end(); ++it) {
                    REAL_T fx = exp.EvaluateDerivative((*it)->id);
                    int j = 0;
                    entry.entries[i] = (std::move(AdjointDerivative<REAL_T>((*it), fx)));
                    
                    if (Variable<REAL_T>::gradient_structure_g.derivative_trace_level == GRADIENT_AND_HESSIAN) {
                        for (jt = ids.begin(); jt != ids.end(); ++jt) {
                            REAL_T fxx = exp.EvaluateDerivative((*it)->id, (*jt)->id);
                            entry.second_order_partials[i * isize + j] = fxx;
                            j++;
                        }
                        
                    }
                    i++;
                }
            }
            this->SetValue(exp.GetValue());
            return *this;
        }
        
        inline Variable operator-() {
            return -1.0 * (*this);
        }
        
        inline Variable& operator+=(const REAL_T & val) {
            *this = *this+val;
            return *this;
        }
        
        inline Variable& operator+=(const Variable & other) {
            *this = *this+other;
            return *this;
        }
        
        template<class A>
        inline Variable& operator+=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this+exp;
            return *this;
        }
        
        inline Variable& operator-=(const REAL_T & val) {
            *this = *this-val;
            return *this;
        }
        
        inline Variable& operator-=(const Variable & other) {
            *this = *this-other;
            return *this;
        }
        
        template<class A>
        inline Variable& operator-=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this-exp;
            return *this;
        }
        
        inline Variable& operator*=(const REAL_T & val) {
            *this = *this*val;
            return *this;
        }
        
        inline Variable& operator*=(const Variable & other) {
            *this = *this*other;
            return *this;
        }
        
        template<class A>
        inline Variable& operator*=(const ExpressionBase<REAL_T, A>& exp) {
            *this = (*this) * exp;
            return *this;
        }
        
        inline Variable& operator/=(const REAL_T & val) {
            *this = *this / val;
            return *this;
        }
        
        inline Variable& operator/=(const Variable & other) {
            *this = *this / other;
            return *this;
        }
        
        template<class A>
        inline Variable& operator/=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this / exp;
            return *this;
        }
        
        inline Variable& operator++() {
            *this = *this+static_cast<REAL_T> (1.0);
            return *this;
        }
        
        inline const Variable operator++(int i) {
            Variable temp = *this;
            *this = static_cast<REAL_T> (1.0)+ (*this);
            return temp;
        }
        
        inline Variable& operator--() {
            *this = *this-static_cast<REAL_T> (1.0);
            return *this;
        }
        
        inline const Variable operator--(int i) {
            Variable temp = *this;
            *this = (*this) - static_cast<REAL_T> (1.0);
            return temp;
        }
        
        inline REAL_T GetValue() const {
            return info->vvalue;
        }
        
        inline REAL_T GetInternalValue() const {
            if (this->IsBounded()) {
                return this->transformation->External2Internal(this->GetValue(), this->GetMinBoundary(), this->GetMaxBoundary());
            } else {
                return this->GetValue();
            }
        }
        
        REAL_T GetScaledGradient(REAL_T x) {
            if (bounded_m) {
                return this->transformation->DerivativeInternal2External(x, min_boundary_m, max_boundary_m);
                //                return std::cos(x) * (max_boundary_m - min_boundary_m) / 2.0;
            } else {
                return 1.0;
            }
        }
        
        /**
         * If this variable is bounded, the value of
         * v will be transformed from internal to external
         * representation and set to the value of this
         * variable.
         * @param v
         */
        inline void UpdateValue(REAL_T v) {
            if (this->IsBounded()) {
                this->SetValue(this->transformation->Internal2External(v, this->GetMinBoundary(), this->GetMaxBoundary()));
            } else {
                this->SetValue(v);
            }
        }
        
        /**
         * Sets the value of this variable. If the variable is bounded,
         * the value will be set between the min and max boundary. If the
         * value is less than the minimum boundary, the variables value is
         * set to minimum boundary. If the  value is greater than the maximum
         * boundary, the variables value is set to maximum boundary. If the
         * value is signaling nan, the value is set to the mid point between
         * the min and max boundary.
         *
         * @param value
         */
        inline void SetValue(const REAL_T & value) {
            
            if (bounded_m) {
                
                if (value != value) {//nan
                    info->vvalue = min_boundary_m + (max_boundary_m - min_boundary_m) / static_cast<REAL_T> (2.0);
                    
                    return;
                }
                
                if (value < min_boundary_m) {
                    
                    info->vvalue = min_boundary_m;
                } else if (value > max_boundary_m) {
                    info->vvalue = max_boundary_m;
                    
                    
                } else {
                    info->vvalue = value;
                }
            } else {
                info->vvalue = value;
            }
        }
        
        ParameterTransformation<REAL_T>& GetParameterTransformation() {
            return this->transformation;
        }
        
        REAL_T GetMaxBoundary() const {
            return max_boundary_m;
        }
        
        void SetMaxBoundary(REAL_T max_boundary) {
            this->max_boundary_m = max_boundary;
        }
        
        REAL_T GetMinBoundary() const {
            return min_boundary_m;
        }
        
        void SetMinBoundary(REAL_T min_boundary) {
            this->min_boundary_m = min_boundary;
        }
        
        inline void SetBounds(REAL_T min_boundary, REAL_T max_boundary) {
            this->bounded_m = true;
            this->SetMinBoundary(min_boundary);
            this->SetMaxBoundary(max_boundary);
            //            this->SetValue((min_boundary+max_boundary)/2.0);
        }
        
        bool IsBounded() const {
            return bounded_m;
        }
        
        std::string GetName() const {
            return name_m;
        }
        
        void SetName(std::string name) {
            this->name_m = name;
        }
        
        inline void VariableCount(uint32_t & count) const {
            count++;
        }
        
        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids,bool include_dependent = true)const {
            ids.insert(this->info);
        }
        
        inline REAL_T EvaluateDerivative(uint32_t a) const {
            
            if (this->info->id == a) {
                return 1.0;
            } else {
                return 0.0;
            }
        }
        
        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            return 0.0;
        }
        
        static void ComputeGradient(GradientStructure<REAL_T>& gs, std::vector<atl::Variable<REAL_T>* >& variables, std::vector<REAL_T>& gradient) {
            gs.Accumulate();
            int size = variables.size();
            gradient.resize(size);
            for (int i = 0; i < size; i++) {
                gradient[i] = variables[i]->info->dvalue;
            }
        }
        
        static void ComputeGradient(GradientStructure<REAL_T>& gs, std::vector<atl::Variable<REAL_T>* >& variables, std::valarray<REAL_T>& gradient) {
            gs.Accumulate();
            int size = variables.size();
            gradient.resize(size);
            for (int i = 0; i < size; i++) {
                gradient[i] = variables[i]->info->dvalue;
            }
        }
        
        static void ComputeGradientAndHessian(GradientStructure<REAL_T>& gs, std::vector<atl::Variable<REAL_T>* >& variables, std::vector<REAL_T>& gradient, std::vector<std::vector<REAL_T> >& hessian) {
            gs.HessianAndGradientAccumulate();
            int size = variables.size();
            gradient.resize(size);
            hessian.resize(size);
            for (int i = 0; i < size; i++) {
                gradient[i] = variables[i]->info->dvalue;
                hessian[i].resize(size);
                for (int j = 0; j < size; j++) {
                    hessian[i][j] = variables[i]->info->GetHessianRowValue(variables[j]->info); //hessian_row[variables[j]->info];
                }
            }
        }
        
        /**
         * Builds an adjoint entry;
         * @param entry
         * @param dependent_variable
         * @param independent_variables
         * @param gradient
         * @param hessian
         */
        static void BuildAdjointEntry(atl::Adjoint<REAL_T>& entry,
                                      atl::Variable<REAL_T>& dependent_variable,
                                      std::vector<atl::Variable<REAL_T>* > independent_variables,
                                      std::vector<REAL_T>& gradient,
                                      std::vector< std::vector<REAL_T> >& hessian) {
            
            entry.w = dependent_variable.info;
            int n = independent_variables.size();
            
            if (n == gradient.size() && hessian.size() == n && hessian[0].size() == n) {
                entry.entries.resize(n);
                entry.second_order_partials.resize(n * n);
                for (int i = 0; i < n; i++) {
                    entry.entries[i] = (std::move(AdjointDerivative<REAL_T>(independent_variables[i]->info, gradient[i])));
                    for (int j = 0; j < n; j++) {
                        entry.second_order_partials[i * n + j] = hessian[i][j];
                    }
                }
            } else {
                std::cout << "Cannot build adjoint entry, dimensions do not match.\n";
                exit(0);
            }
            
        }
        
        VariableInfo<REAL_T>* GetMappedInfo() const {
            return mapped_info;
        }
        
        void SetMappedInfo(VariableInfo<REAL_T>* mapped_info) {
            this->mapped_info = mapped_info;
        }
        
        
    };
    
    template<typename REAL_T, int group>
    /* ATTRIBUTE_TLS */ GradientStructure<REAL_T> Variable<REAL_T, group>::gradient_structure_g;
    
    template<typename REAL_T, int group>
    SinParameterTransformation<REAL_T> Variable<REAL_T, group>::default_transformation;
    
    
    /**
     * Experimental code for forward mode AD. Uses the Transitionally Accumulated
     * Gradients method (TAG) from ET4AD.
     */
    namespace atl_forward_test {
        
        enum AD_Method {
            TAG,
            REVERSE
        };
        
        template<typename REAL_T, AD_Method method>
        struct AD_Variable : public atl::ExpressionBase<REAL_T, AD_Variable<REAL_T, method> > {
            
            AD_Variable() {
                
            }
            
            AD_Variable(const AD_Variable<REAL_T, method>& other) {
                
            }
            
            template<class A>
            inline AD_Variable(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator=(const AD_Variable<REAL_T, method>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, method>& operator=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator+=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator+=(const AD_Variable<REAL_T, method>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, method>& operator+=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator-=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator-=(const AD_Variable<REAL_T, method>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, method>& operator-=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator*=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator*=(const AD_Variable<REAL_T, method>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, method>& operator*=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator/=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, method>& operator/=(const AD_Variable<REAL_T, method>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, method>& operator/=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            
            
            
        };
        
        template<typename REAL_T>
        struct AD_Variable<REAL_T, TAG> {
            int id;
            
            flat_set<uint32_t> ids;
            flat_map<uint32_t, REAL_T> gradient;
            flat_map<uint32_t, flat_map<uint32_t, REAL_T> > hessian;
            typedef typename flat_set<uint32_t>::iterator iditerator;
            typedef typename flat_map<uint32_t, REAL_T>::iterator giterator;
            typedef typename flat_map<uint32_t, flat_map<uint32_t, REAL_T> >::iterator hiterator;
            
            AD_Variable() : id(atl::VariableIdGenerator::instance()->next()) {
                
            }
            
            AD_Variable(const AD_Variable<REAL_T, TAG>& other) : gradient(other.gradient), hessian(other.hessian) {
                
            }
            
            template<class A>
            inline AD_Variable(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator=(const AD_Variable<REAL_T, TAG>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, TAG>& operator=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator+=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator+=(const AD_Variable<REAL_T, TAG>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, TAG>& operator+=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator-=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator-=(const AD_Variable<REAL_T, TAG>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, TAG>& operator-=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator*=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator*=(const AD_Variable<REAL_T, TAG>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, TAG>& operator*=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator/=(const REAL_T& val) {
                
            }
            
            inline AD_Variable<REAL_T, TAG>& operator/=(const AD_Variable<REAL_T, TAG>& other) {
                
            }
            
            template<class A>
            inline AD_Variable<REAL_T, TAG>& operator/=(const ExpressionBase<REAL_T, A>& exp) {
                
            }
            
            inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids)const {
                ids.insert(this->info);
            }
            
            inline REAL_T EvaluateDerivative(uint32_t a) const {
                if (this->id == a) {
                    return 1.0;
                } else {
                    giterator it = this->gradient.find(a);
                    if (it != gradient.end()) {
                        return (*it).second;
                    } else {
                        return 0.0;
                    }
                }
            }
            
            inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
                hiterator it = this->hessian.find(a);
                if (it != hessian.end()) {
                    giterator git = (*it).find(b);
                    if (git != (*it).end()) {
                        return (*git).second;
                    } else {
                        return 0.0;
                    }
                }
            }
            
            
        };
        //
        //    template<typename REAL_T>
        //    struct AD_Variable<REAL_T, REVERSE> {
        //    };
        //
        
    }
}

template<typename REAL_T, int group>
std::istream& operator>>(std::istream& in, atl::Variable<REAL_T, group>& v) {
    REAL_T r;
    in>>r;
    v = r;
    return in;
}




#endif	/* VARIABLE_HPP */
