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
            return std::asin((2.0 * val) / (max_ - min_) - min_ / (max_ - min_) - max_ / (max_ - min_));
        }

        virtual REAL_T Internal2External(REAL_T val, REAL_T min_, REAL_T max_) const {
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
        mutable IDSet<atl::VariableInfo<REAL_T>* > ids_m;
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

        template<typename A>
        Variable(const ExpressionBase<REAL_T, A>& exp) :
        info(new VariableInfo<REAL_T>),
        bounded_m(false),
        min_boundary_m(std::numeric_limits<REAL_T>::min()),
        max_boundary_m(std::numeric_limits<REAL_T>::max()),
        transformation(&default_transformation) {
            mapped_info = (this->info);

            if (Variable<REAL_T>::gradient_structure_g.recording) {

                size_t index = Variable<REAL_T>::gradient_structure_g.NextIndex();
                StackEntry<REAL_T>& entry = Variable<REAL_T>::gradient_structure_g.gradient_stack[index];


                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator jt;

                entry.w = this->info;

                exp.PushIds(entry.ids, false);

                entry.first.resize(entry.ids.size());
                entry.second.resize(entry.ids.size() * entry.ids.size());
                //
                entry.local_size = entry.ids.size();
                int i = 0;
                REAL_T dx, dxx;
                for (it = entry.ids.begin(); it != entry.ids.end(); ++it) {

                    dx = exp.EvaluateDerivative((*it)->id);
                    entry.first[i] = dx;

                    if (Variable<REAL_T>::gradient_structure_g.derivative_trace_level == GRADIENT_AND_HESSIAN) {
                        int j = 0;
                        for (jt = entry.ids.begin(); jt <= it; ++jt) {
                            dxx = exp.EvaluateDerivative((*it)->id, (*jt)->id);
                            entry.second[i * entry.first.size() + j] = dxx;
                            if (i > 0 && i != j) {
                                entry.second[j * entry.first.size() + i] = dxx;
                            }
                            j++;
                        }
                    }
                    i++;
                }

            }
             this->SetValue(exp.GetValue());
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

                size_t index = gs.NextIndex();
                StackEntry<REAL_T>& entry = gs.gradient_stack[index];


                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator jt;

                entry.w = new VariableInfo<REAL_T>();

                exp.PushIds(entry.ids, false);

                entry.first.resize(entry.ids.size());
                entry.second.resize(entry.ids.size() * entry.ids.size());
                //
                entry.local_size = entry.ids.size();
                int i = 0;
                REAL_T dx, dxx;
                for (it = entry.ids.begin(); it != entry.ids.end(); ++it) {

                    dx = exp.EvaluateDerivative((*it)->id);
                    entry.first[i] = dx;

                    if (gs.derivative_trace_level == GRADIENT_AND_HESSIAN) {
                        int j = 0;
                        for (jt = entry.ids.begin(); jt <= it; ++jt) {
                            dxx = exp.EvaluateDerivative((*it)->id, (*jt)->id);
                            entry.second[i * entry.first.size() + j] = dxx;
                            if (i > 0 && i != j) {
                                entry.second[j * entry.first.size() + i] = dxx;
                            }
                            j++;
                        }
                    }
                    i++;
                }

                this->info->Release();
                this->info = entry.w;
            }

            this->SetValue(exp.GetValue());
      
        }

        inline Variable<REAL_T>& operator=(const REAL_T & value) {
            ids_m.clear();
            this->SetValue(value);
            return *this;
        }

        inline Variable<REAL_T>& operator=(const Variable<REAL_T> & other) {
                if (Variable<REAL_T>::gradient_structure_g.recording) {
                    StackEntry<REAL_T>& entry = Variable<REAL_T>::gradient_structure_g.gradient_stack[Variable<REAL_T>::gradient_structure_g.NextIndex()];
                    entry.ids.insert(other.info);
                    entry.w = info;
                    entry.first.push_back(1.0);
                    entry.second.resize(1, 0.0);
                    
                    
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

                size_t index = Variable<REAL_T>::gradient_structure_g.NextIndex();
                StackEntry<REAL_T>& entry = Variable<REAL_T>::gradient_structure_g.gradient_stack[index];


                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator jt;

                entry.w = new VariableInfo<REAL_T>();

                exp.PushIds(entry.ids, false);

                entry.first.resize(entry.ids.size());
                entry.second.resize(entry.ids.size() * entry.ids.size());
                //
                entry.local_size = entry.ids.size();
                int i = 0;
                REAL_T dx, dxx;
                for (it = entry.ids.begin(); it != entry.ids.end(); ++it) {

                    dx = exp.EvaluateDerivative((*it)->id);
                    entry.first[i] = dx;

                    if (Variable<REAL_T>::gradient_structure_g.derivative_trace_level == GRADIENT_AND_HESSIAN) {
                        int j = 0;
                        for (jt = entry.ids.begin(); jt <= it; ++jt) {
                            dxx = exp.EvaluateDerivative((*it)->id, (*jt)->id);
                            entry.second[i * entry.first.size() + j] = dxx;
                            if (i > 0 && i != j) {
                                entry.second[j * entry.first.size() + i] = dxx;
                            }
                            j++;
                        }
                    }
                    i++;
                }

                this->info->Release();
                this->info = entry.w;
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
            //            atl::Variable<REAL_T> temp;
            //            temp = *this+0.0;;
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

        inline void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids, bool include_dependent = true)const {
            ids.insert(this->info);
            //            if (include_dependent) {
            //                if (this->info->ids.size() > 0) {
            //                    typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            //                    for (it = this->info->ids.begin(); it != this->info->ids.end(); ++it) {
            //                        ids.insert((*it));
            //                    }
            //                }
            //            }
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            ids.insert(this->info->id);
            //            if (this->ids_m.size() > 0) {
            //                typename IDSet<atl::VariableInfo<REAL_T>* >::iterator it;
            //                for (it = this->info->ids.begin(); it != this->info->ids.end(); ++it) {
            //                    ids.insert((*it)->id);
            //                }
            //            }

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
         * These structures must be sorted by dependent_variable.
         * @param entry
         * @param dependent_variable
         * @param independent_variables
         * @param gradient
         * @param hessian
         */
        static void BuildAdjointEntry(atl::StackEntry<REAL_T>& entry,
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
                    entry.ids.insert(dependent_variable[i]);
                    entry.first[i] = gradient[i];
                    for (int j = 0; j < n; j++) {
                        entry.second[i * n + j] = hessian[i][j];
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
    //    namespace atl_forward_test {

    //        enum AD_Method {
    //            TAG,
    //            REVERSE
    //        };
    //
    //        template<typename REAL_T, AD_Method method>
    //        struct AD_Variable : public atl::ExpressionBase<REAL_T, AD_Variable<REAL_T, method> > {
    //            //            flat_map<uint32_t, double> gradient;
    //            //            flat_map<uint32_t, flat_map<uint32_t, double> > hessian;
    //            //            bool independent_variable;
    //
    //            AD_Variable() {
    //
    //            }
    //
    //            AD_Variable(const AD_Variable<REAL_T, method>& other) {
    //
    //            }
    //
    //            template<class A>
    //            inline AD_Variable(const ExpressionBase<REAL_T, A>& exp) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator=(const REAL_T& val) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator=(const AD_Variable<REAL_T, method>& other) {
    //
    //            }
    //
    //            template<class A>
    //            inline AD_Variable<REAL_T, method>& operator=(const ExpressionBase<REAL_T, A>& exp) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator+=(const REAL_T& val) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator+=(const AD_Variable<REAL_T, method>& other) {
    //
    //            }
    //
    //            template<class A>
    //            inline AD_Variable<REAL_T, method>& operator+=(const ExpressionBase<REAL_T, A>& exp) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator-=(const REAL_T& val) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator-=(const AD_Variable<REAL_T, method>& other) {
    //
    //            }
    //
    //            template<class A>
    //            inline AD_Variable<REAL_T, method>& operator-=(const ExpressionBase<REAL_T, A>& exp) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator*=(const REAL_T& val) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator*=(const AD_Variable<REAL_T, method>& other) {
    //
    //            }
    //
    //            template<class A>
    //            inline AD_Variable<REAL_T, method>& operator*=(const ExpressionBase<REAL_T, A>& exp) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator/=(const REAL_T& val) {
    //
    //            }
    //
    //            inline AD_Variable<REAL_T, method>& operator/=(const AD_Variable<REAL_T, method>& other) {
    //
    //            }
    //
    //            template<class A>
    //            inline AD_Variable<REAL_T, method>& operator/=(const ExpressionBase<REAL_T, A>& exp) {
    //
    //            }
    //
    //
    //
    //
    //        };

    template<typename REAL_T>
    class TAD_Variable : public atl::ExpressionBase<REAL_T, TAD_Variable<REAL_T> > {
    public:
        int id;
        static bool compute_hessian;
        bool is_independent;
        size_t max_dependent;
        std::vector<REAL_T> gradient;
        std::vector<REAL_T> hessian;

        flat_set<uint32_t> dependent_ids;
        //        flat_map<uint32_t, REAL_T> gradient;
        //        flat_map<uint32_t, flat_map<uint32_t, REAL_T> > hessian;
        typedef typename flat_set<uint32_t>::iterator iditerator;
        typedef typename flat_set<uint32_t>::const_iterator ciditerator;
        //        typedef typename flat_map<uint32_t, REAL_T>::iterator giterator;
        //        typedef typename flat_map<uint32_t, flat_map<uint32_t, REAL_T> >::iterator hiterator;
        //        typedef typename flat_map<uint32_t, REAL_T>::const_iterator cgiterator;
        //        typedef typename flat_map<uint32_t, flat_map<uint32_t, REAL_T> >::const_iterator chiterator;
        REAL_T value;

        TAD_Variable() : id(atl::VariableIdGenerator::instance()->next()) {
            max_dependent = 0;
        }

        TAD_Variable(const TAD_Variable<REAL_T>& other) : gradient(other.gradient), hessian(other.hessian) {

        }

        inline TAD_Variable(const REAL_T& val) : value(val), id(atl::VariableIdGenerator::instance()->next()) {
        }

        template<class A>
        inline TAD_Variable(const ExpressionBase<REAL_T, A>& exp) {
            exp.PushIds(dependent_ids);

            if (dependent_ids.size() > 0) {
                size_t size = *(dependent_ids.end() - 1);
                max_dependent = size + 1;
                this->gradient.resize(size + 1);
                if (TAD_Variable<REAL_T>::compute_hessian) {
                    this->hessian.resize((size + 1)*(size + 1));
                }
                iditerator it;
                iditerator jt;
                REAL_T dx = 0.0;
                REAL_T dx2 = 0.0;
                for (it = dependent_ids.begin(); it != dependent_ids.end(); ++it) {
                    if (TAD_Variable<REAL_T>::compute_hessian) {
                        for (jt = dependent_ids.begin(); jt != it; ++jt) {
                            dx2 = exp.EvaluateDerivative((*it), (*jt));
                            this->hessian[(*it)*(size + 1)+ (*jt)] = dx2;
                            this->hessian[(*jt)*(size + 1)+ (*it)] = dx2;
                        }
                    }
                    dx = exp.EvaluateDerivative((*it));
                    this->gradient[(*it)] = dx;
                }
            }
            this->value = exp.GetValue();
        }

        inline const REAL_T GetValue() const {
            return this->value;
        }

        inline TAD_Variable& operator=(const REAL_T& val) {
            this->value = val;
            return *this;
        }

        inline TAD_Variable& operator=(const TAD_Variable<REAL_T>& other) {
            other.PushIds(dependent_ids);
            if (dependent_ids.size() > 0) {
                size_t size = *(dependent_ids.end() - 1);
                max_dependent = size + 1;
                this->gradient.resize(size + 1);
                if (TAD_Variable<REAL_T>::compute_hessian) {
                    this->hessian.resize((size + 1)*(size + 1));
                }
                iditerator it;
                iditerator jt;
                REAL_T dx = 0.0;
                REAL_T dx2 = 0.0;
                for (it = dependent_ids.begin(); it != dependent_ids.end(); ++it) {
                    if (TAD_Variable<REAL_T>::compute_hessian) {
                        for (jt = dependent_ids.begin(); jt != it; ++jt) {
                            dx2 = other.EvaluateDerivative((*it), (*jt));
                            this->hessian[(*it)*(size + 1)+ (*jt)] = dx2;
                            this->hessian[(*jt)*(size + 1)+ (*it)] = dx2;
                        }
                    }
                    dx = other.EvaluateDerivative((*it));
                    this->gradient[(*it)] = dx;
                }
            }
            this->value = other.GetValue();
            return *this;
        }

        template<class A>
        inline TAD_Variable& operator=(const ExpressionBase<REAL_T, A>& exp) {
            exp.PushIds(dependent_ids);

            if (dependent_ids.size() > 0) {
                size_t size = *(dependent_ids.end() - 1);
                max_dependent = size + 1;
                this->gradient.resize(size + 1);
                if (TAD_Variable<REAL_T>::compute_hessian) {
                    this->hessian.resize((size + 1)*(size + 1));
                }
                iditerator it;
                iditerator jt;
                REAL_T dx = 0.0;
                REAL_T dx2 = 0.0;
                for (it = dependent_ids.begin(); it != dependent_ids.end(); ++it) {
                    if (TAD_Variable<REAL_T>::compute_hessian) {
                        for (jt = dependent_ids.begin(); jt != it; ++jt) {
                            dx2 = exp.EvaluateDerivative((*it), (*jt));
                            this->hessian[(*it)*(size + 1)+ (*jt)] = dx2;
                            this->hessian[(*jt)*(size + 1)+ (*it)] = dx2;
                        }
                    }
                    dx = exp.EvaluateDerivative((*it));
                    this->gradient[(*it)] = dx;
                }
            }
            this->value = exp.GetValue();

            return *this;
        }

        inline TAD_Variable& operator+=(const REAL_T& val) {
            *this = *this+val;
            return *this;
        }

        inline TAD_Variable& operator+=(const TAD_Variable<REAL_T>& other) {
            *this = *this+other;
            return *this;
        }

        template<class A>
        inline TAD_Variable& operator+=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this+exp;
            return *this;
        }

        inline TAD_Variable& operator-=(const REAL_T& val) {
            *this = *this-val;
            return *this;

        }

        inline TAD_Variable& operator-=(const TAD_Variable<REAL_T>& other) {
            *this = *this-other;
            return *this;
        }

        template<class A>
        inline TAD_Variable& operator-=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this-exp;
            return *this;
        }

        inline TAD_Variable& operator*=(const REAL_T& val) {
            *this = *this*val;
            return *this;
        }

        inline TAD_Variable& operator*=(const TAD_Variable<REAL_T>& other) {
            *this = *this*other;
            return *this;
        }

        template<class A>
        inline TAD_Variable& operator*=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this*exp;
            return *this;
        }

        inline TAD_Variable& operator/=(const REAL_T& val) {
            *this = *this / val;
            return *this;
        }

        inline TAD_Variable& operator/=(const TAD_Variable<REAL_T>& other) {
            *this = *this / other;
            return *this;
        }

        template<class A>
        inline TAD_Variable& operator/=(const ExpressionBase<REAL_T, A>& exp) {
            *this = *this / exp;
            return *this;
        }

        inline void PushIds(IDSet<uint32_t >& ids)const {
            if (this->id != 0) {
                ids.insert(this->id);
            } else
                if (this->dependent_ids.size() > 0) {
                ciditerator it;
                for (it = this->dependent_ids.begin(); it != this->dependent_ids.end(); ++it) {
                    ids.insert((*it));
                }
            }
        }

        REAL_T WRT(const TAD_Variable<REAL_T>& v) {
            iditerator it = this->dependent_ids.find(v.id);
            if (it != this->dependent_ids.end()) {
                return gradient[v.id];
            }
            return 0.0;
        }

        REAL_T WRT(const TAD_Variable<REAL_T>& x, const TAD_Variable<REAL_T>& y) {
            uint32_t max_ = std::max(x.id, y.id);
            uint32_t min_ = std::min(x.id, y.id);

            iditerator it = this->dependent_ids.find(max_);
            if (it != this->dependent_ids.end()) {
                return hessian[max_ * gradient.size() + min_];
            }


            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t a) const {
            //            return this->WRT(
            //            if (this->id != 0) {
            //                if (this->id == a) {
            //                    return 1.0;
            //                } else {
            //                    return 0;
            //                }
            //            } else {
            //                if (gradient.size() > 0) {
            //                    cgiterator it = this->gradient.find(a);
            //                    if (it != gradient.end()) {
            //                        //                    std::cout << "found id = a-> " << (*it).second << "\n";
            //                        return (*it).second;
            //                    } else {
            //                        return 0.0;
            //                    }
            //                }
            //            }
            //         
            if (a < gradient.size()) {
                return gradient[a];
            }
            return 0.0;
        }

        inline REAL_T EvaluateDerivative(uint32_t a, uint32_t b) const {
            //
            //
            //            if (this->id == a && this->id == b) {
            //                return 0;
            //            }
            //            uint32_t max_ = std::max(a,b);
            //            uint32_t min_ = std::min(a,b);

            if (a < gradient.size() && b < gradient.size()) {
                //                std::cout<<a<<","<<b<<" "<<max_dependent<<"\n";
                return hessian[a * gradient.size() + b];
            }


            return 0.0;

            //            if (hessian.size() > 0) {
            //                chiterator it = this->hessian.find(a);
            //                if (it != hessian.end()) {
            //                    if ((*it).second.size() > 0) {
            //                        cgiterator git = (*it).second.find(b);
            //                        if (git != (*it).second.end()) {
            //                            return (*git).second;
            //                        } else {
            //                            return 0.0;
            //                        }
            //                    } else {
            //                        return 0.0;
            //                    }
            //                }
            //            }
            return 0.0;
        }


    };
    template<class REAL_T>
    bool TAD_Variable<REAL_T>::compute_hessian = false;
    //
    //    template<typename REAL_T>
    //    struct AD_Variable<REAL_T, REVERSE> {
    //    };
    //

    //    }
}

template<typename REAL_T, int group>
std::istream& operator>>(std::istream& in, atl::Variable<REAL_T, group>& v) {
    REAL_T r;
    in>>r;
    v = r;
    return in;
}




#endif	/* VARIABLE_HPP */
