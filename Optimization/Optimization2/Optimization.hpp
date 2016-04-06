/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FunctionMinimizer.hpp
 * Author: matthewsupernaw
 *
 * Created on February 9, 2016, 1:29 PM
 */

#ifndef FUNCTIONMINIMIZER_HPP
#define FUNCTIONMINIMIZER_HPP

#include "../../AutoDiff/AutoDiff.hpp"
#include "../../Utilities/StringUtil.hpp"
//#include "../ATL/Containers/Containers.hpp"
#include <vector>
#include "DerivativeChecker.hpp"
#include <iomanip>

namespace atl {

    template<class T>
    class OptimizationRoutine;

    template<class T>
    class MCMC;

    template<class T>
    class MinimizerResults {
        T function_value;
        bool converged;
        T convergence_criteria;
        int max_phases;

        std::vector<atl::Variable<T>* > hyper_parameters_m;

        std::vector<T> gradient;
        std::vector<std::vector<T> > hyper_parameter_hessian;
        T log_determinant_of_hyper_parmameter_hessian;
        std::vector<std::vector<T> > hyper_parameter_var_covar_matrix;
        std::vector<T> hyper_parameter_standard_deviation;

        std::vector<atl::Variable<T>* > random_variables_m;
        std::vector<std::vector<T> > random_variable_hessian;
        T log_determinant_of_random_variable_hessian;
        std::vector<std::vector<T> > random_variable_var_covar_matrix;
        std::vector<T> random_variable_standard_deviation;
    };

    template<class T>
    class ObjectiveFunction : public DerivativeChecker<T> {
        std::vector<atl::Variable<T>* > hyper_parameters_m;
        std::vector<int> hyper_parameter_phases_m;
        std::vector<atl::Variable<T>* > random_variables_m;
        std::vector<int> random_variable_phases_m;
        int max_phase_m = 0;
        friend class OptimizationRoutine<T>;
        friend class MCMC<T>;
    public:

        inline operator atl::Variable<T>() {
            return this->Evaluate();
        }

        virtual void Initialize() {
        }

        virtual const atl::Variable<T> Evaluate() = 0;

        void RegisterHyperParameter(atl::Variable<T>& v, int phase = 1) {
            this->hyper_parameters_m.push_back(&v);
            this->hyper_parameter_phases_m.push_back(phase);
            if (phase > max_phase_m) {
                max_phase_m = phase;
            }
        }

        void RegisterRandomVariable(atl::Variable<T>& v, int phase = 1) {
            this->random_variables_m.push_back(&v);
            this->random_variable_phases_m.push_back(phase);
            if (phase > max_phase_m) {
                max_phase_m = phase;
            }
        }

        virtual void Objective_Function(atl::Variable<T>& ff) {
            ff = this->Evaluate();
        }


    };

    template<class T>
    class OptimizationRoutine {
    protected:

        class Cholesky {
            typedef std::vector<std::vector<T> > CMatrix;
            typedef std::vector<T> CVector;
            typedef std::vector<std::vector<unsigned> > CPattern;



        public:
            size_t n;
            CMatrix el;
            CPattern pattern;
            CPattern upper_pattern;
            CPattern lower_pattern;

            bool pd;
            bool has_pattern;

            Cholesky() {

            }

            Cholesky(const CMatrix& a, bool store_pattern = true) {
                el = a;
                pd = true;
                n = a.size();
                has_pattern = store_pattern;

                if (has_pattern) {

                    this->pattern.resize(this->n * this->n);
                    this->lower_pattern.resize(n);
                    this->upper_pattern.resize(n);
                    for (int i = 0; i < n; i++) {
                        for (int j = i; j < n; j++) {
                            T sum = 0;
                            int k;
                            for (sum = el[i][j], k = i - 1; k >= 0; k--) {
                                sum -= el[i][k] * el[j][k];
                                if (el[i][k] != 0.0 && el[j][k] != 0.0) {
                                    pattern[i * n + j].push_back(k);

                                }
                            }
                            if (i == j) {
                                if (sum <= 0.0) {
                                    pd = false;
                                    std::cout << "matrix is not positive definite.\n";
                                }
                                el[i][i] = std::sqrt(sum);
                            } else {
                                el[j][i] = sum / el[i][i];
                            }
                        }
                    }
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < i; j++)
                            el[j][i] = 0.0;

                    for (int i = 0; i < n; i++) {
                        int k;
                        for (k = i - 1; k >= 0; k--) {
                            if (el[i][ k] != 0.0) {
                                this->lower_pattern[i].push_back(k);
                            }
                        }
                    }

                    for (int i = n - 1; i >= 0; i--) {
                        int k;
                        for (k = i + 1; k < n; k++) {
                            if (el[k][i] != 0.0) {
                                this->upper_pattern[i].push_back(k);
                            }
                        }
                    }

                } else {


                    for (int i = 0; i < n; i++) {
                        for (int j = i; j < n; j++) {
                            T sum = 0;
                            int k;
                            for (sum = el[i][j], k = i - 1; k >= 0; k--) {
                                sum -= el[i][k] * el[j][k];
                            }
                            if (i == j) {
                                if (sum <= 0.0) {
                                    pd = false;
                                    std::cout << "matrix is not positive definite.\n";
                                }
                                el[i][i] = std::sqrt(sum);
                            } else {
                                el[j][i] = sum / el[i][i];
                            }
                        }
                    }
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j < i; j++)
                            el[j][i] = 0.0;


                }




            }

            void recompute(const CMatrix& a) {
                el = a;
                pd = true;
                n = a.size();
                this->pattern.resize(0);
                this->pattern.resize(this->n * this->n);
                this->lower_pattern.resize(0);
                this->lower_pattern.resize(n);
                this->upper_pattern.resize(0);
                this->upper_pattern.resize(n);

                for (int i = 0; i < n; i++) {
                    for (int j = i; j < n; j++) {
                        T sum = 0;
                        int k;
                        for (sum = el[i][j], k = i - 1; k >= 0; k--) {
                            sum -= el[i][k] * el[j][k];
                            if (el[i][k] != 0.0 && el[j][k] != 0.0) {
                                pattern[i * n + j].push_back(k);
                            }
                        }
                        if (i == j) {
                            if (sum <= 0.0) {
                                pd = false;
                                std::cout << "matrix is not positive definite.\n";
                            }
                            el[i][i] = std::sqrt(sum);
                        } else {
                            el[j][i] = sum / el[i][i];
                        }
                    }
                }
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < i; j++)
                        el[j][i] = 0.0;

                for (int i = 0; i < n; i++) {
                    int k;
                    for (k = i - 1; k >= 0; k--) {
                        if (el[i][ k] != 0.0) {
                            this->lower_pattern[i].push_back(k);
                        }
                    }
                }

                for (int i = n - 1; i >= 0; i--) {
                    int k;
                    for (k = i + 1; k < n; k++) {
                        if (el[k][i] != 0.0) {
                            this->upper_pattern[i].push_back(k);
                        }
                    }
                }
            }

            void compute(const CMatrix& a) {
                el = a;
                pd = true;
                for (int i = 0; i < n; i++) {
                    for (int j = i; j < n; j++) {
                        T sum = 0;
                        int k;
                        //                sum = el[i][j];
#pragma unroll
                        for (sum = el[i][j], k = pattern[i * n + j].size() - 1; k >= 0; k--) {
                            sum -= el[i][pattern[i * n + j][k]] * el[j][pattern[i * n + j][k]];
                        }
                        if (i == j) {
                            if (sum <= 0.0) {
                                pd = false;
                                std::cout << "matrix is not positive definite.\n";
                            }
                            el[i][i] = std::sqrt(sum);
                        } else {
                            el[j][i] = sum / el[i][i];
                        }
                    }
                }
#pragma unroll
                for (int i = 0; i < n; i++)
#pragma unroll
                    for (int j = 0; j < i; j++)
                        el[j][i] = 0.0;

            }

            const CVector solve(const CMatrix& a, const CVector& b, bool recompute = false) {

                if (b.size() != this->n) {
                    std::cout << "error!!!!\n";
                }

                CVector x(this->n);
                if (!this->has_pattern) {
                    std::cout << "no pattern.\n";
                    return x;
                }

                recompute ? this->recompute(a) : this->compute(a);
                // solve L y = b, storing y in x
                for (int i = 0; i < n; i++) {
                    T sum;
                    int k;
                    for (sum = b[i], k = i - 1; k >= 0; k--)
                        sum -= el[i][k] * x[k];
                    x[i] = (sum / el[i][i]);
                }
                // solve L^T x = y
                for (int i = n - 1; i >= 0; i--) {
                    double sum;
                    int k;
                    for (sum = x[i], k = i + 1; k < n; k++)
                        sum -= el(k, i) * x[k];
                    x[i] = (sum / el[i][i]);
                }
                return x;

            }

            void solve(const CVector& b, CVector& x) {
                if (b.size() != this->n) {
                    std::cout << "error!!!!\n";
                }


                // solve L y = b, storing y in x
#pragma unroll
                for (int i = 0; i < n; i++) {
                    T sum;
                    int k;
#pragma unroll
                    for (sum = b[i], k = 0; k < this->lower_pattern[i].size(); k++)
                        sum -= el[i][this->lower_pattern[i][k]] * x[this->lower_pattern[i][k]];
                    x[i] = (sum / el[i][i]);
                }
                // solve L^T x = y
#pragma unroll
                for (int i = n - 1; i >= 0; i--) {
                    double sum;
                    int k;
#pragma unroll
                    for (sum = x[i], k = 0; k < this->upper_pattern[i].size(); k++)
                        sum -= el[this->upper_pattern[i][k]][i] * x[this->upper_pattern[i][k]];
                    x[i] = (sum / el[i][i]);
                }

            }

            void solve(CMatrix &b, CMatrix &x) {
                int i, j, m = b.size();
                if (b.size() != n || x.size() != n || b[0].size() != x[0].size())
                    throw ("Cholesky::solve bad sizes");
                CVector xx(n);
                for (j = 0; j < m; j++) {
                    for (i = 0; i < n; i++) xx[i] = b[i][j];
                    solve(xx, xx);
                    for (i = 0; i < n; i++) x[i][j] = xx[i];
                }
            }

            void inverse(CMatrix& ainv) {
                for (int i = 0; i < n; i++)
                    for (int j = 0; j <= i; j++) {
                        T sum = i == j ? 1.0 : 0.0;
                        for (int k = i - 1; k >= j; k--)
                            sum -= el[i][k] * ainv[j][k];
                        ainv[j][i] = sum / el[i][i];
                    }
                for (int i = n - 1; i >= 0; i--)
                    for (int j = 0; j <= i; j++) {
                        T sum = i < j ? 0.0 : ainv[j][i];
                        for (int k = i + 1; k < n; k++)
                            sum -= el[k][i] * ainv[j][k];
                        T temp = sum / el[i][i];
                        ainv[i][j] = temp;
                        ainv[j][i] = temp;
                    }
            }

            const T logdet() {
                T sum = 0.0;
                for (int i = 0; i < n; i++)
                    sum += std::log(std::fabs(el[i][i]));
                return static_cast<T> (2.0) * sum;
            }

            const T det() {
                return std::exp(logdet());
            }



        private:

            inline int index(int i, int j) {
                return i * n + j;
            }


        };

        ObjectiveFunction<T>* objective_function_m = NULL;
        std::vector<atl::Variable<T>* > hyper_parameters_m;
        std::vector<atl::Variable<T>* > random_variables_m;
        std::vector<atl::Variable<T>* > all_variables_m;

        //runtime 
        T function_value;
        T inner_function_value;
        T tolerance = 1e-4;
        T maxgc;
        T inner_maxgc;
        int phase_m;
        bool re_active = true;
        uint32_t max_line_searches = 1000;
        uint32_t max_iterations = 10000;
        size_t max_history = 1000;

        std::valarray<T> x;
        std::valarray<T> best;
        std::valarray<T> gradient;
        std::valarray<std::valarray<T> > hessian;

        std::valarray<T> inner_x;
        std::valarray<T> inner_best;
        std::valarray<T> inner_gradient;
        std::valarray<std::valarray<T> > inner_hessian;
        atl::GradientStructure<T> inner_gs;
        atl::Variable<T> log_det;

        //
        std::vector<std::vector<T> > rand_hessian; //(this->random_variables_m.size(), std::vector<T>(this->random_variables_m.size()));
        std::vector<std::vector<T> > rand_hessian_inv; //(this->random_variables_m.size(), std::vector<T>(this->random_variables_m.size()));
        std::vector<std::vector<T> > rand_hessian_dx; //(this->random_variables_m.size(), std::vector<T>(this->random_variables_m.size()));
        std::vector<std::vector<T> > ret;
//        atl::Matrix<T > atl_rand_hessian;
//        std::vector<T> sse_rand_hessian_inv;
//        std::vector<T> sse_rand_hessian_dx;
        std::vector<T> results;
        long outer_iteration;

        Cholesky cholesky_outer;
        Cholesky cholesky_inner;
    public:

        OptimizationRoutine(ObjectiveFunction<T>* objective_function = NULL) :
        objective_function_m(objective_function) {
        }

        ObjectiveFunction<T>* GetObjectiveFunction() const {
            return objective_function_m;
        }

        void SetObjectiveFunction(ObjectiveFunction<T>* objective_function) {
            this->objective_function_m = objective_function;
        }

        bool Run() {
            if (this->objective_function_m == NULL) {
                return false;
            }
            bool success = false;
            for (int i = 1; i <= this->objective_function_m->max_phase_m; i++) {
                this->Prepare(i);
                phase_m = i;
                success = this->Evaluate();
            }
            return success;
        }

        const MinimizerResults<T> GetResults() {

        }

    protected:

        void Prepare(int phase) {
            this->outer_iteration = 0;
            this->hyper_parameters_m.resize(0);
            this->random_variables_m.resize(0);
            this->all_variables_m.resize(0);
            for (int i = 0; i < this->objective_function_m->hyper_parameter_phases_m.size(); i++) {
                if (this->objective_function_m->hyper_parameter_phases_m[i] <= phase) {
                    this->hyper_parameters_m.push_back(this->objective_function_m->hyper_parameters_m[i]);
                    this->all_variables_m.push_back(this->objective_function_m->hyper_parameters_m[i]);
                }
            }

            for (int i = 0; i < this->objective_function_m->random_variable_phases_m.size(); i++) {
                if (this->objective_function_m->random_variable_phases_m[i] <= phase) {
                    this->random_variables_m.push_back(this->objective_function_m->random_variables_m[i]);
                    this->all_variables_m.push_back(this->objective_function_m->random_variables_m[i]);
                }
            }

            if (this->random_variables_m.size() > 0) {
                rand_hessian.resize(this->random_variables_m.size(), std::vector<T>(this->random_variables_m.size()));
                rand_hessian_dx.resize(this->random_variables_m.size(), std::vector<T>(this->random_variables_m.size()));
                ret.resize(this->random_variables_m.size(), std::vector<T>(this->random_variables_m.size()));
            }

        }

        void CallInnerObjectiveFunction(atl::Variable<T>& f) {
            atl::Variable<T>::gradient_structure_g.Reset();
            f = this->objective_function_m->Evaluate();
        }

        void CallObjectiveFunction(atl::Variable<T>& f) {
            atl::Variable<T>::gradient_structure_g.Reset();
            if (this->random_variables_m.size() > 0) {//this is a laplace problem
                this->EvaluateLaplace(f);
            } else {
                f = this->objective_function_m->Evaluate();
            }
        }

        void TransposeMatrix(std::vector<std::vector<T> >& m) {

#pragma unroll
            for (int i = 0; i < m.size(); i++) {
#pragma unroll
                for (int j = i + 1; j < m.size(); j++) {
                    std::swap(m[i][j], m[j][i]);
                }
            }
        }

        void TransposeMatrix(std::vector<T>& m) {
            size_t size = this->random_variables_m.size();
#pragma unroll
            for (int i = 0; i < size; i++) {
#pragma unroll
                for (int j = i + 1; j < size; j++) {
                    std::swap(m[i * size + j], m[j * size + i]);
                }
            }
        }

        void MatrixMultiply(std::vector<std::vector<T> >& a, std::vector<std::vector<T> >& b, std::vector<std::vector<T> >& ret) {
            for (int i = 0; i < a.size(); i++)
                for (int j = 0; j < a.size(); j++)
                    for (int k = 0; k < a.size(); k++) {
                        ret[i][j] += a[i][k] * b[k][j];
                    }
        }

        void MatrixMultiply(std::vector<T>& a, std::vector<T>& b, std::vector<T>& ret) {
            size_t size = this->random_variables_m.size();
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    for (int k = 0; k < size; k++) {
                        ret[i * size + j] += a[i * size + k] * b[k * size + j];
                    }
        }

        void Swap(atl::GradientStructure<T>* a, atl::GradientStructure<T>*b) {
            atl::GradientStructure<T>* temp = a;
            a = b;
            b = temp;
        }

        void ClearReStructures() {
            for (int i = 0; i < this->random_variables_m.size(); i++) {
                for (int j = 0; j < this->random_variables_m.size(); j++) {
                    //                    this->rand_hessian_dx[i][j] = 0.0;
                    this->rand_hessian[i][j] = 0.0;
                    //                    this->rand_hessian_inv[i][j] = 0.0;
                }
            }
        }

    public:

        void EvaluateLaplace(atl::Variable<T>& f) {
            bool recording = atl::Variable<T>::gradient_structure_g.recording;
            size_t PARAMETERS_SIZE = this->hyper_parameters_m.size();
            size_t RANDOM_SIZE = this->random_variables_m.size();
            size_t ALL_SIZE = this->all_variables_m.size();
            if (recording) {
                for (int i = 0; i < RANDOM_SIZE; i++) {
                    this->random_variables_m[i]->SetValue(0.0);
                }
                if (this->NewtonInner(20, 1e-7)) {
                    std::cout << "Inner converged!\n";
                    std::cout << "Inner f = " << this->inner_function_value << "\n";
                    std::cout << "Inner maxg = " << this->inner_maxgc << "\n\n";

                } else {
                    std::cout << "Inner failed!\n";
                    std::cout << "Inner f = " << this->inner_function_value << "\n";
                    std::cout << "Inner maxg = " << this->inner_maxgc << "\n\n";
                }
            }

            std::unordered_map<atl::VariableInfo<T>*, T> derivatives_logdet;
            std::unordered_map<atl::VariableInfo<T>*, T> derivatives_f;

            this->ClearReStructures();
            if (recording) {
                atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::THIRD_ORDER_MIXED_PARTIALS;
            } else {
                atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::SECOND_ORDER_MIXED_PARTIALS;
            }
            atl::Variable<T>::gradient_structure_g.recording = true;
            atl::Variable<T>::gradient_structure_g.Reset();
            atl::Variable<T> innerf2 = this->objective_function_m->Evaluate();
            atl::Variable<T>::gradient_structure_g.AccumulateFirstOrder();
            std::vector<T> g(PARAMETERS_SIZE);

            for (int j = 0; j < PARAMETERS_SIZE; j++) {
                g[j] = this->hyper_parameters_m[j]->info->dvalue;
            }

            atl::Variable<T>::gradient_structure_g.Reset();
            atl::Variable<T> innerf = this->objective_function_m->Evaluate();

            if (innerf.GetValue() != innerf.GetValue()) {
                std::cout << "Inner minimizer returned NaN, bailing out!" << std::endl;
                exit(0);
            }

            innerf.SetName("innerf");
            if (recording && re_active) {
                std::cout << "accumulating third-order mixed derivatives..." << std::flush;
                atl::Variable<T>::gradient_structure_g.AccumulateThirdOrderMixed();
                std::cout << "done\n" << std::flush;
            } else {
                atl::Variable<T>::gradient_structure_g.AccumulateSecondOrderMixed();
            }




            int ii = 0;
            for (int i = 0; i < RANDOM_SIZE; i++) {
                for (int j = 0; j < RANDOM_SIZE; j++) {
                    rand_hessian[i][j] = atl::Variable<T>::gradient_structure_g.Value(this->random_variables_m[i]->info->id, this->random_variables_m[j]->info->id);

                }
            }

            this->outer_iteration ? cholesky_outer.compute(rand_hessian) : cholesky_outer.recompute(rand_hessian);




            //compute logdetH
            T ld = cholesky_outer.logdet();
            log_det = ld;
            log_det.SetName("log_det");



            if (recording) {

                for (int i = 0; i < PARAMETERS_SIZE; i++) {
                    derivatives_f[this->hyper_parameters_m[i]->info] = g[i]; //hr(0, i);
                }


                //compute derivatives of the logdetH using third-order mixed partials
               
                for (int p = 0; p < PARAMETERS_SIZE; p++) {

                    for (int i = 0; i < RANDOM_SIZE; i++) {
                        for (int j = 0; j < RANDOM_SIZE; j++) {
                            this->rand_hessian_dx[i][j] = atl::Variable<T>::gradient_structure_g.Value(this->random_variables_m[i]->info->id,
                                    this->random_variables_m[j]->info->id, this->hyper_parameters_m[p]->info->id);
                            ret[i][j] = 0.0;
                        }
                    }
                    //                    }

                    cholesky_outer.solve(rand_hessian_dx, ret);


                    T tr = 0;
                    for (int i = 0; i < RANDOM_SIZE; i++) {
                        tr += ret[i][i]; 
                    }

                    derivatives_logdet[this->hyper_parameters_m[p]->info] = tr;
                }
            }
            //reset derivative structure
            atl::Variable<T>::gradient_structure_g.Reset();
            atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::FIRST_ORDER;


            atl::Variable<T>::gradient_structure_g.recording = false;
            f = this->objective_function_m->Evaluate();
            f.SetName("F");
            atl::Variable<T>::gradient_structure_g.recording = recording;
            if (recording) {

                //push adjoint entry for log_det
                atl::StackEntry<T>& entry = atl::Variable<T>::gradient_structure_g.NextEntry();
                entry.w = log_det.info;
                typename std::unordered_map<atl::VariableInfo<T>*, T>::iterator it;
                for (it = derivatives_logdet.begin(); it != derivatives_logdet.end(); ++it) {
                    entry.ids.insert((*it).first);
                }
                typename atl::StackEntry<T>::id_itereator id_it;
                entry.first.resize(entry.ids.size());
                size_t in = 0;
                for (id_it = entry.ids.begin(); id_it != entry.ids.end(); ++id_it) {
                    T dx = derivatives_logdet[(*id_it)];
                    entry.first[in] = dx;
                    in++;
                }


                //                push adjoint entry for of
                atl::StackEntry<T>& entry2 = atl::Variable<T>::gradient_structure_g.NextEntry();
                entry2.w = f.info;
                for (it = derivatives_f.begin(); it != derivatives_f.end(); ++it) {
                    entry2.ids.insert((*it).first);
                }
                entry2.first.resize(entry.ids.size());
                in = 0;
                for (id_it = entry2.ids.begin(); id_it != entry2.ids.end(); ++id_it) {

                    T dx = derivatives_f[(*id_it)];
                    entry2.first[in] = dx;
                    in++;
                }
                atl::Variable<T>::gradient_structure_g.recording = true;
            }

            f += static_cast<T> (.5) * log_det;
            f -= static_cast<T> (.5)*(T(RANDOM_SIZE) * std::log((static_cast<T> (2.0 * M_PI))));

        }

        void ComputeGradient(std::vector<atl::Variable<T>* >&p,
                std::valarray<T>&g, T & maxgc) {
            g.resize(p.size());
            atl::Variable<T>::gradient_structure_g.AccumulateFirstOrder();
           
            for (int i = 0; i < g.size(); i++) {
                g[i] = p[i]->info->dvalue;
                if (i == 0) {
                    maxgc = std::fabs(g[i]);
                } else {
                    if (std::fabs(g[i]) > maxgc) {

                        maxgc = std::fabs(g[i]);
                    }
                }
            }
        }

        void Print() {
            std::cout << "F = " << this->function_value << "\n";
            std::cout << " -----------------------------------------------------------------------------------\n";
            std::cout << '|' << util::center("Name", 15) << '|' << util::center("Value", 12) << '|' << util::center("Gradient", 12) << '|'
                    << util::center("Name", 15) << '|' << util::center("Value", 12) << '|' << util::center("Gradient", 12) << '|' << std::endl;
            std::cout << " -----------------------------------------------------------------------------------\n";
            int i = 0;
            std::cout.precision(4);
            std::cout << std::scientific;
            for (i = 0; (i + 2)< this->hyper_parameters_m.size(); i += 2) {
                std::cout << '|' << util::center(this->hyper_parameters_m[i]->GetName(), 15) << '|';
                std::cout << std::setw(12) << this->hyper_parameters_m[i]->GetValue() << '|' << std::setw(12) << this->gradient[i];
                std::cout << '|' << util::center(this->hyper_parameters_m[i + 1]->GetName(), 15) << '|';
                std::cout << std::setw(12) << this->hyper_parameters_m[i + 1]->GetValue() << '|' << std::setw(12) << this->gradient[i + 1] << "|\n";

                //                }
            }

            for (; i< this->hyper_parameters_m.size(); i++) {
                std::cout << '|' << util::center(this->hyper_parameters_m[i]->GetName(), 15) << '|';
                std::cout << std::setw(12) << this->hyper_parameters_m[i]->GetValue() << '|' << std::setw(12) << this->gradient[i];
                std::cout << '|' << util::center("-", 15) << '|';
                std::cout << util::center("-", 12) << '|' << util::center("-", 12) << "|\n";
            }
            std::cout << " -----------------------------------------------------------------------------------\n";
            std::cout << "\n\n";
        }

        virtual bool Evaluate() =0;

        T abs_sum(const std::valarray<T>&array) {
            T sum = 0.0;
            for (int i = 0; i < array.size(); i++) {
                sum += std::fabs(array[i]);
            }
            return sum;
        }

        bool NewtonInner(int max_iter = 20, T tol = 1e-12) {

            atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::SECOND_ORDER_MIXED_PARTIALS;
            atl::Variable<T>::SetRecording(true);
            atl::Variable<T> fx;
            atl::Variable<T>::gradient_structure_g.Reset();

            int nops = this->random_variables_m.size();

            std::vector<T> gradient_(nops);
            std::vector<T> p(nops);
            std::vector<std::vector<T> > hessian_(nops, std::vector<T>(nops));
            std::vector<std::vector<T> > hessian_inv_(nops, std::vector<T>(nops));

            this->inner_x.resize(nops);
            this->inner_best.resize(nops);
            this->inner_gradient.resize(nops);


            for (int i = 0; i < nops; i++) {
                if (this->random_variables_m[i]->IsBounded()) {
                    this->inner_x[i] = this->random_variables_m[i]->GetInternalValue();
                } else {
                    this->inner_x[i] = this->random_variables_m[i]->GetValue();
                }
                this->inner_gradient[i] = 0;
            }
            Cholesky lu;



            for (int iter = 0; iter < max_iter; iter++) {


                std::cout << "Newton raphson " << iter << std::endl;
                atl::Variable<T>::gradient_structure_g.Reset();
                fx = this->objective_function_m->Evaluate();
                this->inner_function_value = fx.GetValue();
                atl::Variable<T>::ComputeGradientAndHessian(
                        atl::Variable<T>::gradient_structure_g,
                        this->random_variables_m, gradient_, hessian_);


                //                Cholesky lu(hessian_);
                iter ? lu.compute(hessian_) : lu.recompute(hessian_);

                bool all_positive = true;
                this->inner_maxgc = std::fabs(gradient_[0]);

                for (int i = 0; i < gradient_.size(); i++) {
                    if (std::fabs(gradient_[i]) > this->inner_maxgc) {
                        this->inner_maxgc = std::fabs(gradient_[i]);
                    }

                }
                if (this->inner_maxgc <= tol /*&& all_positive*/) {
                    return true;
                }
                std::fill(p.begin(), p.end(), 0.0);
                lu.solve(gradient_, p);


                for (int j = 0; j < random_variables_m.size(); j++) {
                    random_variables_m[j]->SetValue(random_variables_m[j]->GetValue() - p[j]);
                }
            }


            return false;

        }

    
        bool line_search(std::vector<atl::Variable<T>* >& parameters,
                T& function_value,
                std::valarray<T>& x,
                std::valarray<T>& best,
                std::valarray<T>& z,
                std::valarray<T>& gradient,
                std::valarray<T>& wg,
                T& maxgc, int& i,
                bool inner = true) {
            T descent = 0;

            int nops = parameters.size();
            std::valarray<T> nwg(nops);
            std::valarray<T> ng(nops);

            for (size_t j = 0; j < nops; j++) {
                descent += z[j] * wg[j];
            }//end for

            T norm_g = this->norm(gradient);
            T relative_tolerance = this->tolerance * std::max<T > (T(1.0), norm_g);

            descent *= T(-1.0); // * Dot(z, g);
            if ((descent > T(-0.00000001) * relative_tolerance /* tolerance relative_tolerance*/)) {
                z = wg + .001;
                this->max_iterations -= i;
                i = 0;
                descent = -1.0 * Dot(z, wg);
            }//end if

            T step = i ? 1.0 : (1.0 / norm_g);

            if (step != step) {
                step = 1.0;
            }

            bool down = false;

            int ls;




            for (int j = 0; j < parameters.size(); j++) {
                best[j] = parameters[j]->GetValue();
            }

            atl::Variable<T> fx;
            for (ls = 0; ls < this->max_line_searches; ++ls) {



                // Tentative solution, gradient and loss
                std::valarray<T> nx = x - step * z;

                for (size_t j = 0; j < nops; j++) {

                    if (nx[j] != nx[j]) {
                    }
                    parameters[j]->UpdateValue(nx[j]);
                }

                //line_search:
                atl::Variable<T>::SetRecording(false);
                if (inner) {
                    this->CallInnerObjectiveFunction(fx);
                } else {
                    this->CallObjectiveFunction(fx);
                }
                //                atl::Variable<T>::gradient_structure_g.Reset();
                if (fx.GetValue() <= function_value + tolerance * T(10e-4) * step * descent) { // First Wolfe condition

                    for (size_t j = 0; j < nops; j++) {
                        best[j] = parameters[j]->GetInternalValue();
                    }


                    atl::Variable<T>::SetRecording(true);
                    atl::Variable<T>::gradient_structure_g.Reset();
                    atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::FIRST_ORDER;
                    if (inner) {
                        this->CallInnerObjectiveFunction(fx);
                    } else {
                        this->CallObjectiveFunction(fx);
                    }

                    this->ComputeGradient(parameters, ng, maxgc);
                    atl::Variable<T>::gradient_structure_g.Reset();
                    if (down || (-1.0 * Dot(z, nwg) >= 0.9 * descent)) { // Second Wolfe condition
                        x = nx;
                        gradient = ng;
                        function_value = fx.GetValue();
                        return true;
                    } else {
                        atl::Variable<T>::SetRecording(false);
                        step *= 2.0; //10.0;
                    }
                } else {
                    step *= .5; /// /= 10.0;
                    down = true;
                }
            }

            for (size_t j = 0; j < nops; j++) {
                parameters[j]->SetValue(best[j]);
            }

            return false;

        }

        /**
         * Compute the Norm of the vector v.
         *  
         * @param v
         * @return 
         */
        const T norm(std::valarray<T> &v) {

            T ret = (T) 0.0;
            unsigned int i;
            for (i = 0; i < v.size(); i++) {

                ret += v[i] * v[i];

            }
            return std::sqrt(ret);
        }

        /**
         * Compute the dot product of two vectors.
         * @param a
         * @param b
         * @return 
         */
        const T Dot(const std::valarray<T> &a, const std::valarray<T> &b) {
            T ret = 0;
            for (size_t i = 0; i < a.size(); i++) {

                ret += a[i] * b[i];
            }
            return ret;
        }

        /**
         * returns the a column of a matrix as a std::valarray.
         * @param matrix
         * @param column
         * @return 
         */
        const std::valarray<T> Column(std::valarray<std::valarray<T> > &matrix, size_t column, size_t length) {

            std::valarray<T> ret(length);

            for (int i = 0; i < ret.size(); i++) {

                ret[i] = matrix[i][column];
            }
            return ret;
        }

    };



    template<typename T>
    class LBFGS : public atl::OptimizationRoutine<T> {
    public:

        virtual bool Evaluate() {
            atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::FIRST_ORDER;
            atl::Variable<T>::SetRecording(true);
            atl::Variable<T>::gradient_structure_g.Reset();
            int nops = this->hyper_parameters_m.size();
            bool has_random_effects = false;
            if (this->random_variables_m.size()) {
                has_random_effects = true;
            }

            atl::Variable<T> pen;

            this->x.resize(nops);
            this->best.resize(nops);
            this->gradient.resize(nops);

            for (int i = 0; i < nops; i++) {
                if (this->hyper_parameters_m[i]->IsBounded()) {
                    this->x[i] = this->hyper_parameters_m[i]->GetInternalValue();
                } else {
                    this->x[i] = this->hyper_parameters_m[i]->GetValue();
                }
                this->gradient[i] = 0;
            }
            //
            //
            std::valarray<T> wg(nops);
            std::valarray<T> nwg(nops);
            std::valarray<T> ng(nops);


            //initial evaluation
            atl::Variable<T> fx(0.0);
            this->CallObjectiveFunction(fx);
            this->function_value = fx.GetValue();
            //
            //Historical evaluations
            std::valarray<T> px(nops);
            std::valarray<T> pg(nops);
            std::valarray<std::valarray<T> > dxs(std::valarray<T > (this->max_history), nops);
            std::valarray<std::valarray<T> > dgs(std::valarray<T > (this->max_history), nops);
            //search direction
            std::valarray<T> z(nops);

            this->ComputeGradient(this->hyper_parameters_m, this->gradient, this->maxgc);

            atl::Variable<T>::gradient_structure_g.Reset();
          

          
            std::valarray<T> p(this->max_history);
            std::valarray<T>a(this->max_history);

            int i;
            for (int iteration = 0; iteration < this->max_iterations; iteration++) {
                i = iteration;
                this->outer_iteration = iteration;
                for (int j = 0; j < nops; j++) {
                    wg[j] = this->hyper_parameters_m[j]->GetScaledGradient(this->hyper_parameters_m[j]->GetInternalValue()) * this->gradient[j];
                }

                if ((i % 10) == 0 || i == 0) {
                    std::cout << "Iteration " << i << "\n";
                    std::cout << "Phase = " << this->phase_m << "\n";

                    this->Print();
                }

                if (this->maxgc < this->tolerance) {
                    std::cout << "Iteration " << i << "\n";
                    std::cout << "Phase = " << this->phase_m << "\n";

                    this->Print();
                    return true;
                }

                z = wg;

                if (i > 0 && this->max_history > 0) {

                    size_t h = std::min<size_t > (i, this->max_history);
                    size_t end = (i - 1) % h;

                    //update histories
                    for (size_t r = 0; r < nops; r++) {
                        dxs[r][end] = this->hyper_parameters_m[r]->GetInternalValue() - px[r];
                        dgs[r][end] = wg[r] - pg[r];
                    }



                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end - j + h) % h;
                        p[k] = 1.0 / this->Dot(this->Column(dxs, k, this->hyper_parameters_m.size()), this->Column(dgs, k, this->hyper_parameters_m.size()));

                        a[k] = p[k] * this->Dot(this->Column(dxs, k, this->hyper_parameters_m.size()), z);
                        z -= a[k] * this->Column(dgs, k, this->hyper_parameters_m.size());
                    }
                    // Scaling of initial Hessian (identity matrix)
                    z *= this->Dot(this->Column(dxs, end, this->hyper_parameters_m.size()), this->Column(dgs, end, this->hyper_parameters_m.size())) / this->Dot(this->Column(dgs, end, this->hyper_parameters_m.size()), Column(dgs, end, this->hyper_parameters_m.size()));

                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end + j + 1) % h;
                        const T b = p[k] * Dot(this->Column(dgs, k, this->hyper_parameters_m.size()), z);
                        z += this->Column(dxs, k, this->hyper_parameters_m.size()) * (a[k] - b);
                    }

                }//end if(i>0)

                for (size_t j = 0; j < nops; j++) {
                    px[j] = this->hyper_parameters_m[j]->GetInternalValue();
                    this->x[j] = px[j];
                    pg[j] = wg[j];


                }//end for




                if (!this->line_search(this->hyper_parameters_m,
                        this->function_value,
                        this->x,
                        this->best,
                        z,
                        this->gradient,
                        wg,
                        this->maxgc,
                        iteration, false)) {
                    std::cout << "Outer Max line searches (" << this->max_line_searches << ").";
                    return false;

                }

            }

            std::cout << "Outer Max iterations!";

            return false;
        }
    private:

        /**
         * Compute the dot product of two vectors.
         * @param a
         * @param b
         * @return 
         */
        const T Dot(const std::valarray<T> &a, const std::valarray<T> &b) {
            T ret = 0;
            for (size_t i = 0; i < a.size(); i++) {

                ret += a[i] * b[i];
            }
            return ret;
        }

        /**
         * returns the a column of a matrix as a std::valarray.
         * @param matrix
         * @param column
         * @return 
         */
        const std::valarray<T> Column(std::valarray<std::valarray<T> > &matrix, size_t column, size_t length) {

            std::valarray<T> ret(length);

            for (int i = 0; i < ret.size(); i++) {

                ret[i] = matrix[i][column];
            }
            return ret;
        }


    };

    template<class T>
    class MCMC {
        ObjectiveFunction<T>* objective_function_m = NULL;
        std::vector<atl::Variable<T>* > hyper_parameters_m;
        std::vector<atl::Variable<T>* > random_variables_m;
        std::vector<atl::Variable<T>* > all_variables_m;
    public:

        MCMC(ObjectiveFunction<T>* objective_function = NULL) :
        objective_function_m(objective_function) {
        }

        ObjectiveFunction<T>* GetObjectiveFunction() const {

            return objective_function_m;
        }

        void SetObjectiveFunction(ObjectiveFunction<T>* objective_function) {

            this->objective_function_m = objective_function;
        }

        bool Run() {
            if (this->objective_function_m == NULL) {
                return false;
            }
            bool success = false;
            for (int i = 1; i <= this->objective_function_m->max_phase_m; i++) {

                this->Prepare(i);
                success = this->MetropoliHastings();
            }
            return success;
        }

        const MinimizerResults<T> GetResults() {

        }

    private:

        void Prepare(int phase) {
            this->hyper_parameters_m.resize(0);
            this->random_variables_m.resize(0);

            for (int i = 0; i < this->objective_function_m->hyper_parameter_phases_m.size(); i++) {
                if (this->objective_function_m->hyper_parameter_phases_m[i] <= phase) {
                    this->hyper_parameters_m.push_back(this->objective_function_m->hyper_parameters_m[i]);
                    this->all_variables_m.push_back(this->objective_function_m->hyper_parameters_m[i]);
                }
            }

            for (int i = 0; i < this->objective_function_m->random_variable_phases_m.size(); i++) {
                if (this->objective_function_m->random_variable_phases_m[i] <= phase) {

                    this->hyper_parameters_m.push_back(this->objective_function_m->random_variables_m[i]);
                    this->all_variables_m.push_back(this->objective_function_m->random_variables_m[i]);
                }
            }

        }

        bool MetropoliHastings() {
            //minimizer wrt hyper-parameters
        }

        bool LBGFS_Inner() {
            //minimize wrt random variables
        }

    };



}


#endif /* FUNCTIONMINIMIZER_HPP */
