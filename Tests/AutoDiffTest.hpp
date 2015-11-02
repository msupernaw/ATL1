/* 
 * File:   AutoDiffTests.hpp
 * Author: matthewsupernaw
 *
 * Created on June 15, 2015, 8:08 AM
 */

#ifndef HESSIANTESTS_HPP
#define	HESSIANTESTS_HPP

#include <valarray>
#include <fstream>
#include <sstream>
#include <iostream>
#include "../ATL.hpp"
namespace atl {
    namespace tests {
        namespace auto_diff {

            int fail = 0;
            int tests = 0;

            template<class T>
            class AutoDiffTest {
                typedef atl::Variable<T> var;
                std::vector<var*> active_parameters_m;
                std::vector<uint32_t> pids;
                static T tolerance;

            public:

                operator int() {
                    return fail;
                }

                AutoDiffTest() {
                }

                AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                static T GetTolerance() {
                    return tolerance;
                }

                static void SetTolerance(T tolerance) {
                    AutoDiffTest::tolerance = tolerance;
                }

                void Register(var& v) {
                    this->active_parameters_m.push_back(&v);
                    this->pids.push_back(v.info->id);
                }

                virtual void Initialize() {

                }

                virtual void ObjectiveFunction(atl::Variable<T>& f) {

                }

                void RunTest(bool exact_only = false) {
                    var::gradient_structure_g.Reset();
                    var::SetRecording(true);
                    var::derivative_trace_level = GRADIENT_AND_HESSIAN;
                    this->Initialize();
                    //                std::cout.precision(50);

                    std::cout << "Hessian Test:\n";
                    std::cout << "Number of Paramters: " << this->active_parameters_m.size() << "\n";


                    var f;
                    std::cout << "evaluating..." << std::flush;
                    auto exact_start = std::chrono::steady_clock::now();

                    auto eval_start = std::chrono::steady_clock::now();
                    this->ObjectiveFunction(f);
                    auto eval_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> eval_time = eval_end - eval_start;
                    std::cout << (eval_time.count()) << " sec...";

                    std::cout << "computing exact gradient and hessian..." << std::flush;

                    var::gradient_structure_g.HessianAndGradientAccumulate();
                    auto exact_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> exact_time = (exact_end - exact_start);
                    std::cout << exact_time.count() << "sec...";
                    std::vector<std::vector<T> > exact_hessian(this->active_parameters_m.size(), std::vector<T> (this->active_parameters_m.size()));
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            exact_hessian[i][j] = this->active_parameters_m[i]->info->hessian_row[this->active_parameters_m[j]->info];
                            //                        std::cout<<exact_hessian[i][j]<<" ";
                        }
                        //                    std::cout<<"\n";
                    }
                    //                std::vector<std::vector<T> > exact_hessian = var::gradient_structure_g.HessianAndGradientAccumulate(pids);
                    std::cout << "done!\n";
                    //     exit(0);
                    std::cout << std::scientific;
                    if (!exact_only) {


                        //                std::cout<<std::fixed;
                        std::cout << "estimating hessian..." << std::flush;
                        auto estimated_start = std::chrono::steady_clock::now();
                        std::valarray<std::valarray<T> > estimated_hessian = this->EstimatedHessian();
                        auto estimated_end = std::chrono::steady_clock::now();
                        std::chrono::duration<double> estimated_time = (estimated_end - estimated_start);

                        std::cout << "done\n";

                        if (this->active_parameters_m.size() < 20) {//dont print a lot of stuff
                            std::cout << "estimated: \n";
                            for (int i = 0; i < estimated_hessian.size(); i++) {
                                for (int j = 0; j < estimated_hessian[0].size(); j++) {
                                    std::cout << estimated_hessian[i][j] << " ";
                                }
                                std::cout << std::endl;
                            }
                        }
                        std::cout << "estimated done\n";
do_exact:

                        double sum = 0;
                        //
                        if (this->active_parameters_m.size() < 20) {//dont print a lot of stuff
                            std::cout << "Exact:\n";
                        }
                        for (int i = 0; i < exact_hessian.size(); i++) {
                            for (int j = 0; j < exact_hessian[0].size(); j++) {
                                if (this->active_parameters_m.size() < 20) {//dont print a lot of stuff
                                    std::cout << exact_hessian[i][j] << " ";
                                }
                                sum += (exact_hessian[i][j] - estimated_hessian[i][j])* (exact_hessian[i][j] - estimated_hessian[i][j]);
                            }
                            if (this->active_parameters_m.size() < 20) {
                                std::cout << std::endl;
                            }
                        }

                        if (this->active_parameters_m.size() < 20) {//dont print a lot of stuff
                            std::cout << "\n\nDifference:\n";
                            for (int i = 0; i < exact_hessian.size(); i++) {
                                for (int j = 0; j < exact_hessian[0].size(); j++) {
                                    std::cout << (exact_hessian[i][j] - estimated_hessian[i][j]) << " ";
                                }
                                std::cout << std::endl;
                            }
                        }

                        T mss = (sum / exact_hessian.size() * exact_hessian.size());

                        if (mss <= AutoDiffTest<T>::tolerance) {
                            std::cout << "Test Passed.\n";
                        } else {
                            std::cout << "Test Failed.\n";
                        }
                        std::cout << "Mean squared error = " << mss << "\n";
                        std::cout << std::fixed;
                        std::cout << "Time to compute exact Hessian: " << exact_time.count() << " sec\n";
                        std::cout << "Time to estimate Hessian: " << estimated_time.count() << " sec\n";
                        std::cout << "Speed up = " << estimated_time / exact_time << "";
                    }
                }

                virtual void Description(std::stringstream& out) {

                }

                void RunTestToFile(std::ofstream& out) {

                    tests += 2;
                    var::gradient_structure_g.Reset();
                    var::SetRecording(true);
                    var::gradient_structure_g.derivative_trace_level = GRADIENT;
                    this->Initialize();
                    var f;

                    T hmin;
                    T hmax;
                    T gmin;
                    T gmax;

                    std::stringstream ss;
                    this->Description(ss);
                    out << "\n\n" << ss.str() << "\n";

                    std::cout << "\n\n" << ss.str() << "\n";
                    T fval;
                    //                std::cout.precision(50);
                    std::vector<T> gradient(this->active_parameters_m.size());

                    out << "Number of Parameters: " << this->active_parameters_m.size() << "\n";


                    std::cout << "evaluating..." << std::flush;
                    ;


                    auto eval_gstart = std::chrono::steady_clock::now();
                    this->ObjectiveFunction(f);
                    auto eval_gend = std::chrono::steady_clock::now();
                    std::chrono::duration<double> eval_gtime = eval_gend - eval_gstart;

                    fval = f.GetValue();

                    std::cout << "computing exact gradient..." << std::flush;
                    auto exact_gstart = std::chrono::steady_clock::now();
                    var::gradient_structure_g.Accumulate();
                    auto exact_gend = std::chrono::steady_clock::now();
                    std::chrono::duration<double> exact_gtime = (exact_gend - exact_gstart);
                    var::gradient_structure_g.Reset();
                    std::cout << "done!\nevaluating..." << std::flush;
                    ;
                    f = 0.0;
                    var::gradient_structure_g.derivative_trace_level = GRADIENT_AND_HESSIAN;

                    //run hessian twice to make sure everything resets properly
                    this->ObjectiveFunction(f);
                    var::gradient_structure_g.HessianAndGradientAccumulate();
                    var::gradient_structure_g.Reset();

                    f = 0.0;
                    auto eval_start = std::chrono::steady_clock::now();
                    this->ObjectiveFunction(f);
                    auto eval_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> eval_time = eval_end - eval_start;

                    fval = f.GetValue();

                    std::cout << "computing exact hessian..." << std::flush;
                    auto exact_start = std::chrono::steady_clock::now();
                    var::gradient_structure_g.HessianAndGradientAccumulate();
                    auto exact_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> exact_time = (exact_end - exact_start);
                    std::vector<std::vector<T> > exact_hessian(this->active_parameters_m.size(), std::vector<T> (this->active_parameters_m.size()));
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        gradient[i] = this->active_parameters_m[i]->info->dvalue;
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            exact_hessian[i][j] = this->active_parameters_m[i]->info->hessian_row[this->active_parameters_m[j]->info->id];
                        }
                    }
                    std::cout << "done!\n";
                    std::cout << std::scientific;

                    var::gradient_structure_g.Reset();
                    std::cout << "estimating hessian..." << std::flush;
                    auto estimated_start = std::chrono::steady_clock::now();
                    std::valarray<std::valarray<T> > estimated_hessian = this->EstimatedHessian();
                    auto estimated_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> estimated_time = (estimated_end - estimated_start);

                    std::cout << "done\n";
                    T sum = 0;

                    for (int i = 0; i < exact_hessian.size(); i++) {
                        for (int j = 0; j < exact_hessian[0].size(); j++) {
                            T diff = (exact_hessian[i][j] - estimated_hessian[i][j]);
                            if (i == 0 && j == 0) {
                                hmax = diff;
                                hmin = diff;
                            } else {
                                if (diff > hmax) {
                                    hmax = diff;
                                }
                                if (diff < hmin) {
                                    hmin = diff;
                                }
                            }

                            sum += (exact_hessian[i][j] - estimated_hessian[i][j])* (exact_hessian[i][j] - estimated_hessian[i][j]);
                        }
                    }

                    T mse = sum / (exact_hessian.size() * exact_hessian.size());

                    var::gradient_structure_g.Reset();
                    auto estimatedg_start = std::chrono::steady_clock::now();
                    std::vector<T> estimated_gradient = this->EstimateGradient();
                    auto estimatedg_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> estimatedg_time = (estimatedg_end - estimatedg_start);

                    T gsum = 0;
                    for (int i = 0; i < gradient.size(); i++) {
                        T diff = (gradient[i] - estimated_gradient[i]);
                        if (i == 0) {
                            gmax = diff;
                            gmin = diff;
                        } else {
                            if (diff > gmax) {
                                gmax = diff;
                            }
                            if (diff < gmin) {
                                gmin = diff;
                            }
                        }

                        gsum += (gradient[i] - estimated_gradient[i])*(gradient[i] - estimated_gradient[i]);
                    }
                    T gmse = gsum / gradient.size();


                    if (gmse <= AutoDiffTest<T>::tolerance) {
                        std::cout << "Gradient Test Passed.\n";
                        out << "Gradient  Test Passed.\n";
                    } else {
                        fail++;
                        std::cout << "Gradient  Test Failed(" << gmse << ">" << AutoDiffTest<T>::tolerance << ")\n";
                        out << "Gradient Test Failed.\n";
                    }

                    if (mse <= AutoDiffTest<T>::tolerance) {
                        std::cout << "Hessian Test Passed.\n";
                        out << "Hessian  Test Passed.\n";
                    } else {
                        fail++;
                        std::cout << "Hessian  Test Failed(" << mse << ">" << AutoDiffTest<T>::tolerance << ")\n";
                        out << "Hessian Test Failed.\n";
                    }
                    out << "Function value: " << fval << "\n";
                    out << "Gradient mse = " << gmse << ", Error Range{" << gmin << " - " << gmax << "}\n";
                    out << "Hessian mse = " << mse << ", Error Range{" << hmin << " - " << hmax << "}\n";
                    out << std::fixed;
                    out << "Time to evaluate objective function with only gradient info: " << eval_gtime.count() << " sec\n";
                    out << "Time to evaluate objective function with Hessian info: " << eval_time.count() << " sec\n";
                    out << "Time to compute exact gradient: " << exact_gtime.count() << " sec\n";
                    out << "Time to compute exact gradient and Hessian: " << exact_time.count() << " sec\n";
                    out << "Time to estimate gradient: " << estimatedg_time.count() << " sec\n";
                    out << "Time to estimate Hessian: " << estimated_time.count() << " sec\n";
                    out << "Gradient Speed up = " << estimatedg_time / exact_gtime << "\n";
                    out << "Hessain Speed up = " << estimated_time / exact_time << "\n\n";



                    out << "estimated gradient:\n";
                    for (int i = 0; i < gradient.size(); i++) {
                        out << estimated_gradient[i] << " ";
                    }
                    out << "\n\n";

                    out << "exact gradient:\n";
                    for (int i = 0; i < gradient.size(); i++) {
                        out << gradient[i] << " ";
                    }
                    out << "\n\n";


                    out << "gradient difference:\n";
                    for (int i = 0; i < gradient.size(); i++) {
                        out << gradient[i] - estimated_gradient[i] << " ";
                    }
                    out << "\n\n";

                    out << "estimated Hessian: \n";
                    for (int i = 0; i < estimated_hessian.size(); i++) {
                        for (int j = 0; j < estimated_hessian[0].size(); j++) {
                            out << estimated_hessian[i][j] << " ";
                        }
                        out << std::endl;
                    }


                    //
                    out << "\nexact Hessian:\n";
                    for (int i = 0; i < exact_hessian.size(); i++) {
                        for (int j = 0; j < exact_hessian[0].size(); j++) {
                            out << exact_hessian[i][j] << " ";
                        }
                        out << std::endl;

                    }

                    out << "\nHessian difference:\n";
                    for (int i = 0; i < exact_hessian.size(); i++) {
                        for (int j = 0; j < exact_hessian[0].size(); j++) {
                            out << (exact_hessian[i][j] - estimated_hessian[i][j]) << " ";
                        }
                        out << std::endl;
                    }




                }

                const std::vector<T> EstimateGradient() {
                    var::gradient_structure_g.Reset();
                    atl::Variable<T>::SetRecording(false);
                    std::vector<T> gradient(active_parameters_m.size());
                    atl::Variable<T> f;
                    this->ObjectiveFunction(f);
                    atl::Variable<T> fh;
                    T fv;
                    T delta = 1.e-8;
                    for (int i = 0; i < active_parameters_m.size(); i++) {
                        active_parameters_m[i]->SetValue(active_parameters_m[i]->GetValue() + delta);
                        this->ObjectiveFunction(fh);
                        fv = fh.GetValue();
                        active_parameters_m[i]->SetValue(active_parameters_m[i]->GetValue() - 2 * delta);
                        this->ObjectiveFunction(fh);
                        gradient[i] = (fv - fh.GetValue()) / (2.0 * delta);
                        active_parameters_m[i]->SetValue(active_parameters_m[i]->GetValue() + delta);
                    }
                    return gradient;
                }

                /**
                 * This function is a port of from admb. 
                 * @return 
                 */
                const std::valarray<std::valarray<T> > EstimatedHessian() {

                    atl::Variable<T>::SetRecording(true);
                    atl::Variable<T>::gradient_structure_g.derivative_trace_level = GRADIENT;
                    atl::Variable<T>::gradient_structure_g.Reset();
                    std::valarray<std::valarray<T> > hessian(
                            std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    T delta = 1.e-5;
                    std::valarray < T> g1(active_parameters_m.size());
                    std::valarray < T> g2(active_parameters_m.size());
                    //            std::valarray < T> hess(active_parameters_m.size());
                    std::valarray < T> hess1(active_parameters_m.size());
                    std::valarray < T> hess2(active_parameters_m.size());
                    T eps = .1;
                    T sdelta1;
                    T sdelta2;

                    atl::Variable<T> f;

                    atl::Variable<T>::gradient_structure_g.Reset();
                    for (int i = 0; i < active_parameters_m.size(); i++) {



                        T xsave = active_parameters_m[i]->GetValue();
                        sdelta1 = /*active_parameters_m[i]->GetValue()*/ +delta;
                        //                sdelta1 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta1);
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();
                        //                Gradient(f, active_parameters_m, g1);

                        //        CallGradient(f, parameters_m, g1);
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g1[j] = active_parameters_m[j]->info->dvalue;
                        }

                        atl::Variable<T>::gradient_structure_g.Reset();

                        f = 0.0;
                        sdelta2 = /*active_parameters_m[i]->GetValue()*/ -delta;
                        //                sdelta2 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta2);
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();

                        //                Gradient(f, active_parameters_m, g1);

                        //        CallGradient(f, parameters_m, g2);
                        //                ObjectiveFunction(f);
                        //                //                Gradient(f, active_parameters_m, g2);
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g2[j] = active_parameters_m[j]->info->dvalue;
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();

                        active_parameters_m[i]->SetValue(xsave);

                        hess1 = (g1 - g2) / (sdelta1 - sdelta2);


                        sdelta1 = /*active_parameters_m[i]->GetValue() +*/ eps*delta;
                        //                sdelta1 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta1);
                        f = 0.0;
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();

                        //                Gradient(f, active_parameters_m, g1);

                        //        CallGradient(f, parameters_m, g1);
                        //                ObjectiveFunction(f);
                        //                //                Gradient(f, active_parameters_m, g1);
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g1[j] = active_parameters_m[j]->info->dvalue;
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();



                        active_parameters_m[i]->SetValue(xsave - eps * delta);
                        sdelta2 = /*active_parameters_m[i]->GetValue()*/ -eps*delta;
                        //                sdelta2 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta2);
                        f = 0.0;
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();

                        //                Gradient(f, active_parameters_m, g1);

                        //        CallGradient(f, parameters_m, g2);
                        //                ObjectiveFunction(f);
                        //   
                        //                Gradient(f, active_parameters_m, g2);
                        //                    std::cout << "g2 = ";
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g2[j] = active_parameters_m[j]->info->dvalue;
                            //                        std::cout << g2[j] << " ";
                        }
                        //                    std::cout << "\n";
                        atl::Variable<T>::gradient_structure_g.Reset();

                        active_parameters_m[i]->SetValue(xsave);


                        T eps2 = eps*eps;
                        hess2 = (g1 - g2) / (sdelta1 - sdelta2);
                        hessian[i] = (eps2 * hess1 - hess2) / (eps2 - 1.);
                    }

                    return hessian;

                }




            };

            template<class T>
            T AutoDiffTest<T>::tolerance = T(1e-5);

            template<class T>
            class AddAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                AddAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a + b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a + b;
                }


            };

            template<class T>
            class Add1AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Add1AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a + 3.1459" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a + 3.1459;
                }


            };

            template<class T>
            class Add2AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Add2AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = 3.1459 + b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = 3.1459 + b;
                }


            };

            template<class T>
            class SubtractAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                SubtractAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a - b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a - b;
                }


            };

            template<class T>
            class Subtract1AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Subtract1AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a - 3.1459" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a - 3.1459;
                }


            };

            template<class T>
            class Subtract2AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Subtract2AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = 3.1459 - b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = 3.1459 - b;
                }


            };

            template<class T>
            class MultiplyAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                MultiplyAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a * b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a * b;
                }


            };

            template<class T>
            class Multiply1AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Multiply1AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a * 3.1459" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a * 3.1459;
                }


            };

            template<class T>
            class Multiply2AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Multiply2AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = 3.1459 * b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = 3.1459 * b;
                }


            };

            template<class T>
            class DivideAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                DivideAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a / b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a / b;
                }


            };

            template<class T>
            class Divide1AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Divide1AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a / 3.1459" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = a / 3.1459;
                }


            };

            template<class T>
            class Divide2AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Divide2AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = 3.1459 / b" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = 3.1459 / b;
                }


            };

            template<class T>
            class ACosAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                ACosAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::acos((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::acos((a * b));
                }


            };

            template<class T>
            class ASinAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                ASinAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::asin((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::asin((a * b));
                }


            };

            template<class T>
            class ATanAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                ATanAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::atan((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::atan((a * b));
                }


            };

            template<class T>
            class CeilAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                CeilAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::ceil((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::ceil((a * b));
                }


            };

            template<class T>
            class CosAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                CosAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::cos((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::cos((a * b));
                }


            };

            template<class T>
            class CoshAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                CoshAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::cosh((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::cosh((a * b));
                }


            };

            template<class T>
            class ExpAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                ExpAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::exp((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::exp((a * b));
                }


            };

            template<class T>
            class FabsAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                FabsAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::fabs((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::fabs((a * b));
                }


            };

            template<class T>
            class FloorAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                FloorAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::floor((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::floor((a * b));
                }


            };

            template<class T>
            class LogAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                LogAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::log((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::log((a * b));
                }


            };

            template<class T>
            class Log10AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                Log10AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::log10((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::log10((a * b));
                }


            };

            template<class T>
            class PowAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                PowAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::pow((a * b),(a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::pow((a * b), (a * b));
                }


            };

            template<class T>
            class PowCAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                PowCAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::pow((a * b),2.0)" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::pow((a * b), 2.0);
                }


            };

            template<class T>
            class PowC2AutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                PowC2AutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::pow(2.0,(a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::pow(2.0, (a * b));
                }


            };

            template<class T>
            class SinAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                SinAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::sin((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::sin((a * b));
                }


            };

            template<class T>
            class SinhAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                SinhAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::sinh((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::sinh((a * b));
                }


            };

            template<class T>
            class SqrtAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                SqrtAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::sqrt((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::sqrt((a * b));
                }


            };

            template<class T>
            class TanAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                TanAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::tan((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::tan((a * b));
                }


            };

            template<class T>
            class TanhAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;

                TanhAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);



                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::tanh((a * b))" << std::endl;

                }

                void ObjectiveFunction(var& f) {
                    f = atl::tanh((a * b));
                }


            };

            template<class T>
            class SumAlotOfParametersAutoDiffTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;
                int nobs;
                std::vector<var> v;

                SumAlotOfParametersAutoDiffTest(std::ofstream& out, int n) : nobs(n) {
                    this->RunTestToFile(out);
                }

                void Initialize() {
                    a = .03434;
                    b = .034230;
                    this->Register(a);
                    this->Register(b);
                    //                int nobs = 50;

                    T small = 1.00001212;
                    for (int i = 0; i < nobs; i++) {
                        small += .000000001;
                        v.push_back(var(small));
                        //                 
                    }
                    //
                    for (int i = 0; i < nobs; i++) {
                        this->Register(v[i]);
                    }

                }

                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << ";\n";
                    out << "Variable b = " << b << ";\n";
                    out << "Variable v[" << nobs << "]\n\n";
                    out << "Variable ff = 0;\n";
                    out << "ff = (a / b); \n"; //// (a * b+a/b+a*a) - (a*(a * b+a/b+a*a));
                    out << "for (int i = 0; i < nobs; i++) {\n";
                    out << "    for (int j = 0; j < nobs; j++) {\n";
                    out << "         ff += atl::log(v[j] * v[i] * atl::tanh(a / b * v[j] * v[i]) * atl::tanh(a / b) * v[j] * v[i] * atl::tanh(a / b));\n";
                    out << "     }\n";
                    out << "}\n";
                    out<<"f.assign(Variable::gradient_structure_g, ff);\n";

                }

                void ObjectiveFunction(var& f) {
                    var ff;
                    ff = (a / b); // (a * b+a/b+a*a) - (a*(a * b+a/b+a*a));
                    for (int i = 0; i < nobs; i++) {
                        for (int j = 0; j < nobs; j++) {
                            ff += (v[j] * v[i] * atl::tanh(a / b * v[j] * v[i]) * atl::tanh(a / b) * v[j] * v[i] * atl::tanh(a / b));
                        }
                    }
                    ff = (nobs / 2.0) * atl::log(ff);
                    f.Assign(var::gradient_structure_g, ff);
                }


            };

            void Run() {
                std::ofstream out;
                out.open("autodiff_tests.txt");
                out << "This test compares the computed gradient and Hessian matrices against\n"
                        "those computed using the differences method.\n\n";

                out << "Test Tolerance: " << atl::tests::auto_diff::AutoDiffTest<double>::GetTolerance() << "\n\n";

                std::cout << "running...\n";

                atl::tests::auto_diff::AddAutoDiffTest<double> add(out);
                atl::tests::auto_diff::Add1AutoDiffTest<double> add1(out);
                atl::tests::auto_diff::Add2AutoDiffTest<double> add2(out);
                atl::tests::auto_diff::SubtractAutoDiffTest<double> subtract1(out);
                atl::tests::auto_diff::Subtract1AutoDiffTest<double> subtract(out);
                atl::tests::auto_diff::Subtract2AutoDiffTest<double> subtract2(out);
                atl::tests::auto_diff::MultiplyAutoDiffTest<double> multiply(out);
                atl::tests::auto_diff::Multiply1AutoDiffTest<double> multiply1(out);
                atl::tests::auto_diff::Multiply2AutoDiffTest<double> multiply2(out);
                atl::tests::auto_diff::DivideAutoDiffTest<double> divide(out);
                atl::tests::auto_diff::Divide1AutoDiffTest<double> divide1(out);
                atl::tests::auto_diff::Divide2AutoDiffTest<double> divide2(out);
                atl::tests::auto_diff::CosAutoDiffTest<double> cos(out);
                atl::tests::auto_diff::ACosAutoDiffTest<double> acos(out);
                atl::tests::auto_diff::SinAutoDiffTest<double> sin(out);
                atl::tests::auto_diff::ASinAutoDiffTest<double> asin(out);
                atl::tests::auto_diff::TanAutoDiffTest<double> tan(out);
                atl::tests::auto_diff::ATanAutoDiffTest<double> atan(out);
                atl::tests::auto_diff::CoshAutoDiffTest<double> cosh(out);
                atl::tests::auto_diff::SinhAutoDiffTest<double> sinh(out);
                atl::tests::auto_diff::TanhAutoDiffTest<double> tanh(out);
                atl::tests::auto_diff::ExpAutoDiffTest<double> exp(out);
                atl::tests::auto_diff::LogAutoDiffTest<double> log(out);
                atl::tests::auto_diff::Log10AutoDiffTest<double> log10(out);
                atl::tests::auto_diff::FabsAutoDiffTest<double> fabs(out);
                atl::tests::auto_diff::SqrtAutoDiffTest<double> sqrt(out);
                atl::tests::auto_diff::PowAutoDiffTest<double> pow(out);
                atl::tests::auto_diff::PowCAutoDiffTest<double> pow2(out);
                atl::tests::auto_diff::PowC2AutoDiffTest<double> pow3(out);
                atl::tests::auto_diff::CeilAutoDiffTest<double> ceil(out);
                atl::tests::auto_diff::FloorAutoDiffTest<double> floor(out);
                atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<double> s1(out, 10);
                //                atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<double> s2(out, 50);
                //                            atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<double> s3(out, 100);
                //                            atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<double> s4(out, 200);
                //                        atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<double> s5(out, 500);
                std::cout << "Test complete.\n";
                if (atl::tests::auto_diff::fail == 0) {
                    std::cout << "All tests passed, review file \"autodiff_test.txt\" for details." << std::endl;

                } else {
                    std::cout << atl::tests::auto_diff::fail << " of " << atl::tests::auto_diff::tests << " tests did not agree within a tolerance of " << atl::tests::auto_diff::AutoDiffTest<double>::GetTolerance() << ", review file \"autodiff_test.txt\" for details." << std::endl;
                }

            }

        }
    }
}



#endif	/* HESSIANTESTS_HPP */