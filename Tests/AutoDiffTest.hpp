/* 
 * File:   AutoDiffTests.hpp
 * Author: matthewsupernaw
 *
 * Created on June 15, 2015, 8:08 AM
 */
/*
 * File:   AutoDiffTests.hpp
 * Author: matthewsupernaw
 *
 * Created on June 15, 2015, 8:08 AM
 */

#ifndef HESSIANTESTS_HPP
#define HESSIANTESTS_HPP

#include <valarray>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>
#include "../ATL.hpp"
#include "../Utilities/BigFloat.hpp"

#define HESSIAN_USE_AD_GRADIENT

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
                static T second_tolerance;
                static T third_tolerance;
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
                
                static T GetSecondOrderTolerance() {
                    return AutoDiffTest::second_tolerance;
                }
                
                static void SetSecondOrderTolerance(T tolerance) {
                    AutoDiffTest::second_tolerance = tolerance;
                }
                
                static T GetThirdOrderTolerance() {
                    return AutoDiffTest::third_tolerance;
                }
                
                static void SetThirdOrderTolerance(T tolerance) {
                    AutoDiffTest::third_tolerance = tolerance;
                }
                
                void Register(var& v) {
                    this->active_parameters_m.push_back(&v);
                    this->pids.push_back(v.info->id);
                }
                
                virtual void Initialize() {
                    
                }
                
                virtual void ObjectiveFunction(atl::Variable<T>& f) {
                    
                }
                
                virtual void Description(std::stringstream& out) {
                    
                }
                
                void RunTestToFile(std::ofstream& out) {
                    
                    tests += 2;
                    var::gradient_structure_g.Reset();
                    var::SetRecording(true);
                    var::gradient_structure_g.derivative_trace_level = FIRST_ORDER;
                    this->Initialize();
                    var f;
                    //
                    T hmin;
                    T hmax;
                    T gmin;
                    T gmax;
                    //                    //                    //
                    std::stringstream ss;
                    this->Description(ss);
                    out << "\n\n" << ss.str() << "\n";
                    
                    std::cout << "\n\n" << ss.str() << "\n";
                    T fval;
                    //                std::cout.precision(50);
                    std::vector<T> gradient(this->active_parameters_m.size());
                    
                    out << "Number of Parameters: " << this->active_parameters_m.size() << "\n";
                    
                    
                    std::cout << "evaluating..." << std::flush;
                    
                    
                    var::gradient_structure_g.Reset();
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
                    
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        gradient[i] = this->active_parameters_m[i]->info->dvalue;
                    }
                    
                    
                    var::gradient_structure_g.Reset();
                    std::cout << "done!\nevaluating..." << std::flush;
                    
                    f = 0.0;
                    var::gradient_structure_g.derivative_trace_level = atl::SECOND_ORDER_MIXED_PARTIALS;
                    
                    //                    //run hessian twice to make sure everything resets properly
                    //                    this->ObjectiveFunction(f);
                    //                    var::gradient_structure_g.HessianAndGradientAccumulate();
                    //                    var::gradient_structure_g.Reset();
                    
                    f = 0.0;
                    auto eval_start = std::chrono::steady_clock::now();
                    this->ObjectiveFunction(f);
                    auto eval_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> eval_time = eval_end - eval_start;
                    
                    fval = f.GetValue();
                    
                    std::cout << "computing exact hessian..." << std::flush;
                    auto exact_start = std::chrono::steady_clock::now();
                    var::gradient_structure_g.AccumulateSecondOrderMixed();
                    auto exact_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> exact_time = (exact_end - exact_start);
                    std::vector<std::vector<T> > exact_hessian(this->active_parameters_m.size(), std::vector<T> (this->active_parameters_m.size()));
                    std::vector<std::vector<T> > exact_hessian2(this->active_parameters_m.size(), std::vector<T> (this->active_parameters_m.size()));
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        //                        gradient[i] = this->active_parameters_m[i]->info->dvalue;
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
#ifdef USE_EIGEN
                            
                            exact_hessian[i][j] = var::gradient_structure_g.Value(this->active_parameters_m[i]->info->id, this->active_parameters_m[j]->info->id); //this->active_parameters_m[i]->info->hessian_row(this->active_parameters_m[j]->info->id);
#else
                            exact_hessian[i][j] = var::gradient_structure_g.Value(this->active_parameters_m[i]->info->id, this->active_parameters_m[j]->info->id); //this->active_parameters_m[i]->info->hessian_row(this->active_parameters_m[j]->info->id);this->active_parameters_m[i]->info->hessian_row[this->active_parameters_m[j]->info->id];
#endif
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
                    std::valarray<T> estimated_gradient = this->EstimateGradient();
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
                    
                    
                    
                    
                    atl::Variable<T> tf;
                    atl::Variable<T>::gradient_structure_g.recording = true;
                    atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::THIRD_ORDER_MIXED_PARTIALS;
                    atl::Variable<T>::gradient_structure_g.Reset();
                    std::cout << "evaluating..." << std::flush;
                    
                    auto eval_to_start = std::chrono::steady_clock::now();
                    this->ObjectiveFunction(tf);
                    auto eval_to_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> eval_to_time = (eval_to_end - eval_to_start);
                    
                    std::cout << "computing exact third order mixed..." << std::flush;
                    
                    auto exact_to_start = std::chrono::steady_clock::now();
                    atl::Variable<T>::gradient_structure_g.AccumulateThirdOrderMixed();
                    auto exact_to_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> exact_to_time = (exact_to_end - exact_to_start);
                    std::cout << "done.\n" << std::flush;
                    
                    std::valarray< std::valarray<std::valarray < T> > > third_mixed_exact(std::valarray<std::valarray < T> > (std::valarray<T > (active_parameters_m.size()), active_parameters_m.size()), active_parameters_m.size());
                    
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        if (std::fabs(gradient[i] - this->active_parameters_m[i]->info->dvalue) > 1e-5) {
                            std::cout << "gradient doesn't match....\n";
                            std::cout << gradient[i] << " != " << this->active_parameters_m[i]->info->dvalue;
                            //                            exit(0);
                        }
                        //                        gradient[i]  =  this->active_parameters_m[i]->info->dvalue;
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            for (int k = 0; k < this->active_parameters_m.size(); k++) {
                                if (i == 0) {
                                    exact_hessian2[j][k] = atl::Variable<T>::gradient_structure_g.Value(this->active_parameters_m[j]->info->id, this->active_parameters_m[k]->info->id);
                                    if (std::fabs(exact_hessian[j][k] - exact_hessian2[j][k]) > 1e-5) {
                                        
                                        std::cout << "hessian doesn't match....\n";
                                        std::cout << exact_hessian[j][k] << " != " << exact_hessian2[j][k];
                                        //                                        exit(0);
                                    }
                                }
                                third_mixed_exact[i][j][k] = atl::Variable<T>::gradient_structure_g.Value(this->active_parameters_m[i]->info->id, this->active_parameters_m[j]->info->id, this->active_parameters_m[k]->info->id);
                            }
                        }
                    }
                    
                    std::valarray< std::valarray<std::valarray < T> > > third_mixed_estimated(std::valarray<std::valarray < T> > (std::valarray<T > (active_parameters_m.size()), active_parameters_m.size()), active_parameters_m.size());
                    std::cout << "estimating third order mixed..." << std::flush;
                    
                    auto estimated_to_start = std::chrono::steady_clock::now();
                    third_mixed_estimated = this->EstimateThirdOrderMixed();
                    auto estimated_to_end = std::chrono::steady_clock::now();
                    std::chrono::duration<double> estimated_to_time = (estimated_to_end - estimated_to_start);
                    std::cout << "done.\n" << std::flush;
                    
                    T to_sum = 0;
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            for (int k = 0; k < this->active_parameters_m.size(); k++) {
                                to_sum += third_mixed_exact[i][j][k] - third_mixed_estimated[i][j][k];
                                //                                if (std::fabs(third_mixed_exact[i][j][k] - third_mixed_estimated[i][j][k]) > 1e-2) {
                                //                                    std::cout << "{" << this->active_parameters_m[i]->info->id << "," << this->active_parameters_m[j]->info->id << "," << this->active_parameters_m[k]->info->id << "}" << "third order failed " << third_mixed_exact[i][j][k] << " != " << third_mixed_estimated[i][j][k] << " <-----\n";
                                //                                }
                                //                                    else {
                                //                                        std::cout << "{" << this->active_parameters_m[i]->info->id << "," << this->active_parameters_m[j]->info->id << "," << this->active_parameters_m[k]->info->id << "}" << "third order PASSED " << third_mixed_exact[i][j][k] << " == " << third_mixed_estimated[i][j][k] << "\n";
                                //                                    }
                            }
                        }
                    }
                    
                    T to_mse = to_sum / (exact_hessian.size() * exact_hessian.size() * exact_hessian.size());
                    
                    
                    
                    if (gmse <= AutoDiffTest<T>::tolerance) {
                        std::cout << "Gradient Test Passed.\n";
                        out << "Gradient  Test Passed.\n";
                    } else {
                        fail++;
                        std::cout << "Gradient  Test Failed(" << gmse << ">" << AutoDiffTest<T>::tolerance << ")\n";
                        out << "Gradient Test Failed.\n";
                    }
                    
                    if (mse <= AutoDiffTest<T>::second_tolerance) {
                        std::cout << "Hessian Test Passed.\n";
                        out << "Hessian  Test Passed.\n";
                    } else {
                        fail++;
                        std::cout << "Hessian  Test Failed(" << mse << ">" << AutoDiffTest<T>::tolerance << ")\n";
                        out << "Hessian Test Failed.\n";
                    }
                    
                    
                    if (to_mse <= AutoDiffTest<T>::third_tolerance) {
                        std::cout << "Third-Order Test Passed.\n";
                        out << "Third-Order Test Passed.\n";
                    } else {
                        fail++;
                        std::cout << "Third-Order Test Failed(" << to_mse << ">" << AutoDiffTest<T>::tolerance << ")\n";
                        out << "Third-Order Test Failed.\n";
                    }
                    
                    out << "Function value: " << fval << "\n";
                    out << "Gradient mse = " << gmse << ", Error Range{" << gmin << " - " << gmax << "}\n";
                    out << "Hessian mse = " << mse << ", Error Range{" << hmin << " - " << hmax << "}\n";
                    out << "Third-Order mse = " << to_mse << "\n"; // ", Error Range{" << hmin << " - " << hmax << "}\n";
                    out << std::fixed;
                    out << "Time to evaluate objective function with only gradient info: " << eval_gtime.count() << " sec\n";
                    out << "Time to evaluate objective function with Hessian info: " << eval_time.count() << " sec\n";
                    out << "Time to evaluate objective function with Third-Order info: " << eval_to_time.count() << " sec\n";
                    out << "Time to compute exact gradient: " << exact_gtime.count() << " sec\n";
                    out << "Time to compute exact gradient and Hessian: " << exact_time.count() << " sec\n";
                    out << "Time to compute exact gradient, Hessian, and Third-Order: " << exact_to_time.count() << " sec\n";
                    out << "Time to estimate gradient: " << estimatedg_time.count() << " sec\n";
                    out << "Time to estimate Hessian: " << estimated_time.count() << " sec\n";
                    out << "Time to estimate Third-Order: " << estimated_to_time.count() << " sec\n";
                    out << "Gradient Speed up = " << estimatedg_time / exact_gtime << "\n";
                    out << "Hessain Speed up = " << estimated_time / exact_time << "\n\n";
                    out << "Third-Order Speed up = " << estimated_to_time / exact_to_time << "\n\n";
                    std::cout << std::fixed;
                    std::cout << "Gradient Speed up = " << estimatedg_time / exact_gtime << "\n";
                    std::cout << "Hessain Speed up = " << estimated_time / exact_time << "\n";
                    std::cout << "Third-Order Speed up = " << estimated_to_time / exact_to_time << "\n\n";
                    
                    out << "Estimated gradient:\n";
                    for (int i = 0; i < gradient.size(); i++) {
                        out << estimated_gradient[i] << " ";
                    }
                    out << "\n\n";
                    
                    out << "Exact gradient:\n";
                    for (int i = 0; i < gradient.size(); i++) {
                        out << gradient[i] << " ";
                    }
                    out << "\n\n";
                    
                    
                    out << "Gradient difference:\n";
                    for (int i = 0; i < gradient.size(); i++) {
                        out << gradient[i] - estimated_gradient[i] << " ";
                    }
                    out << "\n\n";
                    
                    out << "Estimated Hessian: \n";
                    for (int i = 0; i < estimated_hessian.size(); i++) {
                        for (int j = 0; j < estimated_hessian[0].size(); j++) {
                            out << estimated_hessian[i][j] << " ";
                        }
                        out << std::endl;
                    }
                    
                    
                    //
                    out << "\nExact Hessian:\n";
                    for (int i = 0; i < exact_hessian.size(); i++) {
                        for (int j = 0; j < exact_hessian[0].size(); j++) {
                            out << exact_hessian[i][j] << " ";
                        }
                        out << std::endl;
                        
                    }
                    
                    out << "\nHessian difference:\n";
                    for (int i = 0; i < exact_hessian.size(); i++) {
                        for (int j = 0; j < exact_hessian[0].size(); j++) {
                            T diff = (exact_hessian[i][j] - estimated_hessian[i][j]);
                            out << diff << " ";
                        }
                        out << std::endl;
                    }
                    
                    out << "\nEstimated Third-Order:\n";
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            out << "[" << i << "] ";
                            for (int k = 0; k < this->active_parameters_m.size(); k++) {
                                out << third_mixed_estimated[i][j][k] << " ";
                            }
                            out << "\n";
                        }
                    }
                    
                    out << "\nExact Third-Order:\n";
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            out << "[" << i << "] ";
                            for (int k = 0; k < this->active_parameters_m.size(); k++) {
                                out << third_mixed_exact[i][j][k] << " ";
                            }
                            out << "\n";
                        }
                    }
                    
                    
                    out << "\nThird-Order difference:\n";
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        
                        for (int j = 0; j < this->active_parameters_m.size(); j++) {
                            out << "[" << i << "] ";
                            for (int k = 0; k < this->active_parameters_m.size(); k++) {
                                out << (third_mixed_exact[i][j][k] - third_mixed_estimated[i][j][k]) << " ";
                            }
                            out << "\n";
                        }
                    }
                    
                    
                }
                
                const std::valarray<T> EstimateGradient() {
                    var::gradient_structure_g.Reset();
                    atl::Variable<T>::SetRecording(false);
                    std::valarray<T> gradient(active_parameters_m.size());
                    atl::Variable<T> f;
                    this->ObjectiveFunction(f);
                    atl::Variable<T> fh;
                    T fv;
                    T delta = 1.e-5;
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
                    T delta = 1.e-3;
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
#ifdef HESSIAN_USE_AD_GRADIENT
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g1[j] = active_parameters_m[j]->info->dvalue;
                        }
                        
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        g1 = this->EstimateGradient();
#endif
                        f = 0.0;
                        sdelta2 = /*active_parameters_m[i]->GetValue()*/ -delta;
                        //                sdelta2 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta2);
#ifdef HESSIAN_USE_AD_GRADIENT
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g2[j] = active_parameters_m[j]->info->dvalue;
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        g2 = this->EstimateGradient();
#endif
                        active_parameters_m[i]->SetValue(xsave);
                        
                        hess1 = (g1 - g2) / (sdelta1 - sdelta2);
                        
                        
                        sdelta1 = /*active_parameters_m[i]->GetValue() +*/ eps*delta;
                        //                sdelta1 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta1);
#ifdef HESSIAN_USE_AD_GRADIENT
                        f = 0.0;
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g1[j] = active_parameters_m[j]->info->dvalue;
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        g1 = this->EstimateGradient();
#endif
                        
                        active_parameters_m[i]->SetValue(xsave - eps * delta);
                        sdelta2 = /*active_parameters_m[i]->GetValue()*/ -eps*delta;
                        //                sdelta2 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta2);
#ifdef HESSIAN_USE_AD_GRADIENT
                        f = 0.0;
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.Accumulate();
                        
                        
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            g2[j] = active_parameters_m[j]->info->dvalue;
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        g2 = this->EstimateGradient();
#endif
                        active_parameters_m[i]->SetValue(xsave);
                        
                        
                        T eps2 = eps*eps;
                        hess2 = (g1 - g2) / (sdelta1 - sdelta2);
                        hessian[i] = (eps2 * hess1 - hess2) / (eps2 - 1.);
                    }
                    var::gradient_structure_g.Reset();
                    return hessian;
                    
                }
                
                const std::valarray<std::valarray<std::valarray<T> > > EstimateThirdOrderMixed() {
                    atl::Variable<T>::SetRecording(true);
                    atl::Variable<T>::gradient_structure_g.derivative_trace_level = atl::SECOND_ORDER_MIXED_PARTIALS;
                    atl::Variable<T>::gradient_structure_g.Reset();
                    
                    std::valarray< std::valarray<std::valarray < T> > > third_mixed(std::valarray<std::valarray < T> > (std::valarray<T > (active_parameters_m.size()), active_parameters_m.size()), active_parameters_m.size());
                    
                    T delta = 1.e-3;
                    T eps = .1;
                    T sdelta1;
                    T sdelta2;
                    atl::Variable<T> f;
                    
                    atl::Variable<T>::gradient_structure_g.Reset();
                    std::valarray < std::valarray < T> > h1(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T >> h2(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    //            std::valarray < T> hess(active_parameters_m.size());
                    std::valarray < std::valarray < T >> t1(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T >> t2(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    
                    std::valarray < std::valarray < T> > sd1(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T >> sd2(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T> > sd3(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T >> sd4(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T >> eps_squared(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    std::valarray < std::valarray < T >> eps_squared_minus_one(std::valarray<T > (active_parameters_m.size()), active_parameters_m.size());
                    for (int i = 0; i < active_parameters_m.size(); i++) {
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            sd1[i][j] = delta;
                            sd2[i][j] = -1.0 * delta;
                            sd3[i][j] = eps*delta;
                            sd4[i][j] = -1.0 * eps*delta;
                            eps_squared[i][j] = eps*eps;
                            eps_squared_minus_one[i][j] = eps_squared[i][j] - 1.0;
                        }
                    }
                    
                    
                    for (int i = 0; i < this->active_parameters_m.size(); i++) {
                        
                        
                        T xsave = active_parameters_m[i]->GetValue();
                        sdelta1 = /*active_parameters_m[i]->GetValue()*/ +delta;
                        //                sdelta1 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta1);
#ifdef HESSIAN_USE_AD_GRADIENT
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.AccumulateSecondOrderMixed();
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            for (int k = 0; k < active_parameters_m.size(); k++) {
                                h1[j][k] = atl::Variable<T>::gradient_structure_g.Value(active_parameters_m[j]->info->id, active_parameters_m[k]->info->id); //active_parameters_m[j]->info->hessian_row[active_parameters_m[k]->info->id];
                            }
                        }
                        
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        h1 = this->EstimatedHessian();
#endif
                        f = 0.0;
                        sdelta2 = /*active_parameters_m[i]->GetValue()*/ -delta;
                        //                sdelta2 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta2);
#ifdef HESSIAN_USE_AD_GRADIENT
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.AccumulateSecondOrderMixed();
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            for (int k = 0; k < active_parameters_m.size(); k++) {
                                h2[j][k] = atl::Variable<T>::gradient_structure_g.Value(active_parameters_m[j]->info->id, active_parameters_m[k]->info->id); //active_parameters_m[j]->info->hessian_row[active_parameters_m[k]->info->id];
                            }
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        h2 = this->EstimatedHessian();
#endif
                        active_parameters_m[i]->SetValue(xsave);
                        
                        
                        t1 = (h1 - h2) / (sd1 - sd2);
                        
                        
                        sdelta1 = /*active_parameters_m[i]->GetValue() +*/ eps*delta;
                        //                sdelta1 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta1);
#ifdef HESSIAN_USE_AD_GRADIENT
                        f = 0.0;
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.AccumulateSecondOrderMixed();
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            for (int k = 0; k < active_parameters_m.size(); k++) {
                                h1[j][k] = atl::Variable<T>::gradient_structure_g.Value(active_parameters_m[j]->info->id, active_parameters_m[k]->info->id); //active_parameters_m[j]->info->hessian_row[active_parameters_m[k]->info->id];
                            }
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        h1 = this->EstimatedHessian();
#endif
                        
                        active_parameters_m[i]->SetValue(xsave - eps * delta);
                        sdelta2 = /*active_parameters_m[i]->GetValue()*/ -eps*delta;
                        //                sdelta2 -= active_parameters_m[i]->GetValue();
                        active_parameters_m[i]->SetValue(xsave + sdelta2);
#ifdef HESSIAN_USE_AD_GRADIENT
                        f = 0.0;
                        ObjectiveFunction(f);
                        atl::Variable<T>::gradient_structure_g.AccumulateSecondOrderMixed();
                        
                        for (int j = 0; j < active_parameters_m.size(); j++) {
                            for (int k = 0; k < active_parameters_m.size(); k++) {
                                h2[j][k] = atl::Variable<T>::gradient_structure_g.Value(active_parameters_m[j]->info->id, active_parameters_m[k]->info->id); //active_parameters_m[j]->info->hessian_row[active_parameters_m[k]->info->id];
                            }
                        }
                        atl::Variable<T>::gradient_structure_g.Reset();
#else
                        h2 = this->EstimatedHessian();
#endif
                        active_parameters_m[i]->SetValue(xsave);
                        
                        
                        T eps2 = eps*eps;
                        t2 = (h1 - h2) / (sd3 - sd4);
                        third_mixed[i] = (eps_squared * t1 - t2) / (eps_squared_minus_one);
                    }
                    
                    
                    
                    
                    return third_mixed;
                    
                }
                
            };
            
            template<class T>
            T AutoDiffTest<T>::tolerance = T(1e-5);
            
            template<class T>
            T AutoDiffTest<T>::second_tolerance = T(1e-3);
            
            template<class T>
            T AutoDiffTest<T>::third_tolerance = T(1e-1);
            
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
                    //                    this->Register(b);
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a + 3.1459" << std::endl;
                    
                }
                
                void ObjectiveFunction(var& f) {
                    f = a + static_cast<T> (3.1459);
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
                    //                    this->Register(a);
                    this->Register(b);
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    //                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = 3.1459 + b" << std::endl;
                    
                }
                
                void ObjectiveFunction(var& f) {
                    f = static_cast<T> (3.1459) + b;
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
                    //                    this->Register(b);
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a - 3.1459" << std::endl;
                    
                }
                
                void ObjectiveFunction(var& f) {
                    f = a - static_cast<T> (3.1459);
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
                    //                    this->Register(a);
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
                    f = static_cast<T> (3.1459) - b;
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
                    //                    this->Register(b);
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a * 3.1459" << std::endl;
                    
                }
                
                void ObjectiveFunction(var& f) {
                    f = a * static_cast<T> (3.1459);
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
                    //                    this->Register(a);
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
                    f = static_cast<T> (3.1459) * b;
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
                    f = atl::sin(a) / b;
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
                    //                    this->Register(b);
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = a / 3.1459" << std::endl;
                    
                }
                
                void ObjectiveFunction(var& f) {
                    f = a / static_cast<T> (3.1459);
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
                    a = 1.03434;
                    b = 1.034230;
                    //                    this->Register(a);
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
                    f = static_cast<T> (3.1459) / b;
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
                    f = atl::exp((a / b));
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
                    a = .3434;
                    b = .34230;
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
                    f = atl::log((a + b));
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
                    a = .3434;
                    b = .34230;
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
                var c;
                
                PowAutoDiffTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }
                
                void Initialize() {
                    a = .00;
                    b = .00;
                    c = .000;
                    this->Register(a);
                    this->Register(b);
                    this->Register(c);
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << "\n";
                    out << "Variable b = " << b << "\n";
                    out << "f = atl::pow((a * b),(a * b))" << std::endl;
                    
                }
                
                void ObjectiveFunction(var& f) {
                    f = 1.0 - atl::pow(atl::exp(a) / atl::exp(b), atl::exp(c));
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
                    f = atl::pow(a, static_cast<T> (2.0));
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
                    f = atl::pow(static_cast<T> (2.0), (a * b));
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
                    a = 1.03434;
                    b = 2.034230;
                    this->Register(a);
                    this->Register(b);
                    //                int nobs = 50;
                    
                    T small = 1.00001212;
                    for (int i = 0; i < nobs; i++) {
                        small += .01;
                        v.push_back(var(small));
                        //
                    }
                    
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
                    out << "f.assign(Variable::gradient_structure_g, ff);\n";
                    
                }
                
                void ObjectiveFunction(var& f) {
                    //                    var ff;
                    //                    std::cout << "B id:" << b.info->id << "\n";
                    //                    var aa = a+b;
                    f = atl::sin(a * b); //atl::tanh((a / b)+(v[0]*v[1])+(v[0]*v[0])+(v[1]*v[1]));// * (v[0] * v[0]); // (a * b+a/b+a*a) - (a*(a * b+a/b+a*a));
                    //                    f = a*b*v[0];
                    //                    std::cout<<"d(a,a,a) = "<<( a*b*v[0]).EvaluateDerivative(b.info->id,b.info->id,b.info->id)<<"\n\n";
                    for (int i = 0; i < nobs; i++) {
                        for (int j = 0; j < nobs; j++) {
                            //                            int i = (nobs-1)-j;
                            f += (v[j] * v[j]) * (atl::log(a / b * v[j] * v[i]) * atl::log(a / b) * v[j] * v[i] * atl::log(a / b));
                        }
                    }
                    //                    ff*=ff;
                    //                    f = (nobs / 2.0) * (f);
                    //                    f.Assign(var::gradient_structure_g, f);
                }
                
                
            };
            
            template<class T>
            class ThirdTest : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                var a;
                var b;
                var c;
                var d;
                var g;
                
                ThirdTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }
                
                void Initialize() {
                    a = 1.0; //.03434;
                    b = 2.0; //.034230;
                    c = 3.0; //1.0213213;
                    d = 4.0; //3.0123123;
                    //                    g = 5.0; //3.0123123;
                    this->Register(a);
                    this->Register(b);
                    this->Register(c);
                    this->Register(d);
                    //                    this->Register(g);
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << ";\n";
                    out << "Variable b = " << b << ";\n";
                    //                    out << "Variable v[" << nobs << "]\n\n";
                    out << "Variable ff = 0;\n";
                    out << "ff = (a / b); \n"; //// (a * b+a/b+a*a) - (a*(a * b+a/b+a*a));
                    out << "for (int i = 0; i < nobs; i++) {\n";
                    out << "    for (int j = 0; j < nobs; j++) {\n";
                    out << "         ff += atl::log(v[j] * v[i] * atl::tanh(a / b * v[j] * v[i]) * atl::tanh(a / b) * v[j] * v[i] * atl::tanh(a / b));\n";
                    out << "     }\n";
                    out << "}\n";
                    out << "f.assign(Variable::gradient_structure_g, ff);\n";
                    std::cout << "a " << a.info->id << "\n";
                    std::cout << "b " << b.info->id << "\n";
                    std::cout << "c " << c.info->id << "\n";
                    std::cout << "d " << d.info->id << "\n";
                }
                
                void ObjectiveFunction(var& f) {
                    //                    var ff;
                    f = (a / b); //+c * d; //+c+d; //+c*d;
                    //                    f += c * d;// * atl::sin(g);
                    //                    f += c * d;// * atl::sin(g);
                    f *= c / d;
                    //                    f = atl::log((a / b)+(c*d));
                    //                    f = 2.0 * f*c*d;
                    //                    std::cout << "d3f = " << (2.0 * atl::pow(f,2.0)).EvaluateDerivative(f.info->id,f.info->id,f.info->id) << "\n";
                }
                
                
            };
            
            template<class T>
            class Ackley : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                typedef atl::Variable<T> variable;
                std::vector<variable > x;
                
                variable a;
                variable b;
                variable c;
                variable d;
                
                atl::Variable<T> f;
                std::mt19937 generator; //// mt19937 is a standard mersenne_twister_engine
                std::normal_distribution<T> distribution;
                
                Ackley(std::ofstream& out, int dd) {
                    d = (T) dd;
                    this->RunTestToFile(out);
                }
                
                void Initialize() {
                    int dimensions = (int) d.GetValue();
                    a = 20;
                    b = .2;
                    c = 2.0 * M_PI;
                    
                    //random seed
                    generator.seed(4090091);
                    //set some random starting values
                    
                    for (int i = 0; i < dimensions; i++) {
                        T r;
                        r = distribution(generator);
                        this->x.push_back(variable(r));
                    }
                    for (int i = 0; i < x.size(); i++) {
                        //                        std::cout<<x[i]<<" ";
                        this->Register(x[i]);
                    }
                    //                    std::cout<<"\n";
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable a = " << a << ";\n";
                    out << "Variable b = " << b << ";\n";
                    out << "Variable c = " << c << ";\n";
                    out << "Variable x[" << (int) d.GetValue() << "];\n";
                    out << "f = 0;" << "\n";
                    out << "variable sum1;" << ";\n";
                    out << "variable sum2;" << ";\n";
                    out << "variable term1;" << ";\n";
                    out << "variable term2;" << ";\n\n";
                    
                    out << "for (int i = 0; i < x.size(); i++) {" << "\n";
                    out << "    sum1 += atl::pow(x[i], 2.0);" << "\n";
                    out << "    sum2 += atl::cos(x[i] * c);" << "\n";
                    out << "}\n\n";
                    
                    out << "term1 = -a * atl::exp(-b * atl::sqrt(sum1 / d));" << "\n";
                    out << "term2 = -1.0 * atl::exp(sum2 / d);" << "\n\n";
                    
                    out << "f = term1 + term2 + a + exp(1);" << "\n";
                }
                
                void ObjectiveFunction(var& f) {
                    f = 0;
                    variable sum1;
                    variable sum2;
                    variable term1;
                    variable term2;
                    
                    for (int i = 0; i < x.size(); i++) {
                        sum1 += atl::pow(x[i], static_cast<T> (2.0));
                        sum2 += atl::cos(x[i] * c);
                    }
                    
                    term1 = -a * atl::exp(-b * atl::sqrt(sum1 / d));
                    term2 = static_cast<T> (-1.0) * atl::exp(sum2 / d);
                    
                    f = term1 + term2 + a + std::exp(1);
                }
                
                
            };
            
            
            //            template<typename T>
            //        const T Bertalanffy(const T& Linf, const T& a0, const T& a, const T& k) {
            //            return Linf * (1.0 - std::exp(-k * (a - a0)));
            //        }
            
            template<class T>
            class LogTheta : public AutoDiffTest<T> {
            public:
                
                std::vector<T> Y = {3.02315,
                    3.21984, 3.65424, 3.60466, 3.5423, 4.43166, 4.29131, 4.36787, 4.37504, 3.94537, 4.2912,
                    4.59584, 4.30217, 4.15193, 4.42821, 4.30537, 4.54861, 4.63205, 4.78491, 4.79962, 4.92626,
                    5.13627, 5.14859, 4.64859, 5.29066, 4.67895, 5.08076, 5.22367, 5.26112, 4.76097, 4.7894,
                    5.22118, 5.20706, 5.20803, 5.18209, 5.41249, 5.8154, 5.20875, 5.78776, 6.08957, 5.90054,
                    5.94058, 6.21408, 5.84143, 6.42039, 6.46694, 6.73742, 6.7644, 6.80683, 6.33115, 6.7795,
                    6.72164, 6.74198, 7.00933, 6.65849, 6.51895, 6.29572, 6.34795, 6.68318, 6.91444, 6.27536,
                    6.41721, 6.82602, 6.58711, 6.50546, 6.85916, 6.67842, 6.83515, 6.88816, 6.82988, 6.62378,
                    6.8927, 6.60941, 6.80232, 6.88835, 6.92828, 7.12176, 6.66136, 6.7966, 6.42991, 6.59199,
                    6.77583, 6.4914, 6.52556, 6.95658, 6.57174, 6.62285, 6.40535, 6.5243, 6.64547, 6.14812,
                    6.58134, 6.24897, 6.04781, 6.51518, 6.21539, 6.99472, 6.40441, 6.99795, 6.82203, 6.73572,
                    6.74861, 6.30817, 6.53584, 6.37093, 6.47595, 6.29266, 6.35378, 6.53895, 6.24377, 6.49304,
                    6.30005, 6.56823, 6.52514, 6.75102, 6.27068, 6.83394, 6.48833, 6.2357, 6.71897, 6.6146,
                    6.71341, 6.62067, 6.50148, 6.38419, 6.72777, 6.76925, 6.66819, 6.31156, 7.05917, 6.69037,
                    6.09881, 6.8032, 6.30866, 6.56355, 6.39023, 6.59533, 6.69288, 6.68218, 6.91985, 6.90063,
                    6.51458, 6.88214, 7.04983, 7.07688, 7.12766, 6.59529, 6.71444, 6.95628, 6.49304, 6.76456,
                    6.4496, 6.52297, 6.47359, 6.72056, 6.67134, 6.1982, 6.47844, 6.91845, 6.7916, 6.96461,
                    6.7269, 6.49073, 6.9752, 7.01309, 7.25482, 6.79274, 7.05099, 6.95183, 6.98215, 7.10636,
                    6.98072, 7.32764, 7.01756, 6.51239, 6.98783, 6.64356, 6.67166, 7.2179, 6.73284, 6.50916,
                    6.32713, 6.90174, 6.5849, 6.71981, 6.44168, 6.74218, 6.6416, 6.87571, 6.727, 7.13012,
                    7.01468, 6.92291, 6.68062, 7.07283, 6.64972, 7.01744, 6.98964, 6.71007, 6.83583};
                
                std::vector<atl::Variable<T> > X;
                atl::Variable<T> logr0; // = -2.6032947;
                atl::Variable<T> logtheta; // = 0.7625692;
                atl::Variable<T> logK; // = 6.0; //  = 6.7250075;
                atl::Variable<T> logQ; // = -4.7496015;
                atl::Variable<T> logR; // = -3.1889239;
                
                LogTheta(std::ofstream& out) {
                    this->RunTestToFile(out);
                }
                
                void Description(std::stringstream& out) {
                    out << "LogTheta model....\n";
                }
                
                void Initialize() {
                    
                    
                    
                    int size = Y.size();
                    std::cout << "size = " << size;
                    this->X.resize(size);
                    
                    
                    
                    this->Register(logr0);
                    this->Register(logtheta);
                    //                    logK.SetBounds(4.6, 7.6);
                    this->Register(logK);
                    this->Register(logQ);
                    this->Register(logR);
                    
                    for (int i = 0; i < X.size(); i++) {
                        //            X[i] = atl::Variable<T>(.012);
                        //                        this->RegisterRandomVariable(X[i]);
                    }
                }
                
                void step(atl::Variable<T>& jnll, const atl::Variable<T>& x1, const atl::Variable<T>& x2, const atl::Variable<T>& logr0, const atl::Variable<T>& logK, const atl::Variable<T>& logtheta, const atl::Variable<T>& logQ) {
                    atl::Variable<T> var = atl::exp(logQ);
                    atl::Variable<T> m = (x1 + atl::exp(logr0) * (1.0 - atl::pow(atl::exp(x1) / atl::exp(logK), atl::exp(logtheta))));
                    jnll += 0.5 * (atl::log(2.0 * M_PI * atl::exp(logQ))+((x2 - m)*(x2 - m)) / var);
                }
                
                void obs(atl::Variable<T>& jnll, const atl::Variable<T>& x, const atl::Variable<T>& logR, int i) {
                    atl::Variable<T> var = atl::exp(logR);
                    jnll += 0.5 * (atl::log(2.0 * M_PI * var)+((x - Y[i])*(x - Y[i])) / var);
                }
                
                virtual void ObjectiveFunction(atl::Variable<T>& f) {
                    
                    f = 0.0;
                    
                    
                    int timeSteps = Y.size();
                    
                    for (int i = 1; i < timeSteps; i++) {
                        step(f, X[i - 1], X[i], logr0, logK, logtheta, logQ);
                    }
                    
                    for (int i = 0; i < timeSteps; i++) {
                        obs(f, X[i], logR, i);
                    }
                    
                    
                }
                
            };
            
            template<class T>
            class LogTheta2 : public AutoDiffTest<T> {
            public:
                
                std::vector<T> Y = {3.02315,
                    3.21984, 3.65424, 3.60466, 3.5423, 4.43166, 4.29131, 4.36787, 4.37504, 3.94537, 4.2912,
                    4.59584, 4.30217, 4.15193, 4.42821, 4.30537, 4.54861, 4.63205, 4.78491, 4.79962, 4.92626,
                    5.13627, 5.14859, 4.64859, 5.29066, 4.67895, 5.08076, 5.22367, 5.26112, 4.76097, 4.7894,
                    5.22118, 5.20706, 5.20803, 5.18209, 5.41249, 5.8154, 5.20875, 5.78776, 6.08957, 5.90054,
                    5.94058, 6.21408, 5.84143, 6.42039, 6.46694, 6.73742, 6.7644, 6.80683, 6.33115, 6.7795,
                    6.72164, 6.74198, 7.00933, 6.65849, 6.51895, 6.29572, 6.34795, 6.68318, 6.91444, 6.27536,
                    6.41721, 6.82602, 6.58711, 6.50546, 6.85916, 6.67842, 6.83515, 6.88816, 6.82988, 6.62378,
                    6.8927, 6.60941, 6.80232, 6.88835, 6.92828, 7.12176, 6.66136, 6.7966, 6.42991, 6.59199,
                    6.77583, 6.4914, 6.52556, 6.95658, 6.57174, 6.62285, 6.40535, 6.5243, 6.64547, 6.14812,
                    6.58134, 6.24897, 6.04781, 6.51518, 6.21539, 6.99472, 6.40441, 6.99795, 6.82203, 6.73572,
                    6.74861, 6.30817, 6.53584, 6.37093, 6.47595, 6.29266, 6.35378, 6.53895, 6.24377, 6.49304,
                    6.30005, 6.56823, 6.52514, 6.75102, 6.27068, 6.83394, 6.48833, 6.2357, 6.71897, 6.6146,
                    6.71341, 6.62067, 6.50148, 6.38419, 6.72777, 6.76925, 6.66819, 6.31156, 7.05917, 6.69037,
                    6.09881, 6.8032, 6.30866, 6.56355, 6.39023, 6.59533, 6.69288, 6.68218, 6.91985, 6.90063,
                    6.51458, 6.88214, 7.04983, 7.07688, 7.12766, 6.59529, 6.71444, 6.95628, 6.49304, 6.76456,
                    6.4496, 6.52297, 6.47359, 6.72056, 6.67134, 6.1982, 6.47844, 6.91845, 6.7916, 6.96461,
                    6.7269, 6.49073, 6.9752, 7.01309, 7.25482, 6.79274, 7.05099, 6.95183, 6.98215, 7.10636,
                    6.98072, 7.32764, 7.01756, 6.51239, 6.98783, 6.64356, 6.67166, 7.2179, 6.73284, 6.50916,
                    6.32713, 6.90174, 6.5849, 6.71981, 6.44168, 6.74218, 6.6416, 6.87571, 6.727, 7.13012,
                    7.01468, 6.92291, 6.68062, 7.07283, 6.64972, 7.01744, 6.98964, 6.71007, 6.83583};
                
                std::vector<atl::Variable<T> > X;
                atl::Variable<T> logr0 = -2.6032947;
                atl::Variable<T> logtheta = 0.7625692;
                atl::Variable<T> logK = 6.0; //  = 6.7250075;
                atl::Variable<T> logQ = -4.7496015;
                atl::Variable<T> logR = -3.1889239;
                
                LogTheta2(std::ofstream& out) {
                    this->RunTestToFile(out);
                }
                
                void Description(std::stringstream& out) {
                    out << "LogTheta model....\n";
                }
                
                void Initialize() {
                    
                    
                    
                    int size = Y.size();
                    std::cout << "size = " << size;
                    this->X.resize(size);
                    
                    
                    
                    this->Register(logr0);
                    this->Register(logtheta);
                    //                    logK.SetBounds(4.6, 7.6);
                    this->Register(logK);
                    this->Register(logQ);
                    this->Register(logR);
                    
                    for (int i = 0; i < X.size(); i++) {
                        X[i] = atl::Variable<T>(3.012);
                        //                        this->RegisterRandomVariable(X[i]);
                    }
                }
                
                const atl::Variable<T> dnorm(const atl::Variable<T> x,
                                             const atl::Variable<T> mean,
                                             const atl::Variable<T> sd, int give_log = 0) {
                    atl::Variable<T> logres;
                    logres = -1.0 * atl::Variable<T>(atl::Variable<T>(std::sqrt(2 * M_PI)) * sd) - atl::Variable<T>(.5) * atl::pow((x - mean) / sd, 2.0);
                    if (give_log)return logres;
                    else return atl::exp(logres);
                }
                
                virtual void ObjectiveFunction(atl::Variable<T>& f) {
                    
                    size_t timeSteps = Y.size();
                    atl::Variable<T> r0 = atl::exp(logr0);
                    atl::Variable<T> theta = atl::exp(logtheta);
                    atl::Variable<T> K = atl::exp(logK);
                    atl::Variable<T> Q = atl::exp(logQ);
                    atl::Variable<T> R = atl::exp(logR);
                    
                    f = 0.0;
                    
                    for (int i = 1; i < timeSteps; i++) {
                        atl::Variable<T> m = X[i - 1] + r0 * (1.0 - atl::pow(atl::exp(X[i - 1]) / K, theta));
                        //                        std::cout << "f = " << f << "\n";
                        
                        f -= dnorm(X[i], m, atl::pow(Q, .5), true);
                    }
                    for (int i = 0; i < timeSteps; i++) {
                        f -= dnorm(Y[i], X[i], atl::pow(R, .5), true);
                    }
                    //                    return ans;
                }
                
            };
            
            template<class T>
            class DNormTest : public AutoDiffTest<T> {
                atl::Variable<T> x;
                atl::Variable<T> sd;
                atl::Variable<T> mean;
                
                const atl::Variable<T> dnorm(const atl::Variable<T> x,
                                             const atl::Variable<T> mean,
                                             const atl::Variable<T> sd, int give_log = 0) {
                    atl::Variable<T> logres;
                    logres = atl::Variable<T>(-1.0) *
                    std::sqrt(2 * M_PI) * sd -
                    atl::Variable<T>(.5) * (((x - mean) / sd)*((x - mean) / sd));
                    if (give_log)return logres;
                    else return atl::exp(logres);
                }
                
                
            public:
                
                DNormTest(std::ofstream& out) {
                    this->RunTestToFile(out);
                }
                
                void Description(std::stringstream& out) {
                    out << "dnorm model....\n";
                }
                
                void Initialize() {
                    
                    x = .5;
                    mean = .123123;
                    sd = .25;
                    
                    this->Register(x);
                    this->Register(mean);
                    this->Register(sd);
                }
                
                
                virtual void ObjectiveFunction(atl::Variable<T>& f) {
                    f = dnorm(x,mean,sd,false);
                }
                
                
            };
            
            template<class T>
            class Bertalanffy : public AutoDiffTest<T> {
            public:
                typedef atl::Variable<T> var;
                typedef atl::Variable<T> variable;
                std::vector<T > y;
                
                std::vector<T> a = {2, 2, 2, 2, 3, 3, 3, 4, 4, 6, 9, 10, 11, 11, 15, 19, 21, 24, 30, 35};
                std::vector<T> l = {1, 3, 4, 4, 3, 4, 5, 6, 9, 10, 14, 13, 16, 17, 18, 19, 20, 21, 20, 21};
                variable Linf;
                variable a0;
                variable k;
                variable sd;
                
                atl::Variable<T> f;
                std::mt19937 generator; //// mt19937 is a standard mersenne_twister_engine
                std::normal_distribution<T> distribution;
                
                Bertalanffy(std::ofstream& out) {
                    this->RunTestToFile(out);
                }
                
                void Initialize() {
                    Linf = 22.1724234324;
                    this->Register(Linf);
                    
                    a0 = 0.929192342346;
                    this->Register(a0);
                    
                    k = 0.113242353188;
                    this->Register(k);
                    
                    sd = 0.289336;
                    this->Register(sd);
                    
                    
                    
                    
                }
                
                void Description(std::stringstream& out) {
                    out << "Test Problem:\n";
                    out << "Parameters:\n";
                    out << "Variable Linf = " << Linf << ";\n";
                    out << "Variable a0= " << a0 << ";\n";
                    out << "Variable k = " << k << ";\n";
                    out << "Variable sd = " << sd << ";\n";
                    out << "vector a[" << a.size() << "];\n";
                    out << "vector l[" << a.size() << "];\n";
                    out << "f = 0;" << "\n";
                    out << "variable sum;" << ";\n";
                    
                    out << "for (int i = 0; i < y.size(); i++) {" << "\n";
                    out << "    variable predL = (Linf * (static_cast<T> (1.0) - atl::exp(static_cast<T> (-1.0) * k * (a[i] - a0))));\n";
                    out << "    sum += static_cast<T> (.5) * atl::square((std::log(l[i]) - atl::log(predL))) / sd;\n";
                    out << "}\n\n";
                    
                    
                    out << "T(a.size()) * atl::log(sd) * sum" << "\n";
                }
                
                void ObjectiveFunction(var& f) {
                    f = 0;
                    variable sum;
                    for (int i = 0; i < a.size(); i++) {
                        variable predL = (Linf * (static_cast<T> (1.0) - atl::exp(static_cast<T> (-1.0) * k * (a[i] - a0))));
                        sum += static_cast<T> (.5) * ((std::log(l[i]) - atl::log(predL))*(std::log(l[i]) - atl::log(predL))) / sd;
                    }
                    
                    
                    f = T(a.size()) * atl::log(sd) * sum;
                }
                
                
            };
            
            template<typename T>
            struct TestTypeTrait {
                static std::string type;
            };
            template<typename T>
            std::string TestTypeTrait<T>::type = "unknown precision";
            
            template<>
            std::string TestTypeTrait<double>::type = "double";
            
            template<>
            std::string TestTypeTrait<long double>::type = "long double";
            

            
            void Run() {
                typedef double real_t;
                std::ofstream out;
                out.open("autodiff_tests.txt");
                out << "This test compares the computed gradient, second-order, "
                "and third-order mixed \npartial derivatives against "
                "those computed using the differences method.\n\n";
                
                out << "Data type: " << TestTypeTrait<real_t>::type << "\n\n";
                out << "Test Gradient Tolerance: " << atl::tests::auto_diff::AutoDiffTest<real_t>::GetTolerance() << "\n\n";
                out << "Test Second-Order Tolerance: " << atl::tests::auto_diff::AutoDiffTest<real_t>::GetSecondOrderTolerance() << "\n\n";
                out << "Test Third-Order Tolerance: " << atl::tests::auto_diff::AutoDiffTest<real_t>::GetThirdOrderTolerance() << "\n\n";
                
                std::cout << "running...\n";
                atl::tests::auto_diff::AddAutoDiffTest<real_t> add(out);
                atl::tests::auto_diff::Add1AutoDiffTest<real_t> add1(out);
                atl::tests::auto_diff::Add2AutoDiffTest<real_t> add2(out);
                atl::tests::auto_diff::SubtractAutoDiffTest<real_t> subtract1(out);
                atl::tests::auto_diff::Subtract1AutoDiffTest<real_t> subtract(out);
                atl::tests::auto_diff::Subtract2AutoDiffTest<real_t> subtract2(out);
                atl::tests::auto_diff::MultiplyAutoDiffTest<real_t> multiply(out);
                atl::tests::auto_diff::Multiply1AutoDiffTest<real_t> multiply1(out);
                atl::tests::auto_diff::Multiply2AutoDiffTest<real_t> multiply2(out);
                atl::tests::auto_diff::DivideAutoDiffTest<real_t> divide(out);
                atl::tests::auto_diff::Divide1AutoDiffTest<real_t> divide1(out);
                atl::tests::auto_diff::Divide2AutoDiffTest<real_t> divide2(out);
                atl::tests::auto_diff::CosAutoDiffTest<real_t> cos(out);
                atl::tests::auto_diff::ACosAutoDiffTest<real_t> acos(out);
                atl::tests::auto_diff::SinAutoDiffTest<real_t> sin(out);
                atl::tests::auto_diff::ASinAutoDiffTest<real_t> asin(out);
                atl::tests::auto_diff::TanAutoDiffTest<real_t> tan(out);
                atl::tests::auto_diff::ATanAutoDiffTest<real_t> atan(out);
                atl::tests::auto_diff::CoshAutoDiffTest<real_t> cosh(out);
                atl::tests::auto_diff::SinhAutoDiffTest<real_t> sinh(out);
                atl::tests::auto_diff::TanhAutoDiffTest<real_t> tanh(out);
                atl::tests::auto_diff::ExpAutoDiffTest<real_t> exp(out);
                atl::tests::auto_diff::LogAutoDiffTest<real_t> log(out);
                atl::tests::auto_diff::Log10AutoDiffTest<real_t> log10(out);
                atl::tests::auto_diff::FabsAutoDiffTest<real_t> fabs(out);
                atl::tests::auto_diff::SqrtAutoDiffTest<real_t> sqrt(out);
                atl::tests::auto_diff::PowAutoDiffTest<real_t> pow(out);
                atl::tests::auto_diff::PowCAutoDiffTest<real_t> pow2(out);
                atl::tests::auto_diff::PowC2AutoDiffTest<real_t> pow3(out);
                atl::tests::auto_diff::CeilAutoDiffTest<real_t> ceil(out);
                atl::tests::auto_diff::FloorAutoDiffTest<real_t> floor(out);
                
#ifdef  HESSIAN_USE_AD_GRADIENT
                //                atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<real_t> s3(out, 50);
#endif
                atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<real_t> s1(out, 10);
                atl::tests::auto_diff::SumAlotOfParametersAutoDiffTest<real_t> s2(out, 20);
                
                atl::tests::auto_diff::Ackley<real_t> ackley10(out, 10);
                atl::tests::auto_diff::Ackley<real_t> ackley20(out, 20);
                atl::tests::auto_diff::Ackley<real_t> ackley30(out, 30);
                atl::tests::auto_diff::Ackley<real_t> ackley40(out, 40);
                atl::tests::auto_diff::Ackley<real_t> ackley50(out, 50);
#ifdef  HESSIAN_USE_AD_GRADIENT
                atl::tests::auto_diff::Ackley<real_t> ackley100(out, 100);
                atl::tests::auto_diff::Ackley<real_t> ackley200(out, 200);
#endif
                //
                atl::tests::auto_diff::Bertalanffy<real_t> bertalanffy(out);
                atl::tests::auto_diff::LogTheta<real_t> lt(out);
                atl::tests::auto_diff::DNormTest<real_t> dnt(out);
                std::cout << "Test complete.\n";
                if (atl::tests::auto_diff::fail == 0) {
                    std::cout << "All tests passed, review file \"autodiff_test.txt\" for details." << std::endl;
                    
                } else {
                    
                    std::cout << atl::tests::auto_diff::fail << " of " << atl::tests::auto_diff::tests << " tests did not agree within a tolerance of " << atl::tests::auto_diff::AutoDiffTest<real_t>::GetTolerance() << ", review file \"autodiff_test.txt\" for details." << std::endl;
                }
                
            }

            
        }
    }
}




#endif	/* HESSIANTESTS_HPP */