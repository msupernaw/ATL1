
/* 
 * File:   DynamicExpression.hpp
 * Author: matthewsupernaw
 *
 * Created on December 29, 2015, 10:06 PM
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


#ifndef DYNAMICEXPRESSION_HPP
#define DYNAMICEXPRESSION_HPP

#include "GradientStructure.hpp"
#include "Expression.hpp"
#include "VariableInfo.hpp"
#include "GradientStructure.hpp"
#include "../Utilities/BigFloat.hpp"
#include "third_party/dlmalloc/malloc.h"
#include <cmath>
#include <sstream>

#define DYNAMIC_AD_LOG10 2.30258509299404590109361379290930926799774169921875

#define MAKE_UNARY_DYNAMIC_EXPRESSION(NAME, EVALUATE, FIRST_DERIVATIVE_EVAL, SECOND_DERIVATIVE_EVAL, DYNAMIC_DIFFERENTIATE ) \
\
template<class REAL_T> \
class NAME : public atl::DynamicExpression< REAL_T> {\
atl::DynamicExpression< REAL_T> * expr_m;    \
public:  \
\
NAME(atl::DynamicExpression< REAL_T>* expr):expr_m(expr) {}\
\
\
        virtual inline const REAL_T Evaluate() { \
        EVALUATE \
         }\
        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt){ \
        FIRST_DERIVATIVE_EVAL \
        } \
       virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) { \
       SECOND_DERIVATIVE_EVAL \
        } \
       virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) { \
                expr_m->push_ids(ids);    \
        } \
       virtual DynamicExpression<REAL_T>* Differentiate() { \
       DYNAMIC_DIFFERENTIATE \
       }\
       };\

namespace atl {
    template<typename REAL_T>
    class DynamicLog;

    template<typename REAL_T>
    class DynamicSqrt;

    template<typename REAL_T>
    class DynamicSin;

    template<typename REAL_T>
    class DynamicSinh;

    /**
     * Base class for dynamic expression types.
     * 
     * Used to record expressions on a tape. Dynamic Expressions are 
     * dynamically allocated from the heap.
     */
    template<typename REAL_T>
    class DynamicExpression {
    public:

        DynamicExpression() {

        }

        virtual ~DynamicExpression() {

        }

        virtual inline const REAL_T Evaluate() = 0;
        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) = 0;
        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) = 0;
        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) = 0;
        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t id) = 0;
        virtual DynamicExpression<REAL_T>* Differentiate() = 0;
        virtual DynamicExpression<REAL_T>* Clone() = 0;
        virtual std::string ToString() = 0;

    };

    //    MAKE_UNARY_DYNAMIC_EXPRESSION(DACOS,
    //    return std::acos(expr_m->Evaluate());,
    //            REAL_T fx = expr_m->Evaluate();
    //    return -1.0 * expr_m->EvaluateDerivative(wrt) / std::sqrt(1.0 - fx * fx);,
    //            REAL_T fx = expr_m->Evaluate();
    //    return (((-1.0 * (fx * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)))
    //            / std::pow((1.0 - fx * fx), 1.5))
    //            - (expr_m->EvaluateDerivative(wrt_x, wrt_y) / std::sqrt(1.0 - fx * fx)));,
    //    return NULL;
    //            )

    template<typename REAL_T>
    class DynamicScalar : public atl::DynamicExpression<REAL_T> {
        REAL_T value_m;
    public:

        DynamicScalar(REAL_T value = static_cast<REAL_T> (0.0)) : value_m(value) {

        }

        ~DynamicScalar() {

        }

        virtual inline const REAL_T Evaluate() {
            return value_m;
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return static_cast<REAL_T> (0.0);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return static_cast<REAL_T> (0.0);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t id) {
            return new DynamicScalar<REAL_T>();
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicScalar<REAL_T>();
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicScalar<REAL_T>(this->value_m);
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << this->value_m;
            return ss.str();
        }
    };

    template<typename REAL_T>
    class DynamicVariable : public atl::DynamicExpression<REAL_T> {
        atl::VariableInfo<REAL_T>* info_m;
    public:

        DynamicVariable(atl::VariableInfo<REAL_T>* info) : info_m(info) {

        }

        ~DynamicVariable() {

        }

        void* operator new(size_t size) {
            return malloc(size);
        }

        void operator delete(void* ptr) {
            free(ptr);
        }

        virtual inline const REAL_T Evaluate() {
            return info_m->vvalue;
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return info_m->id == wrt ? 1 : 0;
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return static_cast<REAL_T> (0.0);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            ids.insert(this->info_m);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t id) {
            return info_m->id == id ? new DynamicScalar<REAL_T>(1) : new DynamicScalar<REAL_T>(0);
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicScalar<REAL_T>(1.0);
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicVariable<REAL_T>(this->info_m);
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "variable(" << this->info_m->id << "," << this->info_m->vvalue << ")";
            return ss.str();
        }
    };

    template<typename REAL_T>
    class DynamicAdd : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* lhs_m;
        atl::DynamicExpression<REAL_T>* rhs_m;
    public:

        DynamicAdd(atl::DynamicExpression<REAL_T>* lhs_m, atl::DynamicExpression<REAL_T>* rhs_m) :
        DynamicExpression<REAL_T>(), lhs_m(lhs_m), rhs_m(rhs_m) {
        }

        ~DynamicAdd() {
            delete lhs_m;
            delete rhs_m;
        }

        void* operator new(size_t size) {
            return malloc(size);
        }

        void operator delete(void* ptr) {
            free(ptr);
        }

        virtual inline const REAL_T Evaluate() {
            return lhs_m->Evaluate() + rhs_m->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return lhs_m->EvaluateDerivative(wrt) + rhs_m->EvaluateDerivative(wrt);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return lhs_m->EvaluateDerivative(wrt_x, wrt_y) + rhs_m->EvaluateDerivative(wrt_x, wrt_y);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            lhs_m->PushIds(ids);
            rhs_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t id) {
            return new DynamicAdd<REAL_T>(lhs_m->Differentiate(id), rhs_m->Differentiate(id));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicAdd<REAL_T>(lhs_m->Differentiate(), rhs_m->Differentiate());
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicAdd<REAL_T>(lhs_m->Clone(), rhs_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "(" << this->lhs_m->ToString() << " + " << this->rhs_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicSubtract : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* lhs_m;
        atl::DynamicExpression<REAL_T>* rhs_m;
    public:

        DynamicSubtract(atl::DynamicExpression<REAL_T>* lhs_m, atl::DynamicExpression<REAL_T>* rhs_m) :
        DynamicExpression<REAL_T>(), lhs_m(lhs_m), rhs_m(rhs_m) {
        }

        ~DynamicSubtract() {
            delete lhs_m;
            delete rhs_m;
        }

        virtual inline const REAL_T Evaluate() {
            return lhs_m->Evaluate() - rhs_m->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return lhs_m->EvaluateDerivative(wrt) - rhs_m->EvaluateDerivative(wrt);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return lhs_m->EvaluateDerivative(wrt_x, wrt_y) - rhs_m->EvaluateDerivative(wrt_x, wrt_y);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            lhs_m->PushIds(ids);
            rhs_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t id) {
            return new DynamicSubtract<REAL_T>(lhs_m->Differentiate(id), rhs_m->Differentiate(id));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicSubtract<REAL_T>(lhs_m->Differentiate(), rhs_m->Differentiate());
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicSubtract<REAL_T>(lhs_m->Clone(), rhs_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "(" << this->lhs_m->ToString() << " - " << this->rhs_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicMultiply : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* lhs_m;
        atl::DynamicExpression<REAL_T>* rhs_m;
    public:

        DynamicMultiply(atl::DynamicExpression<REAL_T>* lhs_m, atl::DynamicExpression<REAL_T>* rhs_m) :
        DynamicExpression<REAL_T>(), lhs_m(lhs_m), rhs_m(rhs_m) {
        }

        ~DynamicMultiply() {
            delete lhs_m;
            delete rhs_m;
        }

        virtual inline const REAL_T Evaluate() {
            return lhs_m->Evaluate() * rhs_m->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return (lhs_m->Evaluate() * rhs_m->EvaluateDerivative(wrt) +
                    lhs_m->EvaluateDerivative(wrt) * rhs_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return lhs_m->EvaluateDerivative(wrt_x) * rhs_m->EvaluateDerivative(wrt_y) + lhs_m->Evaluate() * rhs_m->EvaluateDerivative(wrt_x, wrt_y) +
                    lhs_m->EvaluateDerivative(wrt_y) * rhs_m->EvaluateDerivative(wrt_x) + rhs_m->Evaluate() * lhs_m->EvaluateDerivative(wrt_x, wrt_y);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            lhs_m->PushIds(ids);
            rhs_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t id) {
            return new DynamicAdd<REAL_T>(
                    new DynamicMultiply<REAL_T>(lhs_m->Clone(), rhs_m->Differentiate(id)),
                    new DynamicMultiply<REAL_T>(lhs_m->Differentiate(id), rhs_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicAdd<REAL_T>(
                    new DynamicMultiply<REAL_T>(lhs_m->Clone(), rhs_m->Differentiate()),
                    new DynamicMultiply<REAL_T>(lhs_m->Differentiate(), rhs_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicMultiply(lhs_m->Clone(), rhs_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "(" << this->lhs_m->ToString() << " * " << this->rhs_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicDivide : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* lhs_m;
        atl::DynamicExpression<REAL_T>* rhs_m;
    public:

        DynamicDivide(atl::DynamicExpression<REAL_T>* lhs_m, atl::DynamicExpression<REAL_T>* rhs_m) :
        DynamicExpression<REAL_T>(), lhs_m(lhs_m), rhs_m(rhs_m) {
        }

        ~DynamicDivide() {
            delete lhs_m;
            delete rhs_m;
        }

        virtual inline const REAL_T Evaluate() {
            return lhs_m->Evaluate() / rhs_m->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return (rhs_m->Evaluate() * lhs_m->EvaluateDerivative(wrt) -
                    rhs_m->EvaluateDerivative(wrt) * lhs_m->Evaluate()) / (rhs_m->Evaluate() * rhs_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T aa = ((2.0 * lhs_m->Evaluate() * rhs_m->EvaluateDerivative(wrt_x) * rhs_m->EvaluateDerivative(wrt_y)) / (rhs_m->Evaluate() * rhs_m->Evaluate() * rhs_m->Evaluate())); //(2*f(a,b)*('diff(g(a,b),a,1))*('diff(g(a,b),b,1)))/g(a,b)^3
            REAL_T bb = ((lhs_m->EvaluateDerivative(wrt_x) * rhs_m->EvaluateDerivative(wrt_y)) / (rhs_m->Evaluate() * rhs_m->Evaluate())); //(('diff(f(wrt_x,b),a,1))*('diff(g(wrt_x,b),b,1)))/g(wrt_x,b)^2
            REAL_T cc = ((lhs_m->Evaluate() * rhs_m->EvaluateDerivative(wrt_x, wrt_y)) / (rhs_m->Evaluate() * rhs_m->Evaluate())); //(f(wrt_x,b)*('diff(g(wrt_x,b),a,1,b,1)))/g(wrt_x,b)^2
            REAL_T dd = ((lhs_m->EvaluateDerivative(wrt_y) * rhs_m->EvaluateDerivative(wrt_x)) / (rhs_m->Evaluate() * rhs_m->Evaluate())); //(('diff(f(wrt_x,b),b,1))*('diff(g(wrt_x,b),a,1)))/g(wrt_x,b)^2
            REAL_T ee = lhs_m->EvaluateDerivative(wrt_x, wrt_y) / rhs_m->Evaluate();
            return aa - bb - cc - dd - ee;

        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            lhs_m->PushIds(ids);
            rhs_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicDivide(
                    new DynamicSubtract<REAL_T>(
                    new DynamicMultiply<REAL_T>(rhs_m->Clone(), lhs_m->Differentiate(wrt)),
                    new DynamicMultiply<REAL_T>(rhs_m->Differentiate(wrt), lhs_m->Clone())),
                    new DynamicMultiply<REAL_T>(rhs_m->Clone(), rhs_m->Clone()));

        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicDivide(
                    new DynamicSubtract<REAL_T>(
                    new DynamicMultiply<REAL_T>(rhs_m->Clone(), lhs_m->Differentiate()),
                    new DynamicMultiply<REAL_T>(rhs_m->Differentiate(), lhs_m->Clone())),
                    new DynamicMultiply<REAL_T>(rhs_m->Clone(), rhs_m->Clone()));

        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicDivide(this->lhs_m->Clone(), this->rhs_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "(" << this->lhs_m->ToString() << " / " << this->rhs_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicPow : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* lhs_m;
        atl::DynamicExpression<REAL_T>* rhs_m;
    public:

        DynamicPow(atl::DynamicExpression<REAL_T>* lhs_m, atl::DynamicExpression<REAL_T>* rhs_m) :
        DynamicExpression<REAL_T>(), lhs_m(lhs_m), rhs_m(rhs_m) {
        }

        ~DynamicPow() {
            delete lhs_m;
            delete rhs_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::pow(lhs_m->Evaluate(), rhs_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            REAL_T g = rhs_m->Evaluate();
            REAL_T f = lhs_m->Evaluate();
            REAL_T fx = lhs_m->EvaluateDerivative(wrt);
            REAL_T gx = rhs_m->EvaluateDerivative(wrt);
            return std::pow(f, g)*(std::log(f) * gx + g * fx / f);

        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T fxy = lhs_m->EvaluateDerivative(wrt_x, wrt_y);
            REAL_T g = rhs_m->Evaluate();
            REAL_T f = lhs_m->Evaluate();
            REAL_T fx = lhs_m->EvaluateDerivative(wrt_x);
            REAL_T gy = rhs_m->EvaluateDerivative(wrt_y);
            REAL_T fy = lhs_m->EvaluateDerivative(wrt_y);
            REAL_T gx = rhs_m->EvaluateDerivative(wrt_x);
            REAL_T gxy = rhs_m->EvaluateDerivative(wrt_x, wrt_y);
            return std::pow(f, g)*(((fx * gy) / f) + std::log(f) * gxy + (fy * gx / f) -
                    (g * fx * fy) / (f * f) + g * fxy / f) + std::pow(f, g)*(std::log(f) * gx +
                    g * fx / f)*(std::log(f) * gy + g * fy / f);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            lhs_m->PushIds(ids);
            rhs_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            DynamicExpression<REAL_T>* g = rhs_m->Clone();
            DynamicExpression<REAL_T>* f = lhs_m->Clone();
            DynamicExpression<REAL_T>* fx = lhs_m->Differentiate(wrt);
            DynamicExpression<REAL_T>* gx = rhs_m->Differentiate(wrt);

            return new DynamicMultiply<REAL_T>(
                    new DynamicPow<REAL_T>(f, g),
                    new DynamicAdd<REAL_T>(
                    new DynamicMultiply<REAL_T>(new DynamicLog<REAL_T>(f), gx),
                    new DynamicDivide<REAL_T>(new DynamicMultiply<REAL_T>(g, fx), f)
                    ));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            DynamicExpression<REAL_T>* g = rhs_m->Clone();
            DynamicExpression<REAL_T>* f = lhs_m->Clone();
            DynamicExpression<REAL_T>* fx = lhs_m->Differentiate();
            DynamicExpression<REAL_T>* gx = rhs_m->Differentiate();

            return new DynamicMultiply<REAL_T>(
                    new DynamicPow<REAL_T>(f, g),
                    new DynamicAdd<REAL_T>(
                    new DynamicMultiply<REAL_T>(new DynamicLog<REAL_T>(f), gx),
                    new DynamicDivide<REAL_T>(new DynamicMultiply<REAL_T>(g, fx), f)
                    ));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicPow(lhs_m->Clone(), rhs_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "pow(" << this->lhs_m->ToString() << "," << this->rhs_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicACos : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicACos(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicACos() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::acos(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            REAL_T fx = expr_m->Evaluate();
            return -1.0 * expr_m->EvaluateDerivative(wrt) / std::sqrt(1.0 - fx * fx);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T fx = expr_m->Evaluate();
            return (((-1.0 * (fx * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)))
                    / std::pow((1.0 - fx * fx), 1.5))
                    - (expr_m->EvaluateDerivative(wrt_x, wrt_y) / std::sqrt(1.0 - fx * fx)));

        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicDivide<REAL_T>(
                    new DynamicMultiply<REAL_T>(new DynamicScalar<REAL_T>(-1.0), expr_m->Differentiate(wrt)),
                    new DynamicSqrt<REAL_T>(new DynamicSubtract<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicMultiply<REAL_T>(expr_m->Clone(), expr_m->Clone()))));

        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicDivide<REAL_T>(
                    new DynamicMultiply<REAL_T>(new DynamicScalar<REAL_T>(-1.0), expr_m->Differentiate()),
                    new DynamicSqrt<REAL_T>(new DynamicSubtract<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicMultiply<REAL_T>(expr_m->Clone(), expr_m->Clone()))));

        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicACos(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "acos(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicASin : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicASin(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicASin() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::asin(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            REAL_T fx = expr_m->Evaluate();
            return expr_m->EvaluateDerivative(wrt) / std::sqrt(1.0 - fx * fx);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T fx = expr_m->Evaluate();
            return ((((fx * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)))
                    / std::pow((1.0 - fx * fx), 1.5))
                    +(expr_m->EvaluateDerivative(wrt_x, wrt_y) / std::sqrt(1.0 - fx * fx)));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            //            DynamicExpression<REAL_T>* fx = expr_m->Clone();
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(wrt),
                    new DynamicSqrt<REAL_T>(new DynamicSubtract<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicMultiply<REAL_T>(expr_m->Clone(), expr_m->Clone()))));

        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            //            DynamicExpression<REAL_T>* fx = expr_m->Clone();
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(),
                    new DynamicSqrt<REAL_T>(new DynamicSubtract<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicMultiply<REAL_T>(expr_m->Clone(), expr_m->Clone()))));

        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicASin(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "asin(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicATan : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicATan(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicATan() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::atan(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            REAL_T fx = expr_m->Evaluate();
            return expr_m->EvaluateDerivative(wrt) / (fx * fx + 1.0);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T fx = expr_m->Evaluate();
            return (expr_m->EvaluateDerivative(wrt_x, wrt_y) / (fx * fx + 1.0)) -
                    (2.0 * fx * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)) / ((fx * fx + 1.0)*(fx * fx + 1.0));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(wrt),
                    new DynamicAdd<REAL_T>(new DynamicMultiply<REAL_T>(expr_m->Clone(), expr_m->Clone()), new DynamicScalar<REAL_T>(1.0)));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(),
                    new DynamicAdd<REAL_T>(new DynamicMultiply<REAL_T>(expr_m->Clone(), expr_m->Clone()), new DynamicScalar<REAL_T>(1.0)));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicATan<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "atan(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicCeil : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicCeil(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicCeil() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::ceil(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            exit(0);
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            exit(0);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            exit(0);
        }
        
         virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            std::cout << __FILE__ << ":" << __LINE__ << ":" << __func__ << " not yet implemented.";
            exit(0);
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicCeil<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "aceil(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicCos : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicCos(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicCos() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::cos(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return -1.0 * expr_m->EvaluateDerivative(wrt) * std::sin(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return ((-1.0 * std::cos(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y))
                    - (std::sin(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x, wrt_y)));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicMultiply<REAL_T>(new DynamicMultiply<REAL_T>(new DynamicScalar<REAL_T>(-1.0), expr_m->Differentiate(wrt)),
                    new DynamicSin<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicMultiply<REAL_T>(new DynamicMultiply<REAL_T>(new DynamicScalar<REAL_T>(-1.0), expr_m->Differentiate()),
                    new DynamicSin<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicCos<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "cos(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicCosh : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicCosh(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicCosh() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::cosh(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) * std::sinh(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return ((std::cosh(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y))
                    + (std::sinh(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x, wrt_y)));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt), new DynamicSinh<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(), new DynamicSinh<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicCosh<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "cosh(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicExp : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicExp(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicExp() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::exp(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) * this->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T fx = this->Evaluate();
            return ((fx * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y))
                    + (fx * expr_m->EvaluateDerivative(wrt_x, wrt_y)));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt), this->Clone());
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(), this->Clone());
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicExp<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "exp(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicFabs : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicFabs(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicFabs() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::fabs(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return (expr_m->EvaluateDerivative(wrt) * expr_m->Evaluate()) / this->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return (expr_m->EvaluateDerivative(wrt_x, wrt_y) * expr_m->Evaluate()) / this->Evaluate();
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicDivide<REAL_T>(new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt), expr_m->Clone()), this->Clone());
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicDivide<REAL_T>(new DynamicMultiply<REAL_T>(expr_m->Differentiate(), expr_m->Clone()), this->Clone());
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicFabs<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "fabs(" << this->expr_m->ToString() << ")";
            return ss.str();
        }
    };

    template<typename REAL_T>
    class DynamicFloor : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicFloor(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicFloor() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::floor(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return (expr_m->EvaluateDerivative(wrt) * this->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return (expr_m->EvaluateDerivative(wrt_x, wrt_y) * this->Evaluate());
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }
#warning this needs review, return clone or dynamic scalar???

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt), new DynamicScalar<REAL_T>(this->Evaluate()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(), new DynamicScalar<REAL_T>(this->Evaluate()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicFloor<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "floor(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicLog : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicLog(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicLog() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::log(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) / expr_m->Evaluate();
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return (expr_m->EvaluateDerivative(wrt_x, wrt_y) / expr_m->Evaluate()) - (expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)) / (expr_m->Evaluate() * expr_m->Evaluate());
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(wrt), expr_m->Clone());
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(), expr_m->Clone());
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicLog(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "log(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicLog10 : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicLog10(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicLog10() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::log10(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return (expr_m->EvaluateDerivative(wrt) / (DYNAMIC_AD_LOG10 * expr_m->Evaluate()));
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T fx = expr_m->Evaluate();
            return (expr_m->EvaluateDerivative(wrt_x, wrt_y) / (DYNAMIC_AD_LOG10 * fx)) -
                    ((expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)) / (DYNAMIC_AD_LOG10 * (fx * fx)));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {

            return new DynamicDivide<REAL_T>(expr_m->Differentiate(wrt),
                    new DynamicMultiply<REAL_T>(
                    new DynamicScalar<REAL_T>(DYNAMIC_AD_LOG10), expr_m->Clone()));

        }

        virtual DynamicExpression<REAL_T>* Differentiate() {

            return new DynamicDivide<REAL_T>(expr_m->Differentiate(),
                    new DynamicMultiply<REAL_T>(
                    new DynamicScalar<REAL_T>(DYNAMIC_AD_LOG10), expr_m->Clone()));

        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicLog10(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "log10(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicSin : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicSin(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicSin() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::sin(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) * std::cos(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return (std::cos(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x, wrt_y))-
                    std::sin(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt),
                    new DynamicCos<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(),
                    new DynamicCos<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicSin<REAL_T>(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "sin(" << this->expr_m->ToString() << ")";
            return ss.str();
        }
    };

    template<typename REAL_T>
    class DynamicSinh : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicSinh(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicSinh() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::sinh(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) * std::cosh(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return ((std::sinh(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y))
                    + (std::cosh(expr_m->Evaluate()) * expr_m->EvaluateDerivative(wrt_x, wrt_y)));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt),
                    new DynamicCosh<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(),
                    new DynamicCosh<REAL_T>(expr_m->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicSinh(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "sinh(" << this->expr_m->ToString() << ")";
            return ss.str();
        }
    };

    template<typename REAL_T>
    class DynamicSqrt : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicSqrt(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicSqrt() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::sqrt(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) / (2.0 * this->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            return (expr_m->EvaluateDerivative(wrt_x, wrt_y) / (2.0 * this->Evaluate())) -
                    (expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y)) / (4.0 * std::pow(expr_m->Evaluate(), 1.5));
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(wrt), new DynamicMultiply<REAL_T>(new DynamicScalar<REAL_T>(2.0), new DynamicSqrt<REAL_T>(this->expr_m->Clone())));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            return new DynamicDivide<REAL_T>(expr_m->Differentiate(), new DynamicMultiply<REAL_T>(new DynamicScalar<REAL_T>(2.0), new DynamicSqrt<REAL_T>(this->expr_m->Clone())));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicSqrt<REAL_T>(this->expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "sqrt(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicTan : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicTan(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicTan() {
            delete expr_m;
        }

        virtual inline const REAL_T Evaluate() {
            return std::tan(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            return expr_m->EvaluateDerivative(wrt) * (1.0 / std::cos(expr_m->Evaluate()))*(1.0 / std::cos(expr_m->Evaluate()));
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T sec2 = (1.0 / std::cos(expr_m->Evaluate())) * (1.0 / std::cos(expr_m->Evaluate()));
            return 2.0 * sec2 * this->Evaluate() * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y) +
                    sec2 * expr_m->EvaluateDerivative(wrt_x, wrt_y);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {
            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            DynamicExpression<REAL_T>* s = new DynamicDivide<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicCos<REAL_T>(expr_m->Clone()));
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt), new DynamicMultiply<REAL_T>(s, s->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            DynamicExpression<REAL_T>* s = new DynamicDivide<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicCos<REAL_T>(expr_m->Clone()));
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(), new DynamicMultiply<REAL_T>(s, s->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicTan(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "tan(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };

    template<typename REAL_T>
    class DynamicTanh : public atl::DynamicExpression<REAL_T> {
        atl::DynamicExpression<REAL_T>* expr_m;
    public:

        DynamicTanh(atl::DynamicExpression<REAL_T>* exp) : expr_m(exp) {

        }

        ~DynamicTanh() {
            delete expr_m;
        }

        void* operator new(size_t size) {
            return malloc(size);
        }

        void operator delete(void* ptr) {
            free(ptr);
        }

        virtual inline const REAL_T Evaluate() {
            return std::tanh(expr_m->Evaluate());
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt) {
            REAL_T sech2 = (1.0 / std::cosh(expr_m->Evaluate()))*(1.0 / std::cosh(expr_m->Evaluate()));

            return expr_m->EvaluateDerivative(wrt) * sech2;
        }

        virtual inline const REAL_T EvaluateDerivative(uint32_t wrt_x, uint32_t wrt_y) {
            REAL_T sech2 = (1.0 / std::cosh(expr_m->Evaluate()))*(1.0 / std::cosh(expr_m->Evaluate()));

            return sech2 * expr_m->EvaluateDerivative(wrt_x, wrt_y) - 2.0 * sech2 * this->Evaluate() * expr_m->EvaluateDerivative(wrt_x) * expr_m->EvaluateDerivative(wrt_y);
        }

        virtual void PushIds(IDSet<atl::VariableInfo<REAL_T>* >& ids) {

            expr_m->PushIds(ids);
        }

        virtual DynamicExpression<REAL_T>* Differentiate(uint32_t wrt) {
            DynamicExpression<REAL_T>* s = new DynamicDivide<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicCosh<REAL_T>(expr_m->Clone()));
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(wrt), new DynamicMultiply<REAL_T>(s, s->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Differentiate() {
            DynamicExpression<REAL_T>* s = new DynamicDivide<REAL_T>(new DynamicScalar<REAL_T>(1.0), new DynamicCosh<REAL_T>(expr_m->Clone()));
            return new DynamicMultiply<REAL_T>(expr_m->Differentiate(), new DynamicMultiply<REAL_T>(s, s->Clone()));
        }

        virtual DynamicExpression<REAL_T>* Clone() {
            return new DynamicTanh(expr_m->Clone());
        }

        virtual std::string ToString() {
            std::stringstream ss;
            ss << "tanh(" << this->expr_m->ToString() << ")";
            return ss.str();
        }

    };



















}


#endif /* DYNAMICEXPRESSION_HPP */

