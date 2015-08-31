/* 
 * File:   ExpressionNode.hpp
 * Author: matthewsupernaw
 *
 * Created on April 27, 2015, 9:10 AM
 */

#ifndef EXPRESSIONNODE_HPP
#define	EXPRESSIONNODE_HPP

#include <cmath>
#include "../AutoDiff/third_party/clfmalloc.h"
namespace atl {

    enum Node_Operator {
        NODE_ADD = 0,
        NODE_MINUS,
        NODE_MULTIPLY,
        NODE_DIVIDE,
        NODE_POW,
        NODE_VARIABLE,
        NODE_SCALAR,
        NODE_LOG,
        NODE_LOG10,
        NODE_EXP,
        //        NODE_MFEXP,
        NODE_ACOS,
        NODE_ASIN,
        NODE_ATAN,
        NODE_ATAN2,
        NODE_COS,
        NODE_COSH,
        NODE_FABS,
        NODE_FLOOR,
        NODE_SIN,
        NODE_SINH,
        NODE_SQRT,
        NODE_TAN,
        NODE_TANH
    };

    template<class REAL_T>
    struct ExpressionNode {
        typedef const REAL_T(*EvalNode)(ExpressionNode<REAL_T>*);
        typedef void (*AccumulateNode)(ExpressionNode<REAL_T>*, const REAL_T&, REAL_T);
        typedef void (*ForwardSweepNode)(ExpressionNode<REAL_T>*, const REAL_T&);
        typedef void (*DestroyNode)(ExpressionNode<REAL_T>*);

        Node_Operator op;

        union ADHolder {
            ExpressionNode<REAL_T>* ptr;
            VariableInfo<REAL_T>* derivative;
            REAL_T value;
        };
        ADHolder left;
        ADHolder right;
        REAL_T value;

        inline const REAL_T Evaluate() {
            value = evaluate_func[op](this);
            return value;
        }

        inline void Accumulate(const REAL_T& w, REAL_T coefficient = 1.0) {
            accumulate_func[op](this, w, coefficient);
        }

        inline void Destroy() {

        }

        ~ExpressionNode() {
            Destroy();
        }

        //evaluators

        static inline const REAL_T Add(ExpressionNode<REAL_T>* exp) {
            return exp->left.ptr->Evaluate() + exp->right.ptr->Evaluate();
        }

        static inline const REAL_T Subtract(ExpressionNode<REAL_T>* exp) {
            return exp->left.ptr->Evaluate() - exp->right.ptr->Evaluate();
        }

        static inline const REAL_T Multiply(ExpressionNode<REAL_T>* exp) {
            return exp->left.ptr->Evaluate() * exp->right.ptr->Evaluate();
        }

        static inline const REAL_T Divide(ExpressionNode<REAL_T>* exp) {
            return exp->left.ptr->Evaluate() / exp->right.ptr->Evaluate();
        }

        static inline const REAL_T Pow(ExpressionNode<REAL_T>* exp) {
            return std::pow(exp->left.ptr->Evaluate(), exp->right.ptr->Evaluate());
        }

        static inline const REAL_T Variable(ExpressionNode<REAL_T>* exp) {
            return exp->value;
        }

        static inline const REAL_T Scalar(ExpressionNode<REAL_T>* exp) {
            return exp->value;
        }

        static inline const REAL_T Log(ExpressionNode<REAL_T>* exp) {
            return std::log(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Log10(ExpressionNode<REAL_T>* exp) {
            return std::log10(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Exp(ExpressionNode<REAL_T>* exp) {
            return std::exp(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Acos(ExpressionNode<REAL_T>* exp) {
            return std::acos(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Asin(ExpressionNode<REAL_T>* exp) {
            return std::asin(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Atan(ExpressionNode<REAL_T>* exp) {
            return std::atan(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Atan2(ExpressionNode<REAL_T>* exp) {
            return std::atan2(exp->left.ptr->Evaluate(), exp->right.ptr->Evaluate());
        }

        static inline const REAL_T Cos(ExpressionNode<REAL_T>* exp) {
            return std::cos(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Cosh(ExpressionNode<REAL_T>* exp) {
            return std::cosh(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Fabs(ExpressionNode<REAL_T>* exp) {
            return std::fabs(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Floor(ExpressionNode<REAL_T>* exp) {
            return std::floor(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Sin(ExpressionNode<REAL_T>* exp) {
            return std::sin(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Sinh(ExpressionNode<REAL_T>* exp) {
            return std::sinh(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Sqrt(ExpressionNode<REAL_T>* exp) {
            return std::sqrt(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Tan(ExpressionNode<REAL_T>* exp) {
            return std::tan(exp->left.ptr->Evaluate());
        }

        static inline const REAL_T Tanh(ExpressionNode<REAL_T>* exp) {
            return std::tanh(exp->left.ptr->Evaluate());
        }

        //accumulators

        static inline void AccumulateAdd(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient);
            exp->right.ptr->Accumulate(w, coefficient);
        }

        static inline void AccumulateSubtract(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient);
            exp->right.ptr->Accumulate(w, static_cast<REAL_T> (-1.0) * coefficient);
        }

        static inline void AccumulateMultiply(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * exp->right.value);
            exp->right.ptr->Accumulate(w, coefficient * exp->left.value);
        }

        static inline void AccumulateDivide(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient / exp->right.value);
            exp->right.ptr->Accumulate(w, -1.0 * coefficient * exp->value / exp->right.value);
        }

        static inline void AccumulatePow(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * exp->right.value * pow(exp->left.value, exp->right.value - static_cast<REAL_T> (1.0)));
            exp->right.ptr->Accumulate(w, -1.0 * coefficient * exp->value * std::log(exp->left.value));
        }

        static inline void AccumulateVariable(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->right.derivative->dvalue += w*coefficient;
        }

        static inline void AccumulateScalar(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
        }

        static inline void AccumulateLog(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient / exp->value);
        }

        static inline void AccumulateLog10(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * REAL_T(1.0) / exp->value *
                    2.30258509299404590109361379290930926799774169921875);
        }

        static inline void AccumulateExp(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * exp->value);
        }

        static inline void AccumulateAcos(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient *
                    static_cast<REAL_T> (-1.0) /
                    std::pow((static_cast<REAL_T> (1.0) -
                    std::pow(exp->left.value, static_cast<REAL_T> (2.0))),
                    static_cast<REAL_T> (0.5)));
        }

        static inline void AccumulateAsin(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient *
                    static_cast<REAL_T> (1.0) /
                    std::pow((static_cast<REAL_T> (1.0) -
                    std::pow(exp->left.value,
                    static_cast<REAL_T> (2.0))), static_cast<REAL_T> (0.5)));
        }

        static inline void AccumulateAtan(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient *
                    static_cast<REAL_T> (1.0) / (exp->value * exp->value + static_cast<REAL_T> (1.0)));
        }

        static inline void AccumulateAtan2(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            std::cout << "AccumulateAtan2 not yey implemented!\n";
        }

        static inline void AccumulateCos(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * static_cast<REAL_T> (-1.0) * std::sin(exp->value));
        }

        static inline void AccumulateCosh(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * std::sinh(exp->value));
        }

        static inline void AccumulateFabs(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * exp->value / std::fabs(exp->value));
        }

        static inline void AccumulateFloor(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            std::cout << "AccumulateFloor not yey implemented!\n";
        }

        static inline void AccumulateSin(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * std::cos(exp->value));
        }

        static inline void AccumulateSinh(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * std::cosh(exp->value));
        }

        static inline void AccumulateSqrt(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
             exp->left.ptr->Accumulate(w, coefficient * (exp->value));
        }

        static inline void AccumulateTan(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            exp->left.ptr->Accumulate(w, coefficient * (static_cast<REAL_T> (1.0) / std::cos(exp->value)) * static_cast<REAL_T> (2.0));
        }

        static inline void AccumulateTanh(ExpressionNode<REAL_T>* exp, const REAL_T& w, REAL_T coefficient = 1.0) {
            REAL_T temp = std::exp(static_cast<REAL_T> (2.0) * exp->left.ptr->value);
            exp->left.ptr->Accumulate(w, coefficient  *
                    static_cast<REAL_T> (4.0) * temp /
                    (static_cast<REAL_T> (2.0 + 2) * temp + static_cast<REAL_T> (1.0)));
        }

        //destroyers

        static inline void DestroyAdd(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroySubtract(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyMultiply(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyDivide(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyPow(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyVariable(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyScalar(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyLog(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyLog10(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyExp(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyAcos(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyAsin(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyAtan(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyAtan2(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyCos(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyCosh(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyFabs(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyFloor(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroySin(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroySinh(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroySqrt(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyTan(ExpressionNode<REAL_T>* exp) {
        }

        static inline void DestroyTanh(ExpressionNode<REAL_T>* exp) {
        }

        static EvalNode evaluate_func[23];
        static AccumulateNode accumulate_func[23];
        static DestroyNode destroy_func[23];

        void* operator new(size_t size) {

            return malloc(size);
        }

        void operator delete(void* ptr) {
            free(ptr);
        }

    };

    template<typename REAL_T>
    typename ExpressionNode<REAL_T>::EvalNode ExpressionNode<REAL_T>::evaluate_func[23] = {&ExpressionNode<REAL_T>::Add,
        &ExpressionNode<REAL_T>::Subtract,
        &ExpressionNode<REAL_T>::Multiply,
        &ExpressionNode<REAL_T>::Divide,
        &ExpressionNode<REAL_T>::Pow,
        &ExpressionNode<REAL_T>::Variable,
        &ExpressionNode<REAL_T>::Scalar,
        &ExpressionNode<REAL_T>::Log,
        &ExpressionNode<REAL_T>::Log10,
        &ExpressionNode<REAL_T>::Exp,
        &ExpressionNode<REAL_T>::Acos,
        &ExpressionNode<REAL_T>::Asin,
        &ExpressionNode<REAL_T>::Atan,
        &ExpressionNode<REAL_T>::Atan2,
        &ExpressionNode<REAL_T>::Cos,
        &ExpressionNode<REAL_T>::Cosh,
        &ExpressionNode<REAL_T>::Fabs,
        &ExpressionNode<REAL_T>::Floor,
        &ExpressionNode<REAL_T>::Sin,
        &ExpressionNode<REAL_T>::Sinh,
        &ExpressionNode<REAL_T>::Sqrt,
        &ExpressionNode<REAL_T>::Tan,
        &ExpressionNode<REAL_T>::Tanh};

    template<typename REAL_T>
    typename ExpressionNode<REAL_T>::AccumulateNode ExpressionNode<REAL_T>::accumulate_func[23] = {&ExpressionNode<REAL_T>::AccumulateAdd,
        &ExpressionNode<REAL_T>::AccumulateSubtract,
        &ExpressionNode<REAL_T>::AccumulateMultiply,
        &ExpressionNode<REAL_T>::AccumulateDivide,
        &ExpressionNode<REAL_T>::AccumulatePow,
        &ExpressionNode<REAL_T>::AccumulateVariable,
        &ExpressionNode<REAL_T>::AccumulateScalar,
        &ExpressionNode<REAL_T>::AccumulateLog,
        &ExpressionNode<REAL_T>::AccumulateLog10,
        &ExpressionNode<REAL_T>::AccumulateExp,
        &ExpressionNode<REAL_T>::AccumulateAcos,
        &ExpressionNode<REAL_T>::AccumulateAsin,
        &ExpressionNode<REAL_T>::AccumulateAtan,
        &ExpressionNode<REAL_T>::AccumulateAtan2,
        &ExpressionNode<REAL_T>::AccumulateCos,
        &ExpressionNode<REAL_T>::AccumulateCosh,
        &ExpressionNode<REAL_T>::AccumulateFabs,
        &ExpressionNode<REAL_T>::AccumulateFloor,
        &ExpressionNode<REAL_T>::AccumulateSin,
        &ExpressionNode<REAL_T>::AccumulateSinh,
        &ExpressionNode<REAL_T>::AccumulateSqrt,
        &ExpressionNode<REAL_T>::AccumulateTan,
        &ExpressionNode<REAL_T>::AccumulateTanh};

    template<typename REAL_T>
    typename ExpressionNode<REAL_T>::DestroyNode ExpressionNode<REAL_T>::destroy_func[23] = {&ExpressionNode<REAL_T>::DestroyAdd,
        &ExpressionNode<REAL_T>::DestroySubtract,
        &ExpressionNode<REAL_T>::DestroyMultiply,
        &ExpressionNode<REAL_T>::DestroyDivide,
        &ExpressionNode<REAL_T>::DestroyPow,
        &ExpressionNode<REAL_T>::DestroyVariable,
        &ExpressionNode<REAL_T>::DestroyScalar,
        &ExpressionNode<REAL_T>::DestroyLog,
        &ExpressionNode<REAL_T>::DestroyLog10,
        &ExpressionNode<REAL_T>::DestroyExp,
        &ExpressionNode<REAL_T>::DestroyAcos,
        &ExpressionNode<REAL_T>::DestroyAsin,
        &ExpressionNode<REAL_T>::DestroyAtan,
        &ExpressionNode<REAL_T>::DestroyAtan2,
        &ExpressionNode<REAL_T>::DestroyCos,
        &ExpressionNode<REAL_T>::DestroyCosh,
        &ExpressionNode<REAL_T>::DestroyFabs,
        &ExpressionNode<REAL_T>::DestroyFloor,
        &ExpressionNode<REAL_T>::DestroySin,
        &ExpressionNode<REAL_T>::DestroySinh,
        &ExpressionNode<REAL_T>::DestroySqrt,
        &ExpressionNode<REAL_T>::DestroyTan,
        &ExpressionNode<REAL_T>::DestroyTanh};





}






#endif	/* EXPRESSIONNODE_HPP */

